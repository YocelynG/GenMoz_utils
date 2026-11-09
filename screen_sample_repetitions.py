#!/usr/bin/env python3
"""Find samples that are in more than one allele table (sequenced at different
moments), and optionally resolve them.

Input
  --alleles   Two or more allele tables (.csv, .tsv, .txt, optionally .gz), or
              folders: in a folder, every file whose name starts with "allele"
              (allele_table_run1.txt, allele_data.tsv.gz, ...) is read. All the
              tables are combined into one. Each needs a sample ID column
              (SampleID, sample_id, sample_name or sample), a locus column and a
              reads column. Column names are matched ignoring case.

Origin of each sample: the allele table it comes from. A sample found in more
than one allele table is a repeat.

Only samples are screened: IDs that standardise to a number, such as
"N2076817_1_S12_L001", "N2076817-1_S3_L001", "N2076817_1_S64__HFS23_NextSEQ01_271025",
"N2076817-1" and "2076817.1", which are all the same sample (2076817.1; a
"-0"/"_0" suffix means no decimal). Other IDs are controls (N3D7-10K,
NDd2_100k, ...): they are not screened and are left unchanged.

Screening table (--out)
  SampleID, Allele_Table, Raw_SampleID, N_Runs, Total_Reads, Total_Loci,
  Avg_Reads_per_Locus, Depth_Read_Loci, Best_Run

  Allele_Table         allele table the sample's rows come from
  Raw_SampleID         ID(s) as written in that table
  N_Runs               number of allele tables the sample is in
  Total_Reads          sum of reads over all alleles
  Total_Loci           loci with at least one read
  Avg_Reads_per_Locus  Total_Reads / Total_Loci
  Depth_Read_Loci      loci whose reads (summed over alleles) are > --min-depth
  Best_Run             "yes" for the table with most Depth_Read_Loci per sample
                       (ties broken by Total_Reads)

  Rows of the same sample are listed together, best table first.

What to do with samples in more than one allele table (use at most one; by
default they are only reported)
  --merge TRUE          one row per sample: reads are added up per locus (and
                        per asv in the allele table) across its tables.
                        Allele_Table and Raw_SampleID list all of them, best
                        first, separated by ";", and Best_Run holds the best table.
  --best_coverage TRUE  keep only the best table of each sample.
  --keep_dups TRUE      keep every table, adding _a, _b, ... to the IDs (a = best
                        table), right after _S*/_S*_L001 or else at the end:
                        N2076817-1_S3_L001_a, N2076817_1_S64_b__HFS23_NextSEQ01_271025,
                        2076817.1_a. The letter goes in the sample ID column and
                        in NIDA if the table has it. The screening table gets a
                        New_SampleID column.

--out-alleles FILE also writes one allele table, all input tables combined,
with that option applied. Rows of samples in only one table, and of controls,
are unchanged. An Allele_Table column is added (or kept, if the tables already
have one) with the name of the table each row comes from; merged rows list the
tables of the sample, best first.

Examples
  python3 screen_sample_repetitions.py --alleles all_runs/ \\
      --out screening.tsv --best_coverage TRUE --out-alleles allele_table_best_run.txt
  python3 screen_sample_repetitions.py \\
      --alleles allele_table_run1.txt allele_table_run2.txt --out screening.tsv

  Write the outputs outside the folder of the allele tables.
"""
import argparse
import gzip
import os
import re
import sys

import pandas as pd

# Column names are matched ignoring case (SampleID = sampleID, Reads = reads)
ID_NAMES = ["sampleid", "sample_id", "sample_name", "sample"]
# Standardised IDs of samples (2076817, 2076817.1); any other ID is a control
SAMPLE_ID_PATTERN = r"\d+(\.\d+)?"
MODES = ("merge", "best_coverage", "keep_dups")
# Files read from a folder given to --alleles
ALLELE_FILE_PATTERN = re.compile(r"^[Aa]llele.*\.(txt|tsv|csv)(\.gz)?$")
TABLE_COL = "Allele_Table"   # column naming the allele table of each row


def true_false(value):
    v = str(value).strip().lower()
    if v in ("true", "t", "yes", "y", "1"):
        return True
    if v in ("false", "f", "no", "n", "0"):
        return False
    raise argparse.ArgumentTypeError("use TRUE or FALSE")


def letter(rank):
    """0 -> a, 1 -> b, ..., 25 -> z, 26 -> aa."""
    s, n = "", rank + 1
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(97 + r) + s
    return s


def add_letter(raw, tag):
    """Put tag right after _S12 / _S3_L001 if the ID has it, otherwise at the end:
    N2076817_1_S64__HFS23_NextSEQ01_271025 -> N2076817_1_S64_a__HFS23_NextSEQ01_271025."""
    m = re.search(r"_S\d+(?:_L\d+)?(?=__|$)", raw)
    return raw[:m.end()] + tag + raw[m.end():] if m else raw + tag


def number_text(v):
    return str(int(v)) if float(v).is_integer() else str(float(v))


def separator(path):
    return "," if path.lower().removesuffix(".gz").endswith(".csv") else "\t"


def shortlist(values, n=20):
    values = list(values)
    return ", ".join(map(str, values[:n])) + (" ..." if len(values) > n else "")


def join_unique(values):
    """'a', 'b;a', 'c', '' -> 'a;b;c' (first appearance order, empty values left out)."""
    return ";".join(dict.fromkeys(x for v in values for x in v.split(";") if x))


def join_sorted(values):
    return ";".join(sorted(set(values)))


def summarise(per_locus, keys, min_depth):
    """Reads, loci and deep loci for each group of per-locus read totals."""
    s = per_locus.groupby(keys).agg(
        Total_Reads=("reads", "sum"),
        Total_Loci=("reads", lambda v: int((v > 0).sum())),
        Depth_Read_Loci=("reads", lambda v: int((v > min_depth).sum()))).reset_index()
    s["Avg_Reads_per_Locus"] = (s.Total_Reads / s.Total_Loci.where(s.Total_Loci > 0)).round(2)
    return s


def normalise_id(raw):
    s = str(raw).strip()
    s = re.sub(r"__.*$", "", s)           # __MULE_NextSEQ03_300725
    s = re.sub(r"_S\d+(_L\d+)?$", "", s)  # _S12_L001
    m = re.match(r"^N?(\d+)(?:[-_.](\d+))?$", s)
    if m:
        base, dec = m.groups()
        return base if dec in (None, "0") else f"{base}.{dec}"
    return s.replace("-", "_")            # controls, e.g. N3D7-10K -> N3D7_10K


def read_table(path):
    if path.lower().endswith((".xlsx", ".xls")):
        sys.exit(f"{path}: save it as .csv or tab-separated .txt first.")
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="latin-1") as fh:
        first = fh.readline()
    sep = "\t" if "\t" in first else ","
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, encoding="latin-1")


def find_col(df, wanted, names, label):
    if wanted:
        if wanted not in df.columns:
            sys.exit(f"Column '{wanted}' not found for {label}. Columns: {list(df.columns)}")
        return wanted
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n in lower:
            return lower[n]
    return None


def allele_files(paths):
    """Allele table files from the files and folders given. In a folder, every file whose name
    starts with "allele" and ends in .txt, .tsv or .csv (optionally .gz) is read."""
    files = {}
    for path in paths:
        if os.path.isdir(path):
            folder = path.rstrip("/") or path
            found = [os.path.join(folder, f) for f in sorted(os.listdir(folder))
                     if ALLELE_FILE_PATTERN.match(f) and not os.path.isdir(os.path.join(folder, f))]
            if not found:
                sys.exit(f"{path}: no allele tables in this folder (files named allele*.txt, .tsv or .csv, "
                         "optionally .gz).")
        elif os.path.isfile(path):
            found = [path]
        else:
            sys.exit(f"Allele table or folder not found: {path}")
        for f in found:
            files.setdefault(os.path.abspath(f), f)
    return list(files.values())


def table_names(files):
    """File name of each allele table, or its path when two tables have the same file name."""
    base = [os.path.basename(f) for f in files]
    return [f if base.count(b) > 1 else b for f, b in zip(files, base)]


def read_alleles(files):
    """All allele tables as one, plus the name of its Allele_Table column.

    Column names that differ only in case (sampleID, SampleID) get the first table's spelling.
    Allele_Table holds the name of each row's table (kept where the tables already have one).
    """
    tables, spelling = [], {}
    for path in files:
        df = read_table(path)
        df = df.rename(columns={c: spelling[c.lower()] for c in df.columns
                                if c.lower() in spelling and spelling[c.lower()] not in df.columns})
        for c in df.columns:
            spelling.setdefault(c.lower(), c)
        tables.append(df)
    table_col = spelling.get(TABLE_COL.lower(), TABLE_COL)
    for df, name in zip(tables, table_names(files)):
        if table_col in df.columns:
            df[table_col] = df[table_col].where(df[table_col].str.strip() != "", name)
        else:
            df[table_col] = name
    # Columns missing from a table are left empty
    return pd.concat(tables, ignore_index=True).fillna(""), table_col


def resolve_alleles(al, c, links, repeated, ranks, mode):
    """Allele table with each sample in several tables merged, cut to its best table, or relabelled.

    c       column names in al: id, locus, reads, asv, run, nida, allele, norm, n_alleles, table
    links   row (index in al), SampleID and Allele_Table of every allele row of a sample
    repeated  rows left out when adding up reads (same sequence listed twice)
    ranks   SampleID, Allele_Table and Rank (0 = best) of the samples in several tables
    """
    dup = links.merge(ranks, on=["SampleID", "Allele_Table"])
    if dup.empty:   # no sample in more than one table
        return al
    if mode == "best_coverage":
        return al.drop(index=dup.row[dup.Rank > 0])

    if mode == "keep_dups":
        out = al.copy()
        tag = pd.Series(("_" + dup.Rank.map(letter)).to_numpy(), index=dup.row.to_numpy())
        for col in dict.fromkeys(x for x in (c["id"], c["nida"]) if x):
            vals = out.loc[tag.index, col]
            out.loc[tag.index, col] = [v if v == "" else add_letter(v, t) for v, t in zip(vals, tag)]
        return out

    # merge: one row per sample, locus and sequence; the other columns come
    # from the best table that has that sequence
    locus, reads = c["locus"], c["reads"]
    keys = ["_sid", locus, c["asv"]]
    use = dup[~dup.row.isin(repeated)]
    part = al.loc[use.row].assign(_sid=use.SampleID.to_numpy(), _rank=use.Rank.to_numpy())
    part["_row"] = part.index
    part["_reads"] = pd.to_numeric(part[reads], errors="coerce").fillna(0)
    part = part.sort_values(["_sid", "_rank", "_row"])
    part["_total"] = part.groupby(keys)._reads.transform("sum")
    # Keep each merged sample together, at the place of its first row
    part["_p"] = part.groupby("_sid")._row.transform("min")
    part["_l"] = part.groupby(["_sid", locus])._row.transform("min")
    m = part.drop_duplicates(keys).copy()

    m[reads] = m._total.map(number_text)
    at_locus = m.groupby(["_sid", locus])._total.transform("sum")
    if c["norm"]:
        m[c["norm"]] = (m._total / at_locus.where(at_locus > 0)).map(lambda v: "" if pd.isna(v) else str(v))
    if c["n_alleles"]:
        m[c["n_alleles"]] = m.groupby(["_sid", locus])._total.transform("size").astype(str)
    m[c["id"]] = m._sid.map(part[part._rank == 0].groupby("_sid")[c["id"]].first())
    # Tables (and runs) of the sample, best first
    m[c["table"]] = m._sid.map(part.groupby("_sid", sort=False)[c["table"]].agg(join_unique))
    if c["run"]:
        m[c["run"]] = m._sid.map(part[c["run"]].str.strip().groupby(part._sid, sort=False).agg(join_unique))
    if c["allele"]:
        clash = m.duplicated(["_sid", locus, c["allele"]], keep=False)
        if clash.any():
            print(f"Note: {m[clash].groupby(['_sid', locus, c['allele']]).ngroups} allele names are used for two "
                  "different sequences in the same merged sample, because each table numbers its alleles "
                  "separately. Match alleles by asv, not by allele name.")

    rest = al.drop(index=dup.row)
    rest = rest.assign(_p=rest.index, _l=rest.index, _row=rest.index)
    out = pd.concat([rest, m]).sort_values(["_p", "_l", "_row"], kind="stable")
    return out[list(al.columns)]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--alleles", nargs="+", required=True, metavar="PATH",
                   help="allele table(s), or folder(s) whose allele* files are all read; combined into one table")
    p.add_argument("--out", required=True, help="screening table (.csv for commas, anything else is tab-separated)")
    p.add_argument("--min-depth", type=float, default=100, help="reads a locus must exceed for Depth_Read_Loci (default 100)")
    p.add_argument("--all", action="store_true", help="also list samples in only one allele table")
    p.add_argument("--merge", type=true_false, default=False, metavar="TRUE|FALSE",
                   help="TRUE: one row per sample, adding up its tables (default FALSE)")
    p.add_argument("--best_coverage", "--best-coverage", dest="best_coverage", type=true_false, default=False,
                   metavar="TRUE|FALSE", help="TRUE: keep only the best table of each sample (default FALSE)")
    p.add_argument("--keep_dups", "--keep-dups", dest="keep_dups", type=true_false, default=False,
                   metavar="TRUE|FALSE", help="TRUE: keep every table, adding _a, _b, ... to the IDs (default FALSE)")
    p.add_argument("--out-alleles", "--out_alleles", dest="out_alleles", metavar="FILE",
                   help="also write the allele table with --merge, --best_coverage or --keep_dups applied")
    p.add_argument("--alleles-id-col", help="sample ID column in --alleles (default: sampleID or sample_name)")
    a = p.parse_args()

    modes = [m for m in MODES if getattr(a, m)]
    if len(modes) > 1:
        p.error("use only one of --merge, --best_coverage and --keep_dups")
    mode = modes[0] if modes else None
    if a.out_alleles and not mode:
        p.error("--out-alleles needs --merge TRUE, --best_coverage TRUE or --keep_dups TRUE")
    outputs = [os.path.abspath(x) for x in (a.out, a.out_alleles) if x]
    files = allele_files(a.alleles)
    if {os.path.abspath(f) for f in files} & set(outputs) or len(set(outputs)) < len(outputs):
        p.error("--out and --out-alleles must be different from each other and from the allele tables "
                "(including the allele* files in a folder given to --alleles)")

    print(f"Allele tables read ({len(files)}): {shortlist(files)}")
    al, table_col = read_alleles(files)
    sid = find_col(al, a.alleles_id_col, ID_NAMES, "--alleles IDs")
    locus = find_col(al, None, ["locus"], "--alleles locus")
    reads = find_col(al, None, ["reads"], "--alleles reads")
    run = find_col(al, None, ["run"], "--alleles Run")
    asv = find_col(al, None, ["asv"], "--alleles asv")
    if not (sid and locus and reads):
        sys.exit(f"--alleles needs sample ID, locus and reads columns. Columns: {list(al.columns)}")
    if a.out_alleles and mode == "merge" and not asv:
        sys.exit("--merge TRUE with --out-alleles needs an asv column in --alleles to match alleles across tables "
                 "(allele names are numbered separately in each table).")

    # Origin of each row: the allele table it comes from
    n_tables = al[table_col].nunique()
    data = pd.DataFrame({"row": al.index, "raw": al[sid].str.strip(), "Allele_Table": al[table_col],
                         "run": al[run].str.strip() if run else "", "locus": al[locus],
                         "reads": pd.to_numeric(al[reads], errors="coerce").fillna(0)})
    if asv:
        data["asv"] = al[asv]
    data["SampleID"] = data.raw.map(normalise_id)
    # Only samples are screened; controls (N3D7-10K, NDd2_100k, ...) are left unchanged
    is_sample = data.SampleID.str.fullmatch(SAMPLE_ID_PATTERN)
    controls = sorted(set(data.SampleID[~is_sample]))
    data = data[is_sample]

    links = data[["row", "SampleID", "Allele_Table"]]
    repeated_rows = set()
    if asv:
        # The same sequence listed twice under two allele names would be counted twice
        repeated = data.duplicated(["raw", "Allele_Table", "run", "locus", "asv"])
        if repeated.any():
            print(f"Note: ignored {repeated.sum()} repeated sequence rows (same sample, allele table, locus and asv).")
        repeated_rows = set(data.row[repeated])
        data = data[~repeated]

    per_locus = data.groupby(["SampleID", "Allele_Table", "locus"], as_index=False).reads.sum()
    stats = summarise(per_locus, ["SampleID", "Allele_Table"], a.min_depth)
    raw_ids = data.groupby(["SampleID", "Allele_Table"]).raw.agg(join_sorted).rename("Raw_SampleID")
    stats = stats.merge(raw_ids.reset_index(), on=["SampleID", "Allele_Table"])
    # Samples sequenced twice inside one table (two sequencing IDs or runs): their reads are added up
    inside = data.drop_duplicates(["SampleID", "Allele_Table", "raw", "run"]).groupby(["SampleID", "Allele_Table"]).size()
    n_inside = inside[inside > 1].index.get_level_values("SampleID").nunique()
    stats["N_Runs"] = stats.groupby("SampleID").Allele_Table.transform("nunique")
    stats = stats.sort_values(["SampleID", "Depth_Read_Loci", "Total_Reads", "Allele_Table"],
                              ascending=[True, False, False, True])
    stats["Rank"] = stats.groupby("SampleID").cumcount()
    stats["Best_Run"] = (stats.Rank == 0).map({True: "yes", False: "no"})
    ranks = stats.loc[stats.N_Runs > 1, ["SampleID", "Allele_Table", "Rank"]]
    if not a.all:
        stats = stats[stats.N_Runs > 1]

    cols = ["SampleID", "Allele_Table", "Raw_SampleID", "N_Runs", "Total_Reads", "Total_Loci",
            "Avg_Reads_per_Locus", "Depth_Read_Loci", "Best_Run"]
    if mode == "merge":
        # Add up each sample's reads per locus over its tables, then recount
        pooled = per_locus[per_locus.SampleID.isin(stats.SampleID)]
        pooled = pooled.groupby(["SampleID", "locus"], as_index=False).reads.sum()
        merged = stats.groupby("SampleID", sort=False).agg(
            Allele_Table=("Allele_Table", join_unique), Raw_SampleID=("Raw_SampleID", ";".join),
            N_Runs=("Allele_Table", "size"), Best_Run=("Allele_Table", "first")).reset_index()
        stats = summarise(pooled, ["SampleID"], a.min_depth).merge(merged, on="SampleID").sort_values("SampleID")
    elif mode == "best_coverage":
        stats = stats[stats.Rank == 0]
    elif mode == "keep_dups":
        stats = stats.assign(New_SampleID=stats.SampleID + ("_" + stats.Rank.map(letter)).where(stats.N_Runs > 1, ""))
        cols.append("New_SampleID")
    stats["Total_Reads"] = stats.Total_Reads.round().astype("int64")
    # Files are read as latin-1 and written as latin-1, so every byte is kept as it was
    stats[cols].to_csv(a.out, sep=separator(a.out), index=False, encoding="latin-1")

    print(f"Allele tables: {n_tables} | samples: {links.SampleID.nunique()} | "
          f"controls (not screened, left unchanged): {len(controls)}")
    if n_tables == 1:
        print("Note: only one allele table, so no sample can be in more than one. Give two or more tables, "
              "or a folder with them.")
    if controls:
        print(f"Controls: {shortlist(controls)}")
    shape = {"merge": "merged into one row each",
             "best_coverage": "only the best allele table of each kept",
             "keep_dups": "all kept, IDs ending in _a, _b, ... (a = best)"}.get(mode, "listed in the screening table")
    print(f"Samples in more than one allele table: {ranks.SampleID.nunique()} ({shape})")
    if n_inside:
        print(f"Note: {n_inside} samples are more than once in the same allele table (two sequencing IDs or runs); "
              "their reads are added up as one (see Raw_SampleID).")
    print(f"Depth_Read_Loci = loci with more than {a.min_depth:g} reads")
    print(f"Written: {a.out}")

    if not a.out_alleles:
        if mode:
            print("Note: only the screening table was written; add --out-alleles FILE to also write the allele table.")
        return
    out = resolve_alleles(al, {"id": sid, "locus": locus, "reads": reads, "asv": asv, "run": run,
                               "nida": find_col(al, None, ["nida"], ""),
                               "allele": find_col(al, None, ["allele"], ""),
                               "norm": find_col(al, None, ["norm.reads.locus"], ""),
                               "n_alleles": find_col(al, None, ["n.alleles"], ""), "table": table_col},
                          links, repeated_rows, ranks, mode)
    out.to_csv(a.out_alleles, sep=separator(a.out_alleles), index=False, encoding="latin-1")
    print(f"Allele table rows: {len(al)} in, {len(out)} out")
    print(f"Written: {a.out_alleles}")


if __name__ == "__main__":
    main()
