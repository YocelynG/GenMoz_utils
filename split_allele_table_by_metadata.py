#!/usr/bin/env python3
"""
Split an allele table into one allele table per value of a metadata variable.

For every distinct value of --variable in the metadata (e.g. each province), the
allele-table rows whose SampleID belongs to a sample with that value are written to
  <outdir>/allele_table_<value>_filtered.txt

Notes:
  - Comma- or tab-separated files are detected automatically (any extension).
  - --variable year: if the metadata has no "year" column, the year is taken from
    the "date" column (first 4-digit number, e.g. 3/17/2022 -> 2022).
  - Allele-table samples absent from the metadata (e.g. lab controls) are not
    assigned to any group and are reported, not written.
  - Metadata samples with an empty value for --variable are skipped and reported.
  - Group values are made filename-safe: whitespace -> "_", other characters except
    letters, digits, "." and "-" are dropped 

Examples:
  python3 split_allele_table_by_metadata.py --input allele_table.txt --metadata metadata.txt --variable province
  python3 split_allele_table_by_metadata.py --input allele_table.txt --metadata metadata.txt --variable region --outdir allele_region
 
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import pandas as pd


def safe_name(value: str) -> str:
    name = re.sub(r"\s+", "_", value.strip())
    name = re.sub(r"[^\w.-]", "", name)
    return name.strip("._") or "NA"


def sniff_format(path: Path) -> tuple[str, int]:
    """Return (delimiter, quoting) from the header line: tab vs comma, quoted or not."""
    with open(path, newline="") as fh:
        header = fh.readline()
    sep = "\t" if header.count("\t") > header.count(",") else ","
    quoting = csv.QUOTE_ALL if header.lstrip().startswith('"') else csv.QUOTE_MINIMAL
    return sep, quoting


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, type=Path,
                   help="Allele table to split (comma- or tab-separated)")
    p.add_argument("--metadata", required=True, type=Path,
                   help="Metadata file (comma- or tab-separated)")
    p.add_argument("--variable", required=True,
                   help="Metadata variable to split by (e.g. province, district, us, year)")
    p.add_argument("--outdir", type=Path, default=None,
                   help="Output directory (default: allele_tables_by_<variable> next to --input)")
    p.add_argument("--sample-col", default="SampleID",
                   help="Sample ID column shared by both files (default: %(default)s)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    for path in (args.input, args.metadata):
        if not path.is_file():
            sys.exit(f"ERROR: file not found: {path}")

    meta_sep, _ = sniff_format(args.metadata)
    in_sep, in_quoting = sniff_format(args.input)
    meta = pd.read_csv(args.metadata, sep=meta_sep, dtype=str)
    # Read every field as text and keep it untouched so outputs match the input byte-for-byte.
    alleles = pd.read_csv(args.input, sep=in_sep, dtype=str, keep_default_na=False)

    for name, df in (("metadata", meta), ("allele table", alleles)):
        if args.sample_col not in df.columns:
            sys.exit(f"ERROR: column '{args.sample_col}' not found in {name}. Columns: {list(df.columns)}")

    if args.variable not in meta.columns:
        if args.variable == "year" and "date" in meta.columns:
            meta["year"] = meta["date"].str.extract(r"(\d{4})", expand=False)
            print("NOTE: no 'year' column in metadata; year taken from the 'date' column")
        else:
            sys.exit(f"ERROR: variable '{args.variable}' not found in metadata. Columns: {list(meta.columns)}")

    dup = meta[args.sample_col][meta[args.sample_col].duplicated()]
    if not dup.empty:
        sys.exit(f"ERROR: duplicated {args.sample_col} in metadata: {sorted(dup.unique())}")

    groups = meta[args.variable].str.strip()
    missing = meta.loc[groups.isna() | (groups == ""), args.sample_col]
    if not missing.empty:
        print(f"WARNING: {len(missing)} metadata samples have no '{args.variable}' value and are skipped: "
              f"{sorted(missing)}")
    sample_to_group = dict(zip(meta[args.sample_col], groups))

    # Make sure two different values don't collapse onto the same filename.
    file_names: dict[str, str] = {}
    for value in sorted(groups.dropna().loc[lambda s: s != ""].unique()):
        fname = f"allele_table_{safe_name(value)}_filtered.txt"
        if fname in file_names:
            sys.exit(f"ERROR: '{value}' and '{file_names[fname]}' both map to {fname}")
        file_names[fname] = value

    sample_group = alleles[args.sample_col].map(sample_to_group)
    unassigned = sorted(alleles.loc[sample_group.isna(), args.sample_col].unique())
    if unassigned:
        print(f"NOTE: {len(unassigned)} allele-table samples are not in the metadata and are excluded: "
              f"{', '.join(unassigned)}")
    no_alleles = sorted(set(meta[args.sample_col]) - set(alleles[args.sample_col]))
    if no_alleles:
        print(f"WARNING: {len(no_alleles)} metadata samples have no rows in the allele table: "
              f"{', '.join(no_alleles)}")

    outdir = args.outdir or args.input.resolve().parent / f"allele_tables_by_{safe_name(args.variable)}"
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\nSplitting by '{args.variable}' -> {outdir}")
    print(f"{'group':<45}{'samples':>8}{'rows':>10}  file")
    for fname, value in file_names.items():
        subset = alleles[sample_group == value]
        if subset.empty:
            print(f"{value:<45}{0:>8}{0:>10}  (no allele data, not written)")
            continue
        subset.to_csv(outdir / fname, sep=in_sep, index=False, quoting=in_quoting, lineterminator="\n")
        print(f"{value:<45}{subset[args.sample_col].nunique():>8}{len(subset):>10}  {fname}")

    written = sample_group.notna() & sample_group.ne("")
    print(f"\nTotal: {alleles.loc[written, args.sample_col].nunique()} samples, {int(written.sum())} rows "
          f"written across {len(file_names)} groups")


if __name__ == "__main__":
    main()
