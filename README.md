# Genomic Surveillance utils

Command-line scripts for processing malaria amplicon allele tables: split samples into metadata-defined groups, estimate read-based **complexity of infection (COI)**, calculate allele frequencies and heterozygosity **in populations** using [MOIRE](https://github.com/EPPIcenter/moire), and **between-host genetic relatedness (IBD)** using [Dcifer](https://github.com/EPPIcenter/dcifer) through [PGEcore](https://github.com/PlasmoGenEpi/PGEcore).

## Scripts

| Script | Purpose | Input | Main output |
| --- | --- | --- | --- |
| `screen_sample_repetitions.py` | Finds samples sequenced in more than one run and keep only the best sample among the duplicates | One or more allele tables (or a directory of `allele*` files) + manifest | Screening table and one combined allele table |
| `split_allele_table_by_metadata.py` | Divide a single allele table by a metadata variable (e.g., province, district or year) | Allele table + sample metadata | `allele_table_<group>_filtered.txt` per group |
| `run_coi.R` | Estimate read-based COI from the number of alleles required to cover a selected fraction of reads | One allele table or a directory | `<input_stem>_COI_Table.txt` |
| `run_dcifer.R` | Estimate pairwise between-host relatedness with Dcifer | Allele table(s); optional COI, allele frequencies and/or metadata | `<input_stem>_Dcifer_IBD_Results.txt` and reproducible intermediate files |


## Requirements and installation

Use a terminal with **Python 3** and **R** (accessible as `python3` and `Rscript`). For current `dcifer` releases, **R 4.4+** is recommended/required by the package; check the version requirements of the version you install


### Prerequisites: `screen_sample_repetitions.py`

- Python 3.9 or newer
- pandas (`pip install pandas` or `conda install pandas`)


### Prerequisites: `split_allele_table_by_metadata.py`

- Python 3.9+
- Python package: `pandas`

### Prerequisites:  `run_coi.R`

- R
- R packages: `optparse`, `dplyr`, `tibble` (`tibble` is used explicitly in the code)

### Prerequisites: `run_dcifer.R` and `run_MOIRE.R`

- R (compatible with your installed `dcifer` version)
- R packages: `optparse`, `dplyr`, `tibble`, `PGEcore`, `dcifer`, `moire`, `foreach`, `iterators`, `doParallel`, `parallelly`
- Multiple CPU workers if using `--threads` greater than 1


`PGEcore` does **not** automatically install the specialised `dcifer` package; install both. See [PGEcore installation instructions](https://plasmogenepi.github.io/PGEcore/articles/getting-started.html).

Check availability before running:

```bash
Rscript -e 'stopifnot(all(vapply(c("optparse", "dplyr", "tibble", "PGEcore", "dcifer", "foreach", "iterators", "doParallel", "parallelly"), requireNamespace, logical(1), quietly=TRUE)))'
```

## Input file formats

All three scripts accept **comma- or tab-separated** allele tables (including `.txt`, `.tsv` and `.csv`). They detect delimiters from the header instead of relying on the extension

### Allele tables generated with Mad4HaTter 

Example (tab-separated):

```text
SampleID	Locus	PseudoCIGAR	Reads
S01	LocusA	.	80
S01	LocusA	2A>G	20
S01	LocusB	.	100
S02	LocusA	.	100


### Metadata

For the splitting script, metadata must contain the same sample-ID column as the allele table and a column named by `--variable`. For Dcifer population mode, metadata must contain `SampleID` (or its normalized equivalent) and the specified `--pop_col` (default: `district`).

Example:

```text
SampleID,province,district,date
S01,Province_A,District_1,2022-03-17
S02,Province_B,District_2,2023-05-10
```

Metadata should have one row per sample.

## Usage

### 1. `screen_sample_repetitions.py`

Finds samples sequenced in more than one run, reports their depth in each run, and can write an allele table with those samples merged, cut to their best run, or renamed. It reads one allele table, several, or a directory of them; several tables are combined and checked together.

Column names are matched ignoring case, so `SampleID` = `sampleID` and `Reads` = `reads`.

**Inputs**

**`--alleles`:** one or more allele tables (`.txt`, `.tsv`, `.csv`, optionally `.gz`), or directories. In a directory, every file whose name starts with `allele` is read (e.g. `allele_table_run1.txt`, `allele_data.tsv.gz`). All the tables are combined into one.
- **Required columns:** the sample ID (`SampleID`, `sample_id`, `sample_name` or `sample`), `locus` and `reads`.
- **Optional columns:** `Run`, `asv`, `NIDA`, `norm.reads.locus` and `n.alleles`.
- **Runs:** the `Run` column tells which run each row comes from. Without it, each allele table counts as one run; with a single allele table, each sequencing ID (e.g. `N2076817-1_S3_L001`) counts as its own run. If only some tables have a `Run` column, the others get their file name in it.
- **Outputs:** write them outside the directory of the allele tables, or a later run will read an old output as one more table.

**`--manifest`:** one or more pipeline manifests with the sample ID (`SampleID`, `sample_id`, `sample_name` or `sample`) and `SampleType` columns. Only IDs with `SampleType = sample` are screened; controls and IDs not in a manifest are left unchanged.

```
SampleID,SampleType
2076817.1,sample
N3D7-10K,positive
```

## Options

| Option | Description |
|---|---|
| `--alleles PATH [PATH ...]` | Allele table(s), or directories of `allele*` files (required) |
| `--manifest FILE [FILE ...]` | Manifest(s) with `SampleID` and `SampleType` (required) |
| `--out FILE` | Screening table (required) |
| `--out-alleles FILE` | Also write one allele table, all inputs combined, with one of the three options below applied |
| `--best_coverage TRUE` | Keep only the best run of each sample |
| `--keep_dups TRUE` | Keep every run and add `_a`, `_b`, ... to the IDs (`a` = best run) |
| `--merge TRUE` | Combine each sample's runs into one, adding up the reads |
| `--min-depth N` | Reads a locus must exceed to count as covered (default 100) |
| `--all` | Also list samples in only one run in the screening table |

- **One option at a time:** use only one of `--best_coverage`, `--keep_dups` and `--merge`.
- **Best run:** the run with the most loci above `--min-depth` reads
- **`--keep_dups`:** the letter goes after `_S*_L001` in `sampleID` and at the end of `NIDA`
- **`--merge`:** needs an `asv` column, because allele names are numbered separately in each run. Downstream, match alleles by `asv`.

```bash
python3 screen_sample_repetitions.py --alleles allele_table.txt --manifest manifest.csv \
    --out screening.tsv 
```

**With options for duplicate samples**, --best_coverage TRUE: keep the sample with the best coverage

```bash
python3 screen_sample_repetitions.py --alleles allele_table.txt --manifest manifest.csv \
    --out screening.tsv --best_coverage TRUE --out-alleles allele_table_best_run.txt
```

**Several allele tables**, listed or as a directory; the output allele table is a single file

```bash
python3 screen_sample_repetitions.py --alleles allele_table_run1.txt allele_table_run2.txt \
    --manifest manifest.csv --out screening.tsv --best_coverage TRUE --out-alleles allele_table_best_run.txt

python3 screen_sample_repetitions.py --alleles allele_tables/ --manifest manifest.csv \
    --out screening.tsv --best_coverage TRUE --out-alleles allele_table_best_run.txt
```

## Outputs

**Screening table (`--out`):** one row per sample and run, best run first.

| Column | Meaning |
|---|---|
| `SampleID` | Standardised ID |
| `Run` | Sequencing run (without a `Run` column: the allele table, or the sequencing ID if there is only one table) |
| `Allele_Table` | Allele table(s) the sample's rows in that run come from |
| `Raw_SampleID` | ID as written in the allele table |
| `N_Runs` | Number of runs the sample is in |
| `Total_Reads` | Reads over all loci |
| `Total_Loci` | Loci with at least one read |
| `Avg_Reads_per_Locus` | `Total_Reads` / `Total_Loci` |
| `Depth_Read_Loci` | Loci with more than `--min-depth` reads |
| `Best_Run` | `yes` / `no` (with `--merge`: name of the best run) |
| `New_SampleID` | Only with `--keep_dups`: the ID with its letter |

With `--best_coverage`, only the best row of each sample is kept. With `--merge`, each sample has one row, and its runs and allele tables are listed best first, separated by `;`.

**Allele table (`--out-alleles`):** Remove duplicate samples, keeping only one or renaming duplicates. It is a single file with all the input tables combined, with the same columns as the input plus `Allele_Table`: the file name of the table each row comes from (its path if two tables have the same file name).

| Option | Change in the allele table |
|---|---|
| `--best_coverage` | Sample with the best read depth coverage is retained |
| `--keep_dups` | All duplicate samples are kept, and IDs get `_a`, `_b`, ... |
| `--merge` | Reads are summed and `norm.reads.locus` and `n.alleles` recalculated. `Run` and `Allele_Table` list all runs and tables, best first. |


### 2. Split allele tables by metadata

| Argument | Required | Description |
| --- | --- | --- |
| `-i`, `--input` | Yes | Allele file or directory of `allele*` files |
| `-o`, `--output` | No | Output directory; defaults to input file's directory or the input directory |
| `-t`, `--threshold` | No | Cumulative proportion of locus reads to cover, in `(0, 1]`; default `0.98` |
| `--count_reference` | No | Treat `PseudoCIGAR` value `.` as a reference allele; otherwise exclude it |
| `-h`, `--help` | No | Show help |

**How COI is calculated:** within each sample/locus, reads for identical alleles are summed and alleles are ranked from most to least abundant. Alleles are retained until the cumulative read fraction reaches the threshold, **including the allele that crosses it**. The number retained is calculated for each locus; the sample's **COI is the maximum across loci**. The mean of the retained allele counts across represented loci is also reported.

**Output:** `<input_stem>_COI_Table.txt`, a tab-separated table containing:

| Column | Meaning |
| --- | --- |
| `SampleID` | Sample identifier |
| `COI` | Maximum number of retained alleles among loci |
| `avg_alleles_per_locus` | Mean number of retained alleles per represented locus |

**Important:** Without `--count_reference`, `.` is excluded as a missing allele. In some allele-table conventions `.` denotes the *reference* allele; specify the flag if that applies to your data. Samples with no alleles remaining are omitted from the COI output, with a warning.


```bash
python3 split_allele_table_by_metadata.py \
  --input allele_table.txt \
  --metadata metadata.csv \
  --variable province \
  --outdir allele_tables_by_province
```

Alternatively, split by year (derived from `date` when needed):

```bash
python3 split_allele_table_by_metadata.py \
  --input allele_table.txt \
  --metadata metadata.csv \
  --variable year \
  --outdir allele_tables_by_year
```

**Running in the workstation**:

```bash
python3 split_allele_table_by_metadata.py --input allele_table.txt --metadata metadata.csv --variable province --outdir allele_tables_by_province 
```

### 3. Estimate complexity of infection (COI)

| Argument | Required | Description |
| --- | --- | --- |
| `-i`, `--input` | Yes | Allele file or directory of `allele*` files |
| `-o`, `--output` | No | Output directory; defaults to input file's directory or the input directory |
| `-t`, `--threshold` | No | Cumulative proportion of locus reads to cover, in `(0, 1]`; default `0.98` |
| `--count_reference` | No | Treat `PseudoCIGAR` value `.` as a reference allele; otherwise exclude it |
| `-h`, `--help` | No | Show help |

**How COI is calculated:** within each sample/locus, reads for identical alleles are summed and alleles are ranked from most to least abundant. Alleles are retained until the cumulative read fraction reaches the threshold, **including the allele that crosses it**. The number retained is calculated for each locus; the sample's **COI is the maximum across loci**. The mean of the retained allele counts across represented loci is also reported.

**Output:** `<input_stem>_COI_Table.txt`, a tab-separated table containing:

| Column | Meaning |
| --- | --- |
| `SampleID` | Sample identifier |
| `COI` | Maximum number of retained alleles among loci |
| `avg_alleles_per_locus` | Mean number of retained alleles per represented locus |

**Important:** Without `--count_reference`, `.` is excluded as a missing allele. In some allele-table conventions `.` denotes the *reference* allele; specify the flag if that applies to your data. Samples with no alleles remaining are omitted from the COI output, with a warning.


```bash
Rscript run_coi.R \
  --input allele_tables_by_province \
  --output coi_output \
  --threshold 0.98 \
```

**Running in the workstation** (default: 0.98):

```bash
Rscript run_coi.R --input allele_table.txt --output Coi_output --threshold 0.98
```

### 4. Estimate allele frequencies and Heterozygosis using MOIRE

Run MOIRE analyses on all allele_*.txt files in a directory.

Input files may be tab- or comma-separated, with either header format:
"SampleID","Locus","PseudoCIGAR","Reads","WSAF"

**Use the COI estimates from `run_MOIRE.R`:**

```
Rscript run_MOIRE.R [input_directory] [output_directory] [threads]
```

**Running in the workstation**:

```bash
nohup Rscript run_MOIRE.R [input_directory] [output_directory] 10 & 
```

### 4. Estimate pairwise relatedness using Dcifer

Run Dcifer to estimate allele frequencies from the data and  also calculates a naive COI, or use allele frequencies and COI previously calculated  

| Argument | Required | Description |
| --- | --- | --- |
| `-i`, `--input` | Yes | Allele table or directory of `allele*` files |
| `-o`, `--output` | Yes | Directory to create/write results |
| `-a`, `--allele_freq` | No | Matching MOIRE `*_AF_Results.txt` file or containing directory; default is Dcifer-estimated frequencies |
| `-c`, `--coi` | No | Matching MOIRE `*_COI_Results.txt` or read-based `*_COI_Table.txt` file, or containing directory; default is naive COI |
| `-m`, `--metadata` | No | Metadata for population-aware estimation |
| `-p`, `--pop_col` | No | Population column in metadata; default `district` |
| `-t`, `--threads` | No | Number of parallel workers; default `4` |
| `-h`, `--help` | No | Show help |

**Important rules and defaults:**

- `--allele_freq` **cannot** be combined with `--metadata` because the PGEcore wrapper calculates allele frequencies by population in metadata mode.
- Without `--coi`, naive COI is the **second-highest** number of distinct alleles observed across the sample's loci (or the only count if there is one locus). This is **different** from the read-threshold COI in `run_coi.R`.
- With `--coi`, a MOIRE file must contain `sample_id` and `post_coi_med`, or a `run_coi.R` file must contain `SampleID` and `COI`. MOIRE values are rounded to integer COI (minimum 1).
- MOIRE frequency files must contain `locus`, `allele`, and `post_allele_freqs_mean`, and must correspond to the alleles in the input table. Frequencies are renormalized within loci. Monomorphic loci absent from MOIRE frequencies are omitted; mismatched polymorphic loci cause an error.
- Dcifer collapses duplicate sample–locus–allele entries and removes masked segments from `PseudoCIGAR` values. Missing/empty allele strings become `.` in Dcifer input.
- Samples without a valid population assignment are omitted when `--metadata` is used.
- If analyzing a directory with supplied COI/AF results, provide a **directory** of matching result files, not one result file. Input stems must match: e.g. `allele_table_Province_A_filtered.txt` → `allele_table_Province_A_filtered_COI_Table.txt`.
- The script uses the fixed test settings `rnull = 0` and `alpha = 0.05`.

**Output:** `<input_stem>_Dcifer_IBD_Results.txt`, including:

| Column | Meaning |
| --- | --- |
| `specimen_name_a`, `specimen_name_b` | Pair of samples |
| `btwn_host_rel` | Estimated between-host relatedness |
| `p_value` | Relatedness test p-value |
| `CI_lower`, `CI_upper` | Confidence interval bounds |

The script also retains the exact files supplied to PGEcore in `<output>/dcifer_inputs/`:

```text
dcifer_output/
├── allele_table_Province_A_filtered_Dcifer_IBD_Results.txt
└── dcifer_inputs/
    ├── allele_table_Province_A_filtered_allele_table.tsv
    ├── allele_table_Province_A_filtered_coi.tsv
    ├── allele_table_Province_A_filtered_allele_freq.tsv  # if MOIRE AF supplied
    └── allele_table_Province_A_filtered_metadata.tsv     # if metadata supplied
```

**Use the COI estimates from `run_coi.R`:**


```bash
Rscript run_dcifer.R \
  --input allele_tables_by_province \
  --coi coi_output \
  --output dcifer_output \
  --threads 4
```

**Use MOIRE allele frequencies and COI** (requires matching files previously generated by `run_moire.R`, which is *not* included in this set of scripts):

```bash
Rscript run_dcifer.R \
  --input allele_tables_by_province \
  --allele_freq Moire_output \
  --coi Moire_output \
  --output dcifer_output \
  --threads 4
```

**Run population-aware analysis** using metadata. This estimates population-specific frequencies and considers within- and between-population pairs:

```bash
Rscript run_dcifer.R \
  --input allele_table.txt \
  --metadata metadata.csv \
  --pop_col district \
  --output dcifer_by_district \
  --threads 4
```

**Running in the workstation** (in the background):

```bash
nohup Rscript run_dcifer.R --input allele_tables_by_district --output Dcifer_output   --allele_freq Moire_output --coi COI_tables/ --threads 4 &
```

## Example end-to-end workflow

Assuming `allele_table.txt` and `metadata.csv` are available in the current directory:

```bash
# 1. Create per-province allele tables
python3 split_allele_table_by_metadata.py \
  --input allele_table.txt \
  --metadata metadata.csv \
  --variable province \
  --outdir allele_tables_by_province

# 2. Estimate COI for each province (count '.' as reference, if appropriate)
Rscript run_coi.R \
  --input allele_tables_by_province \
  --output coi_output \
  --count_reference

# 3. Estimate relatedness independently within each province
Rscript run_dcifer.R \
  --input allele_tables_by_province \
  --coi coi_output \
  --output dcifer_output \
  --threads 4
```

## Troubleshooting

- **Missing column:** Check the headers. The Python script requires an exact sample-column match; the R scripts accept normalized aliases.
- **No allele files found:** For directory mode, filenames must start with `allele` and have one of the supported extensions.
- **Missing COI or allele-frequency file:** Check that the result directory contains one corresponding file for **each** allele table, with the same basename before its output suffix.
- **Missing COI for some samples:** `run_coi.R` can omit samples with no retained alleles, while `run_dcifer.R --coi` requires COI for every included sample.
- **Reference allele interpretation:** Review how your pipeline encodes the reference: `run_coi.R` excludes `.` by default, whereas `run_dcifer.R` treats `.` as an allele.
- **Package installation issues:** Verify installed R versions and package compatibility, especially `dcifer` and `PGEcore`.

For a complete list of command-line options, use `python3 split_allele_table_by_metadata.py --help`, `Rscript run_coi.R --help`, or `Rscript run_dcifer.R --help`.

## Related software

- [Dcifer](https://github.com/EPPIcenter/dcifer) — identity-by-descent-based genetic relatedness for polyclonal infections.
- [PGEcore](https://plasmogenepi.github.io/PGEcore/) — standardised wrappers and input/output formats for malaria genomics analyses.
- MOIRE — optional source of allele-frequency
