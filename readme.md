# ReadTRAIL

ReadTRAIL is a two-step MLVA/VNTR analysis workflow:

1. `ReadTRAIL-Inference.py` processes FASTQ files, filters reads with BBDuk, applies strict primer validation, and writes MLVA analysis tables.
2. `ReadTRAIL-Statistics.py` builds allele/locus/sample statistics, publication tables, interpretation notes, and optional Plotly HTML plots.

## Requirements

- Conda or Mamba.

If Conda is not installed, follow the official Conda installation guide:
<https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>

## Manual Installation

Install the complete ReadTRAIL environment:

```bash
conda env create -f environment.yml
conda activate readtrail
```

This installs Python, all Python dependencies, and `bbmap`, which provides `bbduk.sh`. The user does not need to edit paths in the code.

Quick checks:

```bash
python -c "import Bio, numpy, openpyxl, pandas, plotly, psutil, tqdm"
bbduk.sh --help
```

## Primer File Format

The primer file should contain one primer pair per line, preferably as tab-separated columns:

```text
<locus_name>_<pattern_size>bp_<insert_size_in_reference_genome>bp_<corresponding_allele_coding_convention>U    <forward_primer>    <reverse_primer>
```

The locus name must follow this structure because ReadTRAIL uses it to calculate allele sizes:

- `<locus_name>`: marker or locus name.
- `<pattern_size>bp`: repeat unit size in bp.
- `<insert_size_in_reference_genome>bp`: reference amplicon size in bp.
- `<corresponding_allele_coding_convention>U`: reference allele value in repeat units.

Important: the indicated insert size must include both primer sequences.

Example:

```text
VNTR12_6bp_238bp_14U    GCTTACGACATCGTTGACAA    TCGATGGTACGCTTCTTGAT
```

## Inference Usage

```bash
python ReadTRAIL-Inference.py \
  -i /path/to/fastq_directory \
  -o /path/to/inference_output \
  -p /path/to/primers.txt \
  -m 2 \
  -t 16 \
  --genome-size 1900000
```

Main inference parameters:

| Parameter | Description |
| --- | --- |
| `-i`, `--input` | Directory with FASTQ files (`.fastq`, `.fq`, `.fastq.gz`). |
| `-o`, `--output` | Directory for inference output files. |
| `-p`, `--primer` | Primer file; see "Primer File Format" above. |
| `--genome-size` | Required genome size in bp. |
| `-m`, `--mismatch` | Max mismatches allowed; default is `2`. |
| `-t`, `--threads` | Threads per sample; default is `4`. |
| `-b`, `--binning` | Optional binning file. |
| `--flanking-seq` | Optional flanking sequence length; default is `0`. |
| `--p-min` | Detectability threshold; default is `0.05`. |
| `--obs-threshold` | Observable class threshold; default is `0.5`. |

FASTQ naming examples supported by the inference script:

- `Sample_R1_001.fastq.gz` and `Sample_R2_001.fastq.gz`
- `Sample_R1.fastq.gz` and `Sample_R2.fastq.gz`
- `Sample_1.fastq.gz` and `Sample_2.fastq.gz`

## Statistics Usage

Use the MLVA analysis workbook or CSV produced by inference:

```bash
python ReadTRAIL-Statistics.py \
  -i /path/to/inference_output/MLVA_analysis_fastq.xlsx \
  -o /path/to/statistics_output
```

Statistics output includes:

- `master_allele_to_sample_statistics.xlsx`: main Excel workbook with allele-level, locus-level, sample-level statistics, output dictionary, and column description sheets.
- `publication_tables.xlsx`: compact publication-oriented workbook with reliable allele summaries and major allele interpretation tables.
- `output_dictionary.tsv`: tab-separated dictionary describing output columns and formulas.
- `run_metadata.txt`: run summary with input path, parameters, thresholds, sample/locus counts, and output file list.
- `interpretation_notes.txt`: guide explaining how to interpret `Idealized_Share_%`, `Experimental_Share_%`, and reliability/QC columns.
- `hidden_allele_interpretation.txt`: guide and per-sample summary for hidden-allele and mixture flags.
- `plots_idealized_share.html`: per-sample interactive Plotly bar plots for `Idealized_Share_%`.
- `all_samples_idealized_alleles.html`: one interactive Plotly overview of `Idealized_Share_%` across all samples and loci.
- `plots_experimental_share.html`: per-sample interactive Plotly bar plots for `Experimental_Share_%`.
- `all_samples_experimental_alleles.html`: one interactive Plotly overview of `Experimental_Share_%` across all samples and loci.
