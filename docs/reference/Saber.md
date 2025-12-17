# Saber class for filtered CRISPR barcode analyses

The `Saber` class is an S7 object that encapsulates the results of
CRISPR barcode filtering and quality control, including the filtered
variant count matrix, recurrent (stereotyped) variants, per-sample and
per-variant QC summaries, and the thresholds/parameters used to
construct the object.

## Usage

``` r
Saber(
  description = character(0),
  filtered_counts = NULL,
  recurrent_variants = NULL,
  motif_blacklist = character(0),
  per_sample_uniques = NULL,
  per_sample_qc = NULL,
  per_variant_qc = NULL,
  thresholds = list()
)
```

## Details

Objects of class `Saber` are typically constructed via
[`make.Saber()`](make.Saber.md), which performs empirical threshold
estimation, permutation-based detection of recurrent stereotyped
variants, and downstream QC summarization.

The class is defined with the following properties:

- description:

  Character scalar describing the dataset or analysis.

- filtered_counts:

  Matrix of variant counts with stereotyped junk and "no variant" rows
  removed (variants x samples).

- recurrent_variants:

  Tibble/data frame of per-variant recurrence and entropy metrics,
  including flags for recurrent stereotyped variants and "no variant"
  rows.

- motif_blacklist:

  Character vector of variant identifiers treated as uninformative
  (either recurrent stereotyped variants or "no variant" background).

- per_sample_uniques:

  Tibble/data frame summarizing per-sample counts of unique informative
  barcodes.

- per_sample_qc:

  Tibble/data frame of per-sample QC metrics, including total reads,
  fractions in kept/recurrent/no-variant reads, and sample entropy
  before and after filtering.

- per_variant_qc:

  Tibble/data frame of per-variant collision probabilities and related
  metrics.

- thresholds:

  List of thresholds and analysis parameters used to construct the
  object (e.g., count/VAF thresholds, minimum sample depth, permutation
  settings).

The validator checks that `filtered_counts` is a matrix and that
`recurrent_variants`, `per_sample_uniques`, `per_sample_qc`, and
`per_variant_qc` are tibbles (`"tbl_df"`).

## See also

[`make.Saber()`](make.Saber.md)
