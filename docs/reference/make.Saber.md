# Construct a `Saber` object from CRISPR barcode data

Build a `Saber` object from a `CrispRVariants` `crispr_set` by
estimating empirical noise thresholds, identifying stereotyped recurrent
barcode variants using a permutation-based test and entropy filter, and
summarizing per-sample and per-variant quality metrics.

## Usage

``` r
make.Saber(
  crispr_set,
  description = "No description provided.",
  min_sample_depth = 1000,
  background_count = 100,
  background_freq = 0.001,
  depths = NULL,
  no_variant_pattern = "^no variant$",
  thresholds = NULL,
  min_samples = 2L,
  alpha = 0.05,
  background = c("singletons", "all"),
  noise_quantile = 0.99,
  entropy_thresh = 0.8,
  B = 1000L,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- crispr_set:

  A `CrispRVariants` `crisprSet` (or similar) object from which variant
  count matrices can be extracted via
  [`CrispRVariants::variantCounts()`](https://rdrr.io/pkg/CrispRVariants/man/variantCounts.html).

- description:

  Character scalar describing the analysis or dataset.

- min_sample_depth:

  Numeric scalar; samples with total depth (column sum) below this value
  are excluded prior to analysis.

- background_count:

  Integer count threshold used when computing per-variant collision
  probabilities (background model).

- background_freq:

  Numeric VAF threshold used when computing per-variant collision
  probabilities (background model).

- depths:

  Optional numeric vector of per-sample depths. If `NULL`, computed as
  `colSums(count_mat)`.

- no_variant_pattern:

  Character scalar, regular expression used to identify rows
  corresponding to "no variant" background.

- thresholds:

  Optional list as returned by
  [`estimate_barcode_thresholds()`](estimate_barcode_thresholds.md). If
  `NULL`, thresholds are estimated from the data (recommended).

- min_samples:

  Integer; minimum number of samples in which a variant must be present
  (after thresholding) to be considered for recurrent filtering.

- alpha:

  Numeric scalar; FDR threshold applied to permutation-based empirical
  p-values when calling recurrent variants.

- background:

  Character scalar specifying how to define background when estimating
  noise thresholds; either `"singletons"` or `"all"`.

- noise_quantile:

  Numeric scalar in \\(0, 1)\\; quantile of empirical count and VAF
  distributions used as noise upper bound.

- entropy_thresh:

  Numeric scalar in \\\[0, 1\]\\; variants with normalized entropy
  greater than or equal to this value (i.e. more evenly distributed
  across samples) are more likely to be treated as stereotyped junk when
  also recurrent and significant in the permutation test.

- B:

  Integer; number of permutations used to estimate empirical p-values
  for variant recurrence.

- seed:

  Optional integer random seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) for reproducibility
  of the permutation procedure.

- progress:

  Logical; if `TRUE`, progress messages are printed during the
  permutation loop.

## Value

A `Saber` object (S7/R7 class defined elsewhere) containing:

- `description`: Analysis description.

- `filtered_counts`: Matrix of counts with stereotyped junk and
  no-variant rows removed.

- `recurrent_variants`: Data frame of per-variant recurrence and entropy
  statistics.

- `motif_blacklist`: Character vector of variant IDs treated as
  uninformative (recurrent or no-variant).

- `per_sample_uniques`: Data frame of unique barcodes per sample.

- `per_sample_qc`: Data frame of per-sample QC metrics.

- `per_variant_qc`: Tibble of per-variant collision probabilities based
  on the background model.

- `thresholds`: List of thresholds and analysis parameters.

## Details

The workflow implemented by `make.Saber()` is:

1.  Extract a variant count matrix from `crispr_set` and drop samples
    with depth below `min_sample_depth`.

2.  Estimate empirical `count_thresh` and `freq_thresh` using
    [`estimate_barcode_thresholds()`](estimate_barcode_thresholds.md) if
    `thresholds` is `NULL`.

3.  Define variant presence using these thresholds, and exclude
    `"no variant"` rows.

4.  For each variant, compute the number of samples in which it is
    present and its normalized entropy across samples.

5.  Perform a permutation-based recurrence test by shuffling the
    presence matrix within columns `B` times, deriving empirical
    p-values and FDR-adjusted p-values.

6.  Define stereotyped recurrent variants (`is_recurrent`) as those
    that (i) are not `"no variant"`, (ii) are present in at least
    `min_samples` samples, (iii) are significant at FDR `alpha`,
    and (iv) have entropy greater than or equal to `entropy_thresh`.

7.  Remove `"no variant"` rows and recurrent stereotyped variants to
    obtain `filtered_counts`.

8.  Summarize per-sample metrics (total reads, fractions in kept,
    recurrent, and no-variant reads, unique barcodes, entropy before and
    after filtering).

9.  Compute per-variant collision probabilities on the filtered matrix
    using
    [`estimate_barcode_collision_probs()`](estimate_barcode_collision_probs.md)
    with `background_count` and `background_freq`.

10. Construct and return a `Saber` object with all of these components.

## See also

[`CrispRVariants::variantCounts()`](https://rdrr.io/pkg/CrispRVariants/man/variantCounts.html),
[`estimate_barcode_thresholds()`](estimate_barcode_thresholds.md),
[`estimate_barcode_collision_probs()`](estimate_barcode_collision_probs.md)
