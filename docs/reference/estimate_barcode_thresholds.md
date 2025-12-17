# Estimate empirical count and VAF thresholds from noise

Estimate count and variant allele frequency (VAF) thresholds for
defining barcode presence based on the empirical distribution of
low-prevalence variants (noise) across a count matrix.

## Usage

``` r
estimate_barcode_thresholds(
  count_mat,
  depths = NULL,
  no_variant_pattern = "^no variant$",
  background = c("singletons", "all"),
  noise_quantile = 0.99
)
```

## Arguments

- count_mat:

  Numeric matrix of variant counts (rows = variants, columns = samples).

- depths:

  Optional numeric vector of per-sample depths (column sums). If `NULL`,
  computed as `colSums(count_mat)` after excluding rows matching
  `no_variant_pattern`.

- no_variant_pattern:

  Character scalar, regular expression matching rows that represent "no
  variant" background to be excluded from noise estimation.

- background:

  Character scalar specifying which cells to treat as background when
  estimating thresholds. Either `"singletons"` (only variants present in
  exactly one sample) or `"all"` (all non-zero entries).

- noise_quantile:

  Numeric scalar in \\(0, 1)\\ giving the quantile of the empirical
  count and VAF distributions to treat as the upper bound of noise.
  Values just above this quantile are used as thresholds.

## Value

A list with elements:

- count_thresh:

  Integer count threshold for presence.

- freq_thresh:

  Numeric VAF threshold for presence.

- noise_quantile:

  The noise quantile used.

- background:

  The background mode used.

- n_nonzero_used:

  Number of non-zero observations used to estimate thresholds.

## Details

The function:

1.  Excludes rows matching `no_variant_pattern`.

2.  Optionally restricts to "singleton" variants (present in exactly one
    sample) when `background = "singletons"`.

3.  Computes the specified `noise_quantile` of the empirical non-zero
    counts and VAFs and sets `count_thresh` to one plus the floored
    count quantile and `freq_thresh` to the VAF quantile.
