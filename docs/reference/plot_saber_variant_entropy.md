# Variant recurrence vs entropy diagnostic plot for a `Saber` object

Create a scatter plot of the number of samples in which each variant is
observed versus its normalized entropy across samples. Variants are
colored by their `is_recurrent` flag, and optional threshold lines
indicate the recurrence and entropy cutoffs used when calling
stereotyped junk.

## Usage

``` r
plot_saber_variant_entropy(x, show_thresholds = TRUE)
```

## Arguments

- x:

  A `Saber` object.

- show_thresholds:

  Logical; if `TRUE` (default), draw vertical and horizontal lines at
  `min_samples` and `entropy_thresh` taken from `x@thresholds`, when
  available.

## Value

A `ggplot` object.

## Details

The plot is based on the `recurrent_variants` slot of `x` and only
includes variants where `is_no_variant == FALSE`. The x-axis is the
number of samples in which the variant is present (after thresholding),
and the y-axis is the normalized entropy across samples. Variants
flagged as recurrent stereotyped junk (`is_recurrent == TRUE`) are
typically found in the upper-right region of the plot.
