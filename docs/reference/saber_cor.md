# Correlation between matched Saber samples across timepoints

Given a list of `Saber` objects and a mapping (typically from
[`cross_sabers()`](cross_sabers.md)), compute a correlation matrix
between all matched samples across all timepoints, based on shared
barcodes.

## Usage

``` r
saber_cor(x, mapping, method = c("pearson", "spearman"))
```

## Arguments

- x:

  List of `Saber` objects.

- mapping:

  Tibble/data frame with a `true_id` column and one column per timepoint
  containing sample IDs (e.g. output of
  [`cross_sabers()`](cross_sabers.md)). Optionally pre-filter the tibble
  (e.g. exclude unmatched samples).

- method:

  Correlation method: `"pearson"` or `"spearman"`.

## Value

A numeric matrix of sample–sample correlations.
