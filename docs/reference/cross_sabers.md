# Cross-link barcoded samples across multiple `Saber` timepoints

Use distances in barcode count space to link samples across a series of
timepoints represented by `Saber` objects. The first element of `x` is
treated as the reference timepoint, and all later timepoints are matched
directly to this reference.

## Usage

``` r
cross_sabers(x)
```

## Arguments

- x:

  Named list of `Saber` objects, ordered by timepoint. List names are
  used as timepoint labels in the output.

## Value

A tibble with one row per putative individual (`true_id`) and one column
per timepoint name, containing the sample IDs from each timepoint that
best match that individual (or `NA` if no unambiguous match was found).
