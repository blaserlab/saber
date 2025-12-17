# Extract per-sample QC data from a `Saber` object

Generic accessor for retrieving sample-level quality control metrics
from a `Saber` object. For `Saber`, this returns the `@per_sample_qc`
slot.

## Usage

``` r
sample_data(x, ...)
```

## Arguments

- x:

  A `Saber` object.

## Value

For a `Saber` object, a tibble/data frame with per-sample QC metrics,
including total reads, fractions of reads in kept/recurrent/ no-variant
categories, unique barcodes, and entropy before/after filtering.

## See also

[`Saber`](Saber.md), [`make.Saber()`](make.Saber.md)
