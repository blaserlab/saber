# Extract per-variant QC data from a `Saber` object

Generic accessor for retrieving variant-level QC statistics from a
`Saber` object. For `Saber`, this returns the `@per_variant_qc` slot.

## Usage

``` r
variant_data(x, ...)
```

## Arguments

- x:

  A `Saber` object.

## Value

For a `Saber` object, a tibble/data frame with per-variant QC metrics,
including collision probabilities and related statistics.

## See also

[`Saber`](Saber.md), [`make.Saber()`](make.Saber.md)
