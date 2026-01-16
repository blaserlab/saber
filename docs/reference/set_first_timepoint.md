# Initialise mapping from the first timepoint

Use the first `Saber` in a list as the reference timepoint and define
`true_id` from its sample names.

## Usage

``` r
set_first_timepoint(x)
```

## Arguments

- x:

  Named list of `Saber` objects, in time order.

## Value

A tibble with columns `true_id` and one column named after the first
element of `x`, containing the same sample IDs.
