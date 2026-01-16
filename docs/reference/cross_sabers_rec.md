# Recursive helper to cross-link barcodes across timepoints

Recursive helper to cross-link barcodes across timepoints

## Usage

``` r
cross_sabers_rec(x, result)
```

## Arguments

- x:

  Named list of `Saber` objects, in time order.

- result:

  Current wide mapping tibble (must contain `true_id`).

## Value

A tibble with one row per `true_id` and one column per timepoint name.
