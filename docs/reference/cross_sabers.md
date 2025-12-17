# Cross-link barcoded samples across multiple `Saber` objects

Use shared top-`n` barcodes to link samples across groups of `Saber`
objects, typically representing different experimental runs or time
points from the same set of animals (e.g., zebrafish). The first element
in `saber_list` is treated as the reference group, and samples in
subsequent groups are matched to this reference by overlapping top
barcodes.

## Usage

``` r
cross_sabers(saber_list, n_top = 10L, min_match_score = 5L)
```

## Arguments

- saber_list:

  A named list of `Saber` objects. The names of the list are used as
  group identifiers.

- n_top:

  Integer scalar; number of top barcodes per sample (per group) to
  consider when computing matches. Defaults to `10`.

- min_match_score:

  Integer scalar; minimum number of shared top barcodes required to call
  a match between a sample and a reference sample. Defaults to `5`.

## Value

Invisibly, the function either:

- raises an error if some samples match more than one reference, or

- (if integrated into your workflow) yields a tibble of matches as the
  last evaluated expression before the consistency check. The tibble
  contains `group`, `sample_name`, `true_id` (reference sample), `top_*`
  columns, and `match_score`.

## Details

For each `Saber` object, the function:

1.  Extracts the filtered count matrix via
    [`barcodes`](barcodes.md)`()`.

2.  Uses [`top_n_rownames_by_column()`](top_n_rownames_by_column.md) to
    identify the top `n_top` barcodes per sample.

3.  Adds a `group` column indicating the originating list element.

All groups are then combined and passed to [`rec_name()`](rec_name.md),
which recursively matches samples to the first group's samples
(reference) using the number of shared top barcodes as a score. After
matching, a check is performed to ensure that each sample is matched to
at most one reference sample; if this condition is violated, an error is
raised suggesting increasing stringency or improving barcode filtering.

## See also

[`Saber`](Saber.md), [`barcodes`](barcodes.md)`()`
