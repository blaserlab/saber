# Match current timepoint samples to the first timepoint using distances

Compute a one-to-one mapping from samples at a later timepoint to
samples at the first timepoint based on Euclidean distances in scaled
barcode count space.

## Usage

``` r
match_to_first_timepoint(x)
```

## Arguments

- x:

  A list of length at least 2 of `Saber` objects, where `x[[1]]` is the
  reference (first timepoint) and `x[[2]]` is the timepoint to be
  matched.

## Value

A tibble with columns `true_id` (sample from first timepoint) and one
column named after `names(x)[2]` containing matched sample IDs from that
timepoint. Only one-to-one matches are retained.
