# Recursive matching of samples across groups based on top-`n` barcodes

Internal helper used by [`cross_sabers`](cross_sabers.md)`()` to
recursively link samples across groups by shared top-`n` barcodes. The
first group in `input$group` is treated as the reference group.

## Usage

``` r
rec_name(input, group_order, current_idx, n_top, min_match_score)
```

## Arguments

- input:

  Tibble/data frame containing at least a `group` column, a
  `sample_name` column, and a set of `top_*` columns produced by
  [`top_n_rownames_by_column()`](top_n_rownames_by_column.md).

- n_top:

  Integer scalar; number of `top_*` columns used (passed through from
  [`cross_sabers()`](cross_sabers.md)).

- min_match_score:

  Integer scalar; minimum number of shared top barcodes required to call
  a match (passed through from [`cross_sabers()`](cross_sabers.md)).

## Value

A tibble of matched samples, with one or more rows per sample depending
on ties, including `sample_name`, `first_name` (reference ID), `top_*`
columns, and `match_score`.

## Details

The function:

1.  Identifies the first group in the data as the reference group.

2.  Builds a table of reference samples and their `top_*` barcodes.

3.  Performs a cross join between all samples and reference samples,
    computing, for each pair, the number of shared non-`NA` `top_*`
    barcodes (`match_score`).

4.  For each sample, retains the reference with the maximum
    `match_score`, and filters matches below `min_match_score`.

5.  Removes matched samples and recurses on the remaining input until no
    unmatched samples remain.
