# Visualize cross-timepoint barcode matching as a heatmap

Given the output of [`cross_sabers()`](cross_sabers.md), produce a
heatmap-style plot showing which samples at each timepoint are linked to
each `true_id`.

## Usage

``` r
plot_cross_sabers_heatmap(mapping, show_labels = FALSE)
```

## Arguments

- mapping:

  A tibble as returned by [`cross_sabers()`](cross_sabers.md),
  containing at least `true_id` and one or more timepoint columns.

- show_labels:

  Logical; if `TRUE`, sample IDs are shown as text labels inside tiles.
  Defaults to `FALSE`.

## Value

A `ggplot` object.
