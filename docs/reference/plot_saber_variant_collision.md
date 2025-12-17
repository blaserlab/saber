# Variant collision risk plot for a `Saber` object

## Usage

``` r
plot_saber_variant_collision(x, p_ge2_thresh = 0.05)
```

## Arguments

- x:

  A `Saber` object.

- p_ge2_thresh:

  Numeric scalar; variants with `P_ge2_given_present >= p_ge2_thresh`
  are highlighted. Defaults to `0.05` (5\\

A `ggplot` object. Visualize per-variant collision probabilities to
identify barcode variants that are likely to arise more than once
independently in a sample. This function uses the `per_variant_qc` slot
of `x` (via [`variant_data()`](variant_data.md)). The primary plot is a
scatter plot of `n_samples_present` (x-axis) versus
`P_ge2_given_present` (y-axis), where `P_ge2_given_present` is the
model-based probability that a variant arises at least twice in a
sample, conditional on being observed at least once.A horizontal dashed
line is drawn at `p_ge2_thresh`. Axis limits are set to the observed
ranges of `n_samples_present` and `P_ge2_given_present` so that all data
points are included in the plotting region.
