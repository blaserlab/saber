# Sample-level QC plot for a `Saber` object

Visualize per-sample quality control metrics, highlighting samples with
a low fraction of reads retained after filtering recurrent and
no-variant barcodes. Outliers are defined as samples whose
`frac_reads_kept` is at least 3 median absolute deviations (MADs) below
the median.

## Usage

``` r
plot_saber_sample_qc(x)
```

## Arguments

- x:

  A `Saber` object.

## Value

A `ggplot` object.

## Details

This function uses the `per_sample_qc` slot of `x` (via
[`sample_data()`](sample_data.md)) and produces a bar plot of
`frac_reads_kept` per sample, ordered by this fraction.

The outlier threshold is computed as: \$\$ T =
\mathrm{median}(\mathrm{frac\\reads\\kept}) - 3 \times
\mathrm{MAD}(\mathrm{frac\\reads\\kept}) \$\$ using the usual MAD
definition. Samples with `frac_reads_kept < T` are highlighted in a
different color, and a horizontal dashed line is drawn at `T`.
