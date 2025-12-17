# Column-wise normalized entropy (sample-level)

Compute the normalized Shannon entropy of a numeric vector, typically
representing counts of all variants within a single sample. Entropy is
computed on the non-zero entries only and normalized to \\\[0, 1\]\\.

## Usage

``` r
col_entropy_normalized(x)
```

## Arguments

- x:

  Numeric vector of counts for all variants in one sample.

## Value

A single numeric scalar in \\\[0, 1\]\\ giving the normalized entropy.

## Details

This helper is used internally to summarize per-sample diversity of
barcodes before and after filtering stereotyped recurrent variants.
