# Row-wise normalized entropy (variant-level)

Compute the normalized Shannon entropy of a numeric vector, typically
representing counts of a single variant across samples. Entropy is
computed on the non-zero entries only and normalized to lie in \\\[0,
1\]\\ by dividing by \\\log(m)\\, where \\m\\ is the number of non-zero
entries.

## Usage

``` r
row_entropy_normalized(x)
```

## Arguments

- x:

  Numeric vector of counts for a single variant across samples.

## Value

A single numeric scalar in \\\[0, 1\]\\ giving the normalized entropy.

## Details

This helper is intended for internal use when quantifying how evenly a
given barcode/variant is distributed across samples. Values near 0
indicate that counts are concentrated in a single sample (clone-like),
whereas values near 1 indicate counts are spread evenly across many
samples (stereotyped motif-like).
