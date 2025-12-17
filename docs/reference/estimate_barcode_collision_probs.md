# Estimate per-variant barcode collision probabilities

Estimate, for each variant, the probability that it arises multiple
times independently in a single sample under a simple Poisson model,
using the observed sharing of that variant across samples.

## Usage

``` r
estimate_barcode_collision_probs(count_mat, count_thresh, freq_thresh)
```

## Arguments

- count_mat:

  Numeric matrix of variant counts (rows = variants, columns = samples).

- count_thresh:

  Integer minimum count threshold for presence, used to construct the
  presence matrix.

- freq_thresh:

  Numeric minimum VAF threshold for presence, used to construct the
  presence matrix.

## Value

A tibble with one row per variant containing:

- variant:

  Variant identifier.

- n_samples_present:

  Number of samples with the variant present.

- p_sample_present:

  Estimated probability a random sample contains the variant at least
  once.

- lambda_per_sample:

  Estimated Poisson mean per sample.

- P_ge2_per_sample:

  Unconditional probability \\P(K \ge 2)\\.

- P_ge2_given_present:

  Conditional probability \\P(K \ge 2 \mid K \ge 1)\\.

## Details

A logical presence matrix is first constructed using
[`generate_present_mat()`](generate_present_mat.md). For each variant
\\i\\, the number of samples in which it is present \\m_i\\ is used to
estimate \\\hat{p}\_i = m_i / S\\, where \\S\\ is the number of samples.
Under a Poisson model for the number of independent origins per sample,
\\K_i \sim \mathrm{Poisson}(\lambda_i)\\, we have \\\hat{p}\_i \approx
P(K_i \ge 1) = 1 - e^{-\lambda_i}\\, so \\\hat{\lambda}\_i = -\log(1 -
\hat{p}\_i)\\. From this, `P_ge2_per_sample` and `P_ge2_given_present`
are computed.
