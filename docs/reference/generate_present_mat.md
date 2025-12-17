# Generate a logical presence matrix from counts and thresholds

Create a logical matrix indicating presence/absence of variants in
samples based on count and VAF thresholds.

## Usage

``` r
generate_present_mat(count_mat, count_thresh, freq_thresh)
```

## Arguments

- count_mat:

  Numeric matrix of variant counts (rows = variants, columns = samples).

- count_thresh:

  Integer minimum count threshold for presence.

- freq_thresh:

  Numeric minimum VAF threshold for presence.

## Value

A logical matrix of the same dimension as `count_mat` with `TRUE`
indicating presence.

## Details

Variant allele frequencies are computed as `count / depth`, where
`depth` is the column sum of `count_mat`. A cell is considered present
if its count is greater than or equal to `count_thresh` and its VAF is
greater than or equal to `freq_thresh`.
