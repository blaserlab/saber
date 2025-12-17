# Identify top-`n` variants per sample

For each column (sample) of a numeric matrix, return the row names
corresponding to the top `n` values. The result is a tibble with one row
per sample and columns `top_1`, `top_2`, ..., `top_n`.

## Usage

``` r
top_n_rownames_by_column(mat, n)
```

## Arguments

- mat:

  Numeric matrix of values (rows = variants, columns = samples).

- n:

  Integer scalar; number of top variants to select per sample.

## Value

A tibble with one row per sample, a `sample_name` column, and `top_1` to
`top_n` columns containing variant IDs.

## Details

If row names are missing, they are created as character indices (`"1"`,
`"2"`, ...). If column names are missing, they are created as `"col1"`,
`"col2"`, etc. For each sample, the row names of the top `n` entries are
returned in descending order of the column values. If a sample has fewer
than `n` non-`NA` entries, the remaining positions are padded with `NA`.
