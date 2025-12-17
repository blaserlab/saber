# Extract filtered barcode count matrix from a `Saber` object

Generic accessor for retrieving the filtered barcode count matrix
(variants x samples) from a `Saber` object. For `Saber`, this returns
the `@filtered_counts` slot.

## Usage

``` r
barcodes(x, ...)
```

## Arguments

- x:

  A `Saber` object.

## Value

For a `Saber` object, a numeric matrix of filtered barcode counts with
stereotyped junk and "no variant" rows removed.

## See also

[`Saber`](Saber.md), [`make.Saber()`](make.Saber.md)
