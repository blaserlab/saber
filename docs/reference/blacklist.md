# Extract motif/variant blacklist from a `Saber` object

Generic accessor for retrieving the blacklist of uninformative barcodes
(variants) from a `Saber` object. For `Saber`, this returns the
`@motif_blacklist` slot.

## Usage

``` r
blacklist(x, ...)
```

## Arguments

- x:

  A `Saber` object.

## Value

For a `Saber` object, a character vector of variant IDs corresponding to
stereotyped recurrent and/or "no variant" motifs that should be treated
as background.

## See also

[`Saber`](Saber.md), [`make.Saber()`](make.Saber.md)
