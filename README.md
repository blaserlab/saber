
# Selecting Amplicon Barcodes using Experimental Replicates

The *saber* package provides functions for genotyping GESTALT barcodes, barcode quality control, and quantifying barcodes from an experiment with dozens of samples or more.

This is an updated/modernized version of the code published in 

Teets, E.M., et al., Quantifying Hematopoietic Stem Cell Clonal Diversity by Selecting Informative Amplicon Barcodes. Sci Rep, 2020. 10(1): p. 2153.

Currently only the GESTALT genotyping function is implemented in this package.

## Installation

You can install saber like so:

``` r
pak::pkg_install("blaserlab/saber")
```

## Genotyping Bulk Gestalt Amplicon Data

You need to generate a configuration file in csv format.  Columns should be:

```
sample,group,sample_fp_1,sample_fp_2
```
Where group can identify your experimental group or be "not_applicable", but must be a character string.  "samples_fp_1" and "samples_fp_2" are the file paths to the R1 and R2 fastq files for that sample.

A blank config file can be found by running:

``` r
fs::path_package(package = "saber", "extdata", "gestalt_config.csv")
```

Then run:

``` r
library(saber)

gestalt_typing(config = "<path/to/config/file>", output_folder = "<path/to/output>")

```

