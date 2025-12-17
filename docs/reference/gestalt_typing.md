# Type Gestalt Barcoding DAta

Generates typing results in the specified directory. A time stamp is
appended to the directory path provided.

## Usage

``` r
gestalt_typing(
  output_folder,
  config,
  multiplex = FALSE,
  trimmomatic_path = "/opt/Trimmomatic-0.39/trimmomatic-0.39.jar"
)
```

## Arguments

- output_folder:

  Folder for output

- config:

  Config file for setting up the analysis. See
  extdata/gestalt_config.csv

- trimmomatic_path:

  PARAM_DESCRIPTION, Default:
  '/opt/Trimmomatic-0.39/trimmomatic-0.39.jar'

## Value

Does not return anything
