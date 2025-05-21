md <- readr::read_csv("inst/extdata/references/blank_metadata.csv")
gd <- rtracklayer::import("inst/extdata/references/gestalt2.bed")
reference0 <- readr::read_file("inst/extdata/references/gestalt2_ref.txt")
usethis::use_data(md, gd, reference0, internal = TRUE, overwrite = TRUE)
