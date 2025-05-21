#' @title Type Gestalt Barcoding DAta
#' @description Generates typing results in the specified directory.  A time stamp is appended to the directory path provided.
#' @param output_folder Folder for output
#' @param config Config file for setting up the analysis.  See extdata/gestalt_config.csv
#' @param trimmomatic_path PARAM_DESCRIPTION, Default: '/opt/Trimmomatic-0.39/trimmomatic-0.39.jar'
#' @return Does not return anything
#' @rdname gestalt_typing
#' @export
#' @importFrom readr read_csv write_csv
#' @importFrom stringr str_replace_all
#' @importFrom fs dir_create path
gestalt_typing <- function(output_folder,
                           config,
                           trimmomatic_path = "/opt/Trimmomatic-0.39/trimmomatic-0.39.jar") {
  ts <-
    stringr::str_replace_all(Sys.time(), "[:punct:]|[:alpha:]|[:space:]", "")
  of <- fs::path(paste0(output_folder, "_", ts))
  fs::dir_create(of)
  cf <- readr::read_csv(config)
  tp <- trimmomatic_path
  make_tree(output_folder = of)
  load_gestalt_fastqs(config = cf)
  merge_read_pairs(output_folder = of, conf = cf)
  trim_reads(trimmomatic_path = tp, conf = cf)
  run_cutadapt(conf = cf)
  run_fastqc(output_folder = of)
  align_reads(output_folder = of, conf = cf)
  fix_sam_header()
  sam_to_bam(conf = cf)
  sort_and_index(output_folder = of)
  dat <- make_crispr_set(output_folder = of, conf = cf)
  plot_variants(output_folder = of,
                conf = cf,
                data = dat)
  readr::write_csv(cf, file = fs::path(of, fs::path_file(config)))

}
