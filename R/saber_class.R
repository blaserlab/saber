#' Saber class for filtered CRISPR barcode analyses
#'
#' The \code{Saber} class is an S7 object that encapsulates the results of
#' CRISPR barcode filtering and quality control, including the filtered
#' variant count matrix, recurrent (stereotyped) variants, per-sample and
#' per-variant QC summaries, and the thresholds/parameters used to construct
#' the object.
#'
#' @details
#' Objects of class \code{Saber} are typically constructed via
#' \code{\link{make.Saber}()}, which performs empirical threshold estimation,
#' permutation-based detection of recurrent stereotyped variants, and
#' downstream QC summarization.
#'
#' The class is defined with the following properties:
#' \describe{
#'   \item{description}{Character scalar describing the dataset or analysis.}
#'   \item{filtered_counts}{Matrix of variant counts with stereotyped junk and
#'     "no variant" rows removed (variants x samples).}
#'   \item{recurrent_variants}{Tibble/data frame of per-variant recurrence and
#'     entropy metrics, including flags for recurrent stereotyped variants and
#'     "no variant" rows.}
#'   \item{motif_blacklist}{Character vector of variant identifiers treated as
#'     uninformative (either recurrent stereotyped variants or "no variant"
#'     background).}
#'   \item{per_sample_uniques}{Tibble/data frame summarizing per-sample counts
#'     of unique informative barcodes.}
#'   \item{per_sample_qc}{Tibble/data frame of per-sample QC metrics, including
#'     total reads, fractions in kept/recurrent/no-variant reads, and sample
#'     entropy before and after filtering.}
#'   \item{per_variant_qc}{Tibble/data frame of per-variant collision
#'     probabilities and related metrics.}
#'   \item{thresholds}{List of thresholds and analysis parameters used to
#'     construct the object (e.g., count/VAF thresholds, minimum sample depth,
#'     permutation settings).}
#' }
#'
#' The validator checks that \code{filtered_counts} is a matrix and that
#' \code{recurrent_variants}, \code{per_sample_uniques}, \code{per_sample_qc},
#' and \code{per_variant_qc} are tibbles (\code{"tbl_df"}).
#'
#' @seealso \code{\link{make.Saber}()}
#'
#' @keywords internal
#' @importFrom S7 new_class class_character class_any class_list
Saber <- S7::new_class(
  name = "Saber",
  properties = list(
    description         = class_character,
    filtered_counts     = class_any,
    recurrent_variants  = class_any,
    motif_blacklist     = class_character,
    per_sample_uniques  = class_any,
    per_sample_qc       = class_any,
    per_variant_qc      = class_any,
    thresholds          = class_list
  ),
  validator = function(self) {
    if (!inherits(self@filtered_counts, "matrix")) "@filtered_counts must be a matrix."
    if (!inherits(self@recurrent_variants, "tbl_df")) "@recurrent_variants must be a tibble."
    if (!inherits(self@per_sample_uniques, "tbl_df")) "@per_sample_uniques must be a tibble."
    if (!inherits(self@per_sample_qc, "tbl_df")) "@per_sample_qc must be a tibble."
    if (!inherits(self@per_variant_qc, "tbl_df")) "@per_variant_qc must be a tibble."
  }
)

#' Print method for \code{Saber} objects
#'
#' A simple print method for \code{Saber} objects that displays a short
#' summary including the description, number of samples, and the number of
#' filtered variants.
#'
#' @param x A \code{Saber} object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{Saber}}, \code{\link{make.Saber}()}
#'
#' @keywords internal
#' @name print
#' @importFrom S7 method
S7::method(print, Saber) <- function(x, ...) {
  cat("An S7 Saber object:\n")
  cat("  Description:  ", x@description, "\n")
  cat("  Samples:", ncol(x@filtered_counts), "\n")
  cat("  Filtered variants:", nrow(x@filtered_counts),  "\n")
  invisible(x)
}

#' Extract filtered barcode count matrix from a \code{Saber} object
#'
#' Generic accessor for retrieving the filtered barcode count matrix
#' (variants x samples) from a \code{Saber} object. For \code{Saber},
#' this returns the \code{@filtered_counts} slot.
#'
#' @param x An object from which barcodes can be extracted.
#'
#' @return For a \code{Saber} object, a numeric matrix of filtered barcode
#'   counts with stereotyped junk and "no variant" rows removed.
#'
#' @seealso \code{\link{Saber}}, \code{\link{make.Saber}()}
#'
#' @keywords internal
#' @importFrom S7 new_generic method
barcodes <- S7::new_generic("barcodes", "x")

#' @rdname barcodes
#'
#' @param x A \code{Saber} object.
#'
#' @keywords internal
#' @name barcodes
S7::method(barcodes, Saber) <- function(x) {
  x@filtered_counts
}

#' Extract per-sample QC data from a \code{Saber} object
#'
#' Generic accessor for retrieving sample-level quality control metrics from
#' a \code{Saber} object. For \code{Saber}, this returns the
#' \code{@per_sample_qc} slot.
#'
#' @param x An object from which sample-level QC data can be extracted.
#'
#' @return For a \code{Saber} object, a tibble/data frame with per-sample
#'   QC metrics, including total reads, fractions of reads in kept/recurrent/
#'   no-variant categories, unique barcodes, and entropy before/after
#'   filtering.
#'
#' @seealso \code{\link{Saber}}, \code{\link{make.Saber}()}
#'
#' @keywords internal
sample_data <- S7::new_generic("sample_data", "x")

#' @rdname sample_data
#'
#' @param x A \code{Saber} object.
#'
#' @keywords internal
#' @name sample_data
S7::method(sample_data, Saber) <- function(x) {
  x@per_sample_qc
}

#' Extract per-variant QC data from a \code{Saber} object
#'
#' Generic accessor for retrieving variant-level QC statistics from a
#' \code{Saber} object. For \code{Saber}, this returns the
#' \code{@per_variant_qc} slot.
#'
#' @param x An object from which variant-level QC data can be extracted.
#'
#' @return For a \code{Saber} object, a tibble/data frame with per-variant
#'   QC metrics, including collision probabilities and related statistics.
#'
#' @seealso \code{\link{Saber}}, \code{\link{make.Saber}()}
#'
#' @keywords internal
variant_data <- S7::new_generic("variant_data", "x")

#' @rdname variant_data
#'
#' @param x A \code{Saber} object.
#'
#' @keywords internal
#' @name variant_data
S7::method(variant_data, Saber) <- function(x) {
  x@per_variant_qc
}

#' Extract motif/variant blacklist from a \code{Saber} object
#'
#' Generic accessor for retrieving the blacklist of uninformative barcodes
#' (variants) from a \code{Saber} object. For \code{Saber}, this returns the
#' \code{@motif_blacklist} slot.
#'
#' @param x An object from which a blacklist of variants can be extracted.
#'
#' @return For a \code{Saber} object, a character vector of variant IDs
#'   corresponding to stereotyped recurrent and/or "no variant" motifs that
#'   should be treated as background.
#'
#' @seealso \code{\link{Saber}}, \code{\link{make.Saber}()}
#'
#' @keywords internal
blacklist <- S7::new_generic("blacklist", "x")

#' @rdname blacklist
#'
#' @param x A \code{Saber} object.
#'
#' @keywords internal
#' @name blacklist
S7::method(blacklist, Saber) <- function(x) {
  x@motif_blacklist
}
