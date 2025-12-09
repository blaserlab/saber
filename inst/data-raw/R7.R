Saber <- new_class(
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

method(print, Saber) <- function(x, ...) {
  cat("An R7 Saber object:\n")
  cat("  Description:  ", x@description, "\n")
  cat("  Samples:", ncol(x@filtered_counts), "\n")
  cat("  Filtered variants:", nrow(x@filtered_counts),  "\n")
  invisible(x)
}
