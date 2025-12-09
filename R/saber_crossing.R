#' Identify top-\code{n} variants per sample
#'
#' For each column (sample) of a numeric matrix, return the row names
#' corresponding to the top \code{n} values. The result is a tibble with one
#' row per sample and columns \code{top_1}, \code{top_2}, ..., \code{top_n}.
#'
#' @param mat Numeric matrix of values (rows = variants, columns = samples).
#' @param n Integer scalar; number of top variants to select per sample.
#'
#' @details
#' If row names are missing, they are created as character indices
#' (\code{"1"}, \code{"2"}, ...). If column names are missing, they are created
#' as \code{"col1"}, \code{"col2"}, etc. For each sample, the row names of the
#' top \code{n} entries are returned in descending order of the column values.
#' If a sample has fewer than \code{n} non-\code{NA} entries, the remaining
#' positions are padded with \code{NA}.
#'
#' @return A tibble with one row per sample, a \code{sample_name} column, and
#'   \code{top_1} to \code{top_n} columns containing variant IDs.
#'
#' @keywords internal
#' @importFrom tibble tibble
top_n_rownames_by_column <- function(mat, n) {
  # Basic checks
  stopifnot(is.matrix(mat), is.numeric(mat))

  # Ensure rownames exist
  if (is.null(rownames(mat))) {
    rownames(mat) <- as.character(seq_len(nrow(mat)))
  }

  # Ensure colnames exist
  col_names <- colnames(mat)
  if (is.null(col_names)) {
    col_names <- paste0("col", seq_len(ncol(mat)))
  }

  # For each column, get rownames of the top n values
  top_list <- lapply(seq_len(ncol(mat)), function(j) {
    col_vals <- mat[, j]
    ord <- order(col_vals, decreasing = TRUE, na.last = NA)
    top_idx <- head(ord, n)
    top_rn <- rownames(mat)[top_idx]

    # Pad with NA if there are fewer than n non-NA values
    if (length(top_rn) < n) {
      top_rn <- c(top_rn, rep(NA_character_, n - length(top_rn)))
    }

    top_rn
  })

  # Bind to a matrix and convert to tibble
  res_mat <- do.call(rbind, top_list)
  colnames(res_mat) <- paste0("top_", seq_len(n))

  tibble::tibble(
    sample_name = col_names,
    as.data.frame(res_mat, stringsAsFactors = FALSE)
  )
}

#' Recursive matching of samples across groups based on top-\code{n} barcodes
#'
#' Internal helper used by \code{\link{cross_sabers}()} to recursively link
#' samples across groups by shared top-\code{n} barcodes. The first group in
#' \code{input$group} is treated as the reference group.
#'
#' @param input Tibble/data frame containing at least a \code{group} column,
#'   a \code{sample_name} column, and a set of \code{top_*} columns produced by
#'   \code{top_n_rownames_by_column()}.
#' @param n_top Integer scalar; number of \code{top_*} columns used (passed
#'   through from \code{cross_sabers()}).
#' @param min_match_score Integer scalar; minimum number of shared top barcodes
#'   required to call a match (passed through from \code{cross_sabers()}).
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Identifies the first group in the data as the reference group.
#'   \item Builds a table of reference samples and their \code{top_*} barcodes.
#'   \item Performs a cross join between all samples and reference samples,
#'         computing, for each pair, the number of shared non-\code{NA}
#'         \code{top_*} barcodes (\code{match_score}).
#'   \item For each sample, retains the reference with the maximum
#'         \code{match_score}, and filters matches below \code{min_match_score}.
#'   \item Removes matched samples and recurses on the remaining input until
#'         no unmatched samples remain.
#' }
#'
#' @return A tibble of matched samples, with one or more rows per sample
#'   depending on ties, including \code{sample_name}, \code{first_name}
#'   (reference ID), \code{top_*} columns, and \code{match_score}.
#'
#' @keywords internal
#' @importFrom dplyr pull filter select rename_with cross_join rowwise mutate
#' @importFrom dplyr ungroup group_by anti_join bind_rows all_of
rec_name <- function(input,
                     group_order,
                     current_idx,
                     n_top,
                     min_match_score) {
  # Base cases -------------------------------------------------------------
  if (nrow(input) == 0L) {
    # Nothing left to match
    return(input[0, , drop = FALSE])
  }
  if (current_idx > length(group_order)) {
    # No more reference groups to use
    return(input[0, , drop = FALSE])
  }

  # Current reference group (fixed order from group_order)
  ref_group <- group_order[current_idx]

  top_cols <- grep("^top_", names(input), value = TRUE)

  # Reference ID table from the current reference group
  id_table <- input |>
    dplyr::filter(group == ref_group) |>
    dplyr::select(first_name = sample_name, dplyr::all_of(top_cols))

  # If this reference group has no remaining samples, skip to next group
  if (nrow(id_table) == 0L) {
    return(rec_name(
      input       = input,
      group_order = group_order,
      current_idx = current_idx + 1L,
      n_top       = n_top,
      min_match_score = min_match_score
    ))
  }

  # Rename reference top_* columns so we can keep both sets in the join -----
  id_table_ref <- id_table |>
    dplyr::rename_with(\(nm) paste0(nm, "_ref"), dplyr::all_of(top_cols))

  # Cross join: compare every sample to every reference sample --------------
  result <- dplyr::cross_join(input, id_table_ref) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      match_score = {
        # current sample's top clones
        a <- c_across(dplyr::all_of(top_cols))
        # reference sample's top clones
        b <- c_across(dplyr::all_of(paste0(top_cols, "_ref")))
        a <- a[!is.na(a)]
        b <- b[!is.na(b)]
        length(intersect(a, b))
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(sample_name) |>
    dplyr::filter(match_score == max(match_score, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::filter(match_score >= min_match_score)

  # Remove matched samples and recurse on the remainder ---------------------
  newinput <- dplyr::anti_join(input, result, by = "sample_name")

  # Collect matches from this reference group and then move to the next one
  dplyr::bind_rows(
    result,
    rec_name(
      input       = newinput,
      group_order = group_order,
      current_idx = current_idx + 1L,
      n_top       = n_top,
      min_match_score = min_match_score
    )
  )
}


#' Cross-link barcoded samples across multiple \code{Saber} objects
#'
#' Use shared top-\code{n} barcodes to link samples across groups of
#' \code{Saber} objects, typically representing different experimental runs or
#' time points from the same set of animals (e.g., zebrafish). The first
#' element in \code{saber_list} is treated as the reference group, and samples
#' in subsequent groups are matched to this reference by overlapping top
#' barcodes.
#'
#' @param saber_list A named list of \code{Saber} objects. The names of the
#'   list are used as group identifiers.
#' @param n_top Integer scalar; number of top barcodes per sample (per group)
#'   to consider when computing matches. Defaults to \code{10}.
#' @param min_match_score Integer scalar; minimum number of shared top
#'   barcodes required to call a match between a sample and a reference
#'   sample. Defaults to \code{5}.
#'
#' @details
#' For each \code{Saber} object, the function:
#' \enumerate{
#'   \item Extracts the filtered count matrix via \code{\link{barcodes}()}.
#'   \item Uses \code{top_n_rownames_by_column()} to identify the top
#'         \code{n_top} barcodes per sample.
#'   \item Adds a \code{group} column indicating the originating list element.
#' }
#' All groups are then combined and passed to \code{rec_name()}, which
#' recursively matches samples to the first group's samples (reference) using
#' the number of shared top barcodes as a score. After matching, a check is
#' performed to ensure that each sample is matched to at most one reference
#' sample; if this condition is violated, an error is raised suggesting
#' increasing stringency or improving barcode filtering.
#'
#' @return Invisibly, the function either:
#' \itemize{
#'   \item raises an error if some samples match more than one reference, or
#'   \item (if integrated into your workflow) yields a tibble of matches as the
#'         last evaluated expression before the consistency check. The tibble
#'         contains \code{group}, \code{sample_name}, \code{true_id} (reference
#'         sample), \code{top_*} columns, and \code{match_score}.
#' }
#'
#' @seealso \code{\link{Saber}}, \code{\link{barcodes}()}
#'
#' @export
#' @importFrom purrr map2
#' @importFrom dplyr mutate relocate bind_rows select starts_with count
#' @importFrom cli cli_abort
cross_sabers <- function(saber_list,
                         n_top = 10L,
                         min_match_score = 5L) {
  # Build input table: top n barcodes per sample, per group
  input <- purrr::map2(
    .x = saber_list,
    .y = names(saber_list),
    .f = \(x, y) {
      top_n_rownames_by_column(barcodes(x), n = n_top) |>
        dplyr::mutate(group = y) |>
        dplyr::relocate(group)
    }
  ) |>
    dplyr::bind_rows()

  # Freeze group order once, based on the order they appear in input
  # (which reflects the order in saber_list)
  group_order <- input |>
    dplyr::pull(group) |>
    unique()

  # Run recursive matcher using explicit group_order and starting at index 1
  res <- rec_name(
    input       = input,
    group_order = group_order,
    current_idx = 1L,
    n_top       = n_top,
    min_match_score = min_match_score
  ) |>
    dplyr::select(
      group,
      sample_name,
      true_id = first_name,
      dplyr::starts_with("top_"),
      match_score
    )

  # Check that no sample matches more than one reference
  check <- res |>
    dplyr::count(sample_name) |>
    dplyr::pull(n)

  if (max(check) > 1) {
    cli::cli_abort(
      "Some samples matched to more than one reference. Try increasing {.arg n_top} and {.arg min_match_score}, increasing the stringency of barcode calling in the Saber objects, or filtering out samples with poor barcoding."
    )
  }

  res
}
