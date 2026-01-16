#' Initialise mapping from the first timepoint
#'
#' Use the first \code{Saber} in a list as the reference timepoint and
#' define \code{true_id} from its sample names.
#'
#' @param x Named list of \code{Saber} objects, in time order.
#'
#' @return A tibble with columns \code{true_id} and one column named after
#'   the first element of \code{x}, containing the same sample IDs.
#'
#' @importFrom tibble tibble
set_first_timepoint <- function(x) {
  nm  <- names(x)[1]
  bc1 <- barcodes(x[[1]])

  tibble::tibble(
    true_id = colnames(bc1),
    !!nm    := colnames(bc1)
  )
}

#' Match current timepoint samples to the first timepoint using distances
#'
#' Compute a one-to-one mapping from samples at a later timepoint to samples
#' at the first timepoint based on Euclidean distances in scaled barcode
#' count space.
#'
#' @param x A list of length at least 2 of \code{Saber} objects, where
#'   \code{x[[1]]} is the reference (first timepoint) and \code{x[[2]]}
#'   is the timepoint to be matched.
#'
#' @return A tibble with columns \code{true_id} (sample from first timepoint)
#'   and one column named after \code{names(x)[2]} containing matched sample
#'   IDs from that timepoint. Only one-to-one matches are retained.
#'
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter group_by slice_min mutate rename select ungroup n
#' @importFrom stats dist
match_to_first_timepoint <- function(x) {
  nm2 <- names(x)[2]

  bc1 <- barcodes(x[[1]])
  bc2 <- barcodes(x[[2]])

  new_rownames <- intersect(rownames(bc1), rownames(bc2))

  if (length(new_rownames) < 2L) {
    # Not enough shared barcodes to compute meaningful distances
    return(tibble::tibble(true_id = character(), !!nm2 := character()))
  }

  distances <- stats::dist(t(scale(cbind(
    bc1[new_rownames, , drop = FALSE],
    bc2[new_rownames, , drop = FALSE]
  ))))

  distances |>
    as.matrix() |>
    tibble::as_tibble(rownames = "sample") |>
    tidyr::pivot_longer(-sample) |>
    dplyr::filter(sample %in% colnames(bc1)) |>
    dplyr::filter(name   %in% colnames(bc2)) |>
    dplyr::group_by(name) |>
    dplyr::slice_min(order_by = value, n = 1) |>
    dplyr::group_by(sample) |>
    dplyr::mutate(n = dplyr::n()) |>
    dplyr::filter(n == 1) |>
    dplyr::rename(true_id = sample, !!nm2 := name) |>
    dplyr::select(-value, -n) |>
    dplyr::ungroup()
}

#' Recursive helper to cross-link barcodes across timepoints
#'
#' @param x Named list of \code{Saber} objects, in time order.
#' @param result Current wide mapping tibble (must contain \code{true_id}).
#'
#' @return A tibble with one row per \code{true_id} and one column per
#'   timepoint name.
#'
#' @importFrom dplyr left_join
cross_sabers_rec <- function(x, result) {
  # x is a list of Saber objects ordered by timepoint
  # result is the current wide mapping tibble (must have 'true_id')

  if (length(x) == 1L) {
    # Only the first timepoint left; nothing more to match
    return(result)
  }

  # Match current first timepoint (x[[1]]) to the next one (x[[2]])
  new_result <- match_to_first_timepoint(x)

  # If there were no usable matches (e.g. too few shared barcodes), just skip
  if (nrow(new_result) > 0L) {
    result <- dplyr::left_join(result, new_result, by = "true_id")
  }

  # Recurse with the first timepoint and the remaining later timepoints
  cross_sabers_rec(x[-2], result)
}

#' Cross-link barcoded samples across multiple \code{Saber} timepoints
#'
#' Use distances in barcode count space to link samples across a series of
#' timepoints represented by \code{Saber} objects. The first element of
#' \code{x} is treated as the reference timepoint, and all later timepoints
#' are matched directly to this reference.
#'
#' @param x Named list of \code{Saber} objects, ordered by timepoint. List
#'   names are used as timepoint labels in the output.
#'
#' @return A tibble with one row per putative individual (\code{true_id}) and
#'   one column per timepoint name, containing the sample IDs from each
#'   timepoint that best match that individual (or \code{NA} if no unambiguous
#'   match was found).
#'
#' @importFrom cli cli_abort
#' @export
cross_sabers <- function(x) {
  if (is.null(names(x))) cli::cli_abort("You must supply a named list to this function.")
  base <- set_first_timepoint(x)
  cross_sabers_rec(x, base)
}

#' Visualize cross-timepoint barcode matching as a heatmap
#'
#' Given the output of \code{\link{cross_sabers}()}, produce a heatmap-style
#' plot showing which samples at each timepoint are linked to each
#' \code{true_id}.
#'
#' @param mapping A tibble as returned by \code{\link{cross_sabers}()},
#'   containing at least \code{true_id} and one or more timepoint columns.
#' @param show_labels Logical; if \code{TRUE}, sample IDs are shown as text
#'   labels inside tiles. Defaults to \code{FALSE}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_manual
#' @importFrom ggplot2 labs theme_bw theme element_text
plot_cross_sabers_heatmap <- function(mapping, show_labels = FALSE) {
  if (!"true_id" %in% colnames(mapping)) {
    stop("`mapping` must contain a 'true_id' column (output of cross_sabers()).")
  }

  tp_cols <- setdiff(colnames(mapping), "true_id")
  if (length(tp_cols) == 0L) {
    stop("`mapping` must contain at least one timepoint column besides 'true_id'.")
  }

  df <- mapping |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(tp_cols),
      names_to  = "timepoint",
      values_to = "sample_id"
    ) |>
    dplyr::mutate(
      matched = !is.na(sample_id),
      true_id = factor(true_id, levels = sort(unique(true_id)))
    )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = timepoint, y = true_id, fill = matched)
  ) +
    ggplot2::geom_tile(color = "grey90") +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "white", "TRUE" = "steelblue"),
      labels = c("FALSE" = "no match", "TRUE" = "matched"),
      name   = NULL
    ) +
    ggplot2::labs(
      x = "Timepoint",
      y = "true_id (reference sample)",
      title = "Cross-timepoint barcode linking",
      subtitle = "Tiles indicate samples matched to each reference (true_id)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if (show_labels) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = sample_id),
        size = 2
      )
  }

  p
}



#' Correlation between matched Saber samples across timepoints
#'
#' Given a list of \code{Saber} objects and a mapping (typically from
#' \code{cross_sabers()}), compute a correlation matrix between all matched
#' samples across all timepoints, based on shared barcodes.
#'
#' @param x List of \code{Saber} objects.
#' @param mapping Tibble/data frame with a \code{true_id} column and one
#'   column per timepoint containing sample IDs (e.g. output of
#'   \code{cross_sabers()}).  Optionally pre-filter the tibble (e.g. exclude unmatched samples).
#' @param method Correlation method: \code{"pearson"} or \code{"spearman"}.
#'
#' @return A numeric matrix of sample–sample correlations.
#'
#' @export
#' @importFrom purrr map reduce
#' @importFrom dplyr select filter pull
#' @importFrom tidyr pivot_longer
#' @importFrom stats cor
saber_cor <- function(x,
                      mapping,
                      method = c("pearson", "spearman")) {

  method    <- match.arg(method)

  # Common barcodes across all Saber objects ------------------------------
  rn_list <- purrr::map(x, \(s) rownames(barcodes(s)))
  common_rn <- purrr::reduce(rn_list, intersect)

  if (length(common_rn) == 0L) {
    stop("No common barcodes across all Saber objects.")
  }

  # Build combined matrix: rows = common barcodes, cols = all samples -------
  mat_list <- purrr::map(x, \(s, rn = common_rn) {
    bc <- barcodes(s)
    bc[rn, , drop = FALSE]
  })

  mat <- purrr::reduce(mat_list, cbind)

  # Get the list of mapped samples from mapping -----------------------------
  samples <- mapping |>
    dplyr::select(-true_id) |>
    tidyr::pivot_longer(dplyr::everything()) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::pull(value)

  samples <- unique(samples)

  if (length(samples) == 0L) {
    stop("No mapped samples found in `mapping` (all entries are NA?).")
  }

  # Keep only columns that actually exist in the combined matrix ------------
  samples <- intersect(samples, colnames(mat))
  if (length(samples) == 0L) {
    stop("None of the mapped sample names are present in the combined matrix.")
  }

  mat <- mat[, samples, drop = FALSE]


  # Correlation -------------------------------------------------------------
  stats::cor(mat, use = "pairwise.complete.obs", method = method)
}




