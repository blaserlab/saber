#' Variant recurrence vs entropy diagnostic plot for a \code{Saber} object
#'
#' Create a scatter plot of the number of samples in which each variant is
#' observed versus its normalized entropy across samples. Variants are colored
#' by their \code{is_recurrent} flag, and optional threshold lines indicate
#' the recurrence and entropy cutoffs used when calling stereotyped junk.
#'
#' @param x A \code{Saber} object.
#' @param show_thresholds Logical; if \code{TRUE} (default), draw vertical and
#'   horizontal lines at \code{min_samples} and \code{entropy_thresh} taken
#'   from \code{x@thresholds}, when available.
#'
#' @details
#' The plot is based on the \code{recurrent_variants} slot of \code{x} and
#' only includes variants where \code{is_no_variant == FALSE}. The x-axis
#' is the number of samples in which the variant is present (after thresholding),
#' and the y-axis is the normalized entropy across samples. Variants flagged
#' as recurrent stereotyped junk (\code{is_recurrent == TRUE}) are typically
#' found in the upper-right region of the plot.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline scale_color_manual
#' @importFrom ggplot2 labs theme_bw
plot_saber_variant_entropy <- function(x, show_thresholds = TRUE) {
  rv <- x@recurrent_variants

  required_cols <- c("n_samples", "entropy", "is_recurrent", "is_no_variant")
  missing_cols <- setdiff(required_cols, colnames(rv))
  if (length(missing_cols) > 0L) {
    stop("recurrent_variants is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Thresholds, if available
  min_samples <- NA_real_
  entropy_thresh <- NA_real_
  if (!is.null(x@thresholds)) {
    if (!is.null(x@thresholds$min_samples)) {
      min_samples <- x@thresholds$min_samples
    }
    if (!is.null(x@thresholds$entropy_thresh)) {
      entropy_thresh <- x@thresholds$entropy_thresh
    }
  }

  df <- rv[!rv$is_no_variant, , drop = FALSE]

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = n_samples, y = entropy, color = is_recurrent)
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "grey50", "TRUE" = "red"),
      labels = c("FALSE" = "kept", "TRUE" = "recurrent junk"),
      name   = "Filter status"
    ) +
    ggplot2::labs(
      x = "Number of samples with variant (n_samples)",
      y = "Normalized entropy across samples",
      title = "Variant recurrence vs entropy",
      subtitle = "Points in red are flagged as recurrent stereotyped variants"
    ) +
    ggplot2::theme_bw()

  if (show_thresholds && !is.na(min_samples)) {
    p <- p + ggplot2::geom_vline(xintercept = min_samples, linetype = "dashed")
  }
  if (show_thresholds && !is.na(entropy_thresh)) {
    p <- p + ggplot2::geom_hline(yintercept = entropy_thresh, linetype = "dashed")
  }

  p
}

#' Sample-level QC plot for a \code{Saber} object
#'
#' Visualize per-sample quality control metrics, highlighting samples with a
#' low fraction of reads retained after filtering recurrent and no-variant
#' barcodes. Outliers are defined as samples whose \code{frac_reads_kept} is
#' at least 3 median absolute deviations (MADs) below the median.
#'
#' @param x A \code{Saber} object.
#'
#' @details
#' This function uses the \code{per_sample_qc} slot of \code{x} (via
#' \code{\link{sample_data}()}) and produces a bar plot of
#' \code{frac_reads_kept} per sample, ordered by this fraction.
#'
#' The outlier threshold is computed as:
#' \deqn{
#'   T = \mathrm{median}(\mathrm{frac\_reads\_kept}) - 3 \times \mathrm{MAD}(\mathrm{frac\_reads\_kept})
#' }
#' using the usual MAD definition. Samples with
#' \code{frac_reads_kept < T} are highlighted in a different color, and a
#' horizontal dashed line is drawn at \code{T}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_col coord_flip geom_hline scale_fill_manual
#' @importFrom ggplot2 labs theme_bw
#' @importFrom stats median mad
plot_saber_sample_qc <- function(x) {
  qc <- sample_data(x)

  required_cols <- c("sample", "frac_reads_kept")
  missing_cols <- setdiff(required_cols, colnames(qc))
  if (length(missing_cols) > 0L) {
    stop("per_sample_qc is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  df <- qc

  med <- stats::median(df$frac_reads_kept, na.rm = TRUE)
  m_ad <- stats::mad(df$frac_reads_kept, na.rm = TRUE)
  thresh <- med - 3 * m_ad

  df$low_kept <- df$frac_reads_kept < thresh
  df$sample <- factor(df$sample, levels = df$sample[order(df$frac_reads_kept)])

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = sample, y = frac_reads_kept, fill = low_kept)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = thresh, linetype = "dashed") +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "grey70", "TRUE" = "red"),
      labels = c("FALSE" = "within 3 MAD of median", "TRUE" = "< median - 3 MAD"),
      name   = NULL
    ) +
    ggplot2::labs(
      x = "Sample",
      y = "Fraction of reads kept",
      title = "Sample-level QC: fraction of reads retained",
      subtitle = paste0(
        "Red samples have frac_reads_kept < median - 3 MAD (threshold = ",
        signif(thresh, 3), ")"
      )
    ) +
    ggplot2::theme_bw()
}

#' Variant collision risk plot for a \code{Saber} object
#'
#' Visualize per-variant collision probabilities to identify barcode variants
#' that are likely to arise more than once independently in a sample.
#'
#' @param x A \code{Saber} object.
#' @param p_ge2_thresh Numeric scalar; variants with
#'   \code{P_ge2_given_present >= p_ge2_thresh} are highlighted. Defaults to
#'   \code{0.05} (5\%).
#'
#' @details
#' This function uses the \code{per_variant_qc} slot of \code{x} (via
#' \code{\link{variant_data}()}). The primary plot is a scatter plot of
#' \code{n_samples_present} (x-axis) versus \code{P_ge2_given_present}
#' (y-axis), where \code{P_ge2_given_present} is the model-based probability
#' that a variant arises at least twice in a sample, conditional on being
#' observed at least once.
#'
#' A horizontal dashed line is drawn at \code{p_ge2_thresh}. Axis limits are
#' set to the observed ranges of \code{n_samples_present} and
#' \code{P_ge2_given_present} so that all data points are included in the
#' plotting region.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_hline labs theme_bw
#' @importFrom ggplot2 scale_y_continuous scale_x_continuous
plot_saber_variant_collision <- function(x, p_ge2_thresh = 0.05) {
  vc <- variant_data(x)

  required_cols <- c("n_samples_present", "P_ge2_given_present")
  missing_cols <- setdiff(required_cols, colnames(vc))
  if (length(missing_cols) > 0L) {
    stop("per_variant_qc is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  df <- vc
  df$high_collision <- df$P_ge2_given_present >= p_ge2_thresh

  # Ensure axes cover all observed values
  x_range <- range(df$n_samples_present, na.rm = TRUE)
  y_range <- range(df$P_ge2_given_present, na.rm = TRUE)

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = n_samples_present,
      y = P_ge2_given_present,
      color = high_collision
    )
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = p_ge2_thresh, linetype = "dashed") +
    ggplot2::scale_x_continuous(limits = x_range, expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(limits = y_range, expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::labs(
      x = "Number of samples with variant present",
      y = "P(K ≥ 2 | K ≥ 1)",
      title = "Variant collision risk",
      subtitle = paste0(
        "Variants in red have P_ge2_given_present ≥ ",
        signif(p_ge2_thresh, 3)
      ),
      color = "High collision\nrisk"
    ) +
    ggplot2::theme_bw()
}
