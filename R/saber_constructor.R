#' Row-wise normalized entropy (variant-level)
#'
#' Compute the normalized Shannon entropy of a numeric vector, typically
#' representing counts of a single variant across samples. Entropy is computed
#' on the non-zero entries only and normalized to lie in \eqn{[0, 1]} by
#' dividing by \eqn{\log(m)}, where \eqn{m} is the number of non-zero entries.
#'
#' @param x Numeric vector of counts for a single variant across samples.
#'
#' @details
#' This helper is intended for internal use when quantifying how evenly a
#' given barcode/variant is distributed across samples. Values near 0 indicate
#' that counts are concentrated in a single sample (clone-like), whereas values
#' near 1 indicate counts are spread evenly across many samples (stereotyped
#' motif-like).
#'
#' @return A single numeric scalar in \eqn{[0, 1]} giving the normalized entropy.
#'
#' @keywords internal
row_entropy_normalized <- function(x) {
  x <- as.numeric(x)
  x <- x[x > 0]
  if (length(x) <= 1L)
    return(0)
  p <- x / sum(x)
  H <- -sum(p * log(p))          # natural log
  H / log(length(p))             # normalize to [0, 1]
}

#' Column-wise normalized entropy (sample-level)
#'
#' Compute the normalized Shannon entropy of a numeric vector, typically
#' representing counts of all variants within a single sample. Entropy is
#' computed on the non-zero entries only and normalized to \eqn{[0, 1]}.
#'
#' @param x Numeric vector of counts for all variants in one sample.
#'
#' @details
#' This helper is used internally to summarize per-sample diversity of barcodes
#' before and after filtering stereotyped recurrent variants.
#'
#' @return A single numeric scalar in \eqn{[0, 1]} giving the normalized entropy.
#'
#' @keywords internal
col_entropy_normalized <- function(x) {
  x <- as.numeric(x)
  x <- x[x > 0]
  if (length(x) <= 1L)
    return(0)
  p <- x / sum(x)
  H <- -sum(p * log(p))
  H / log(length(p))
}

#' Estimate empirical count and VAF thresholds from noise
#'
#' Estimate count and variant allele frequency (VAF) thresholds for defining
#' barcode presence based on the empirical distribution of low-prevalence
#' variants (noise) across a count matrix.
#'
#' @param count_mat Numeric matrix of variant counts (rows = variants,
#'   columns = samples).
#' @param depths Optional numeric vector of per-sample depths (column sums).
#'   If \code{NULL}, computed as \code{colSums(count_mat)} after excluding
#'   rows matching \code{no_variant_pattern}.
#' @param no_variant_pattern Character scalar, regular expression matching
#'   rows that represent "no variant" background to be excluded from noise
#'   estimation.
#' @param background Character scalar specifying which cells to treat as
#'   background when estimating thresholds. Either \code{"singletons"} (only
#'   variants present in exactly one sample) or \code{"all"} (all non-zero
#'   entries).
#' @param noise_quantile Numeric scalar in \eqn{(0, 1)} giving the quantile
#'   of the empirical count and VAF distributions to treat as the upper bound
#'   of noise. Values just above this quantile are used as thresholds.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Excludes rows matching \code{no_variant_pattern}.
#'   \item Optionally restricts to "singleton" variants (present in exactly one
#'         sample) when \code{background = "singletons"}.
#'   \item Computes the specified \code{noise_quantile} of the empirical
#'         non-zero counts and VAFs and sets \code{count_thresh} to one plus
#'         the floored count quantile and \code{freq_thresh} to the VAF
#'         quantile.
#' }
#'
#' @return A list with elements:
#' \describe{
#'   \item{count_thresh}{Integer count threshold for presence.}
#'   \item{freq_thresh}{Numeric VAF threshold for presence.}
#'   \item{noise_quantile}{The noise quantile used.}
#'   \item{background}{The background mode used.}
#'   \item{n_nonzero_used}{Number of non-zero observations used to estimate
#'     thresholds.}
#' }
#'
#' @keywords internal
#' @importFrom stats quantile
estimate_barcode_thresholds <- function(count_mat,
                                        depths             = NULL,
                                        no_variant_pattern = "^no variant$",
                                        background         = c("singletons", "all"),
                                        noise_quantile     = 0.99) {
  if (!is.matrix(count_mat))
    count_mat <- as.matrix(count_mat)
  if (!is.numeric(count_mat))
    stop("count_mat must be numeric.")

  background <- match.arg(background)

  rn <- rownames(count_mat)
  if (!is.null(rn)) {
    no_var_rows <- grepl(no_variant_pattern, rn)
  } else {
    no_var_rows <- rep(FALSE, nrow(count_mat))
  }
  mat_use <- count_mat[!no_var_rows, , drop = FALSE]

  n_var  <- nrow(mat_use)
  n_samp <- ncol(mat_use)
  if (n_var < 1L ||
      n_samp < 1L)
    stop("Empty matrix after removing no-variant rows.")

  if (is.null(depths)) {
    depths <- colSums(mat_use)
  } else if (length(depths) != n_samp) {
    stop("length(depths) must equal ncol(count_mat).")
  }

  vaf_mat <- sweep(mat_use, 2, depths, "/")

  nz_mask    <- mat_use > 0
  counts_vec <- mat_use[nz_mask]
  vaf_vec    <- vaf_mat[nz_mask]

  if (background == "singletons") {
    nz_per_row <- rowSums(nz_mask)
    singleton_rows <- which(nz_per_row == 1L)

    if (length(singleton_rows) >= 50L) {
      nz_singleton <- nz_mask[singleton_rows, , drop = FALSE]
      counts_vec <- mat_use[singleton_rows, , drop = FALSE][nz_singleton]
      vaf_vec    <- vaf_mat[singleton_rows, , drop = FALSE][nz_singleton]
    }
  }

  counts_vec <- as.numeric(counts_vec)
  vaf_vec    <- as.numeric(vaf_vec)
  ok <- is.finite(vaf_vec) & !is.na(vaf_vec) & !is.na(counts_vec)
  counts_vec <- counts_vec[ok]
  vaf_vec    <- vaf_vec[ok]

  if (!length(counts_vec))
    stop("No non-zero counts available for noise estimation.")

  count_noise <- stats::quantile(counts_vec, probs = noise_quantile, names = FALSE)
  freq_noise  <- stats::quantile(vaf_vec, probs = noise_quantile, names = FALSE)

  count_thresh <- max(1L, as.integer(floor(count_noise)) + 1L)
  freq_thresh  <- as.numeric(freq_noise)

  list(
    count_thresh   = count_thresh,
    freq_thresh    = freq_thresh,
    noise_quantile = noise_quantile,
    background     = background,
    n_nonzero_used = length(counts_vec)
  )
}

#' Generate a logical presence matrix from counts and thresholds
#'
#' Create a logical matrix indicating presence/absence of variants in samples
#' based on count and VAF thresholds.
#'
#' @param count_mat Numeric matrix of variant counts (rows = variants,
#'   columns = samples).
#' @param count_thresh Integer minimum count threshold for presence.
#' @param freq_thresh Numeric minimum VAF threshold for presence.
#'
#' @details
#' Variant allele frequencies are computed as \code{count / depth}, where
#' \code{depth} is the column sum of \code{count_mat}. A cell is considered
#' present if its count is greater than or equal to \code{count_thresh} and
#' its VAF is greater than or equal to \code{freq_thresh}.
#'
#' @return A logical matrix of the same dimension as \code{count_mat} with
#'   \code{TRUE} indicating presence.
#'
#' @keywords internal
generate_present_mat <- function(count_mat, count_thresh, freq_thresh) {
  depths  <- colSums(count_mat)
  vaf_mat <- sweep(count_mat, 2, depths, "/")

  present_mat <- (count_mat >= count_thresh) &
    (vaf_mat >= freq_thresh)


}

#' Estimate per-variant barcode collision probabilities
#'
#' Estimate, for each variant, the probability that it arises multiple times
#' independently in a single sample under a simple Poisson model, using the
#' observed sharing of that variant across samples.
#'
#' @param count_mat Numeric matrix of variant counts (rows = variants,
#'   columns = samples).
#' @param count_thresh Integer minimum count threshold for presence, used to
#'   construct the presence matrix.
#' @param freq_thresh Numeric minimum VAF threshold for presence, used to
#'   construct the presence matrix.
#'
#' @details
#' A logical presence matrix is first constructed using
#' \code{generate_present_mat()}. For each variant \eqn{i}, the number of
#' samples in which it is present \eqn{m_i} is used to estimate
#' \eqn{\hat{p}_i = m_i / S}, where \eqn{S} is the number of samples. Under a
#' Poisson model for the number of independent origins per sample,
#' \eqn{K_i \sim \mathrm{Poisson}(\lambda_i)}, we have
#' \eqn{\hat{p}_i \approx P(K_i \ge 1) = 1 - e^{-\lambda_i}}, so
#' \eqn{\hat{\lambda}_i = -\log(1 - \hat{p}_i)}. From this,
#' \code{P_ge2_per_sample} and \code{P_ge2_given_present} are computed.
#'
#' @return A tibble with one row per variant containing:
#' \describe{
#'   \item{variant}{Variant identifier.}
#'   \item{n_samples_present}{Number of samples with the variant present.}
#'   \item{p_sample_present}{Estimated probability a random sample contains
#'     the variant at least once.}
#'   \item{lambda_per_sample}{Estimated Poisson mean per sample.}
#'   \item{P_ge2_per_sample}{Unconditional probability \eqn{P(K \ge 2)}.}
#'   \item{P_ge2_given_present}{Conditional probability \eqn{P(K \ge 2 \mid K \ge 1)}.}
#' }
#'
#' @keywords internal
#' @importFrom tibble tibble
estimate_barcode_collision_probs <- function(count_mat, count_thresh, freq_thresh) {
  present_mat <- generate_present_mat(count_mat, count_thresh, freq_thresh)

  n_var  <- nrow(present_mat)
  n_samp <- ncol(present_mat)

  var_ids <- rownames(present_mat)
  if (is.null(var_ids))
    var_ids <- paste0("var", seq_len(n_var))

  # Number of samples each variant appears in
  m <- rowSums(present_mat)

  # Empirical P(appears in a sample at least once)
  p_hat <- m / n_samp

  lambda <- rep(NA_real_, n_var)
  P_ge2  <- rep(NA_real_, n_var)
  P_ge2_given_ge1 <- rep(NA_real_, n_var)

  # Handle variants that never appear
  nonzero <- p_hat > 0 & p_hat < 1

  # For 0 or 1, handle separately to avoid -log(0) issues:
  # p_hat = 0 -> lambda = 0, P_ge2 = 0
  zero_idx <- which(p_hat == 0)
  if (length(zero_idx)) {
    lambda[zero_idx]         <- 0
    P_ge2[zero_idx]          <- 0
    P_ge2_given_ge1[zero_idx] <- NA_real_
  }

  # p_hat = 1 (barcode appears in all samples) -> lambda -> Inf under this model
  one_idx <- which(p_hat == 1)
  if (length(one_idx)) {
    lambda[one_idx]          <- Inf
    P_ge2[one_idx]           <- 1  # essentially guaranteed multiple occurrences
    P_ge2_given_ge1[one_idx] <- 1
  }

  # For 0 < p_hat < 1: normal case
  mid_idx <- which(p_hat > 0 & p_hat < 1)
  if (length(mid_idx)) {
    lambda[mid_idx] <- -log(1 - p_hat[mid_idx])

    # Unconditional P(K >= 2)
    P_ge2[mid_idx] <- 1 - (1 + lambda[mid_idx]) * exp(-lambda[mid_idx])

    # Conditional on being seen at least once
    P_ge1_mid <- 1 - exp(-lambda[mid_idx])
    P_ge2_given_ge1[mid_idx] <- P_ge2[mid_idx] / P_ge1_mid
  }

  tibble::tibble(
    variant            = var_ids,
    n_samples_present  = m,
    p_sample_present   = p_hat,
    lambda_per_sample  = lambda,
    P_ge2_per_sample   = P_ge2,
    P_ge2_given_present = P_ge2_given_ge1
  )
}

#' Construct a \code{Saber} object from CRISPR barcode data
#'
#' Build a \code{Saber} object from a \code{CrispRVariants} \code{crispr_set}
#' by estimating empirical noise thresholds, identifying stereotyped recurrent
#' barcode variants using a permutation-based test and entropy filter, and
#' summarizing per-sample and per-variant quality metrics.
#'
#' @param crispr_set A \code{CrispRVariants} \code{crisprSet} (or similar)
#'   object from which variant count matrices can be extracted via
#'   \code{CrispRVariants::variantCounts()}.
#' @param description Character scalar describing the analysis or dataset.
#' @param min_sample_depth Numeric scalar; samples with total depth (column
#'   sum) below this value are excluded prior to analysis.
#' @param background_count Integer count threshold used when computing
#'   per-variant collision probabilities (background model).
#' @param background_freq Numeric VAF threshold used when computing
#'   per-variant collision probabilities (background model).
#' @param depths Optional numeric vector of per-sample depths. If \code{NULL},
#'   computed as \code{colSums(count_mat)}.
#' @param no_variant_pattern Character scalar, regular expression used to
#'   identify rows corresponding to "no variant" background.
#' @param thresholds Optional list as returned by
#'   \code{estimate_barcode_thresholds()}. If \code{NULL}, thresholds are
#'   estimated from the data (recommended).
#' @param min_samples Integer; minimum number of samples in which a variant
#'   must be present (after thresholding) to be considered for recurrent
#'   filtering.
#' @param alpha Numeric scalar; FDR threshold applied to permutation-based
#'   empirical p-values when calling recurrent variants.
#' @param background Character scalar specifying how to define background when
#'   estimating noise thresholds; either \code{"singletons"} or \code{"all"}.
#' @param noise_quantile Numeric scalar in \eqn{(0, 1)}; quantile of empirical
#'   count and VAF distributions used as noise upper bound.
#' @param entropy_thresh Numeric scalar in \eqn{[0, 1]}; variants with
#'   normalized entropy greater than or equal to this value (i.e. more evenly
#'   distributed across samples) are more likely to be treated as stereotyped
#'   junk when also recurrent and significant in the permutation test.
#' @param B Integer; number of permutations used to estimate empirical
#'   p-values for variant recurrence.
#' @param seed Optional integer random seed passed to \code{set.seed()} for
#'   reproducibility of the permutation procedure.
#' @param progress Logical; if \code{TRUE}, progress messages are printed
#'   during the permutation loop.
#'
#' @details
#' The workflow implemented by \code{make.Saber()} is:
#' \enumerate{
#'   \item Extract a variant count matrix from \code{crispr_set} and drop
#'         samples with depth below \code{min_sample_depth}.
#'   \item Estimate empirical \code{count_thresh} and \code{freq_thresh} using
#'         \code{estimate_barcode_thresholds()} if \code{thresholds} is
#'         \code{NULL}.
#'   \item Define variant presence using these thresholds, and exclude
#'         \code{"no variant"} rows.
#'   \item For each variant, compute the number of samples in which it is
#'         present and its normalized entropy across samples.
#'   \item Perform a permutation-based recurrence test by shuffling the
#'         presence matrix within columns \code{B} times, deriving empirical
#'         p-values and FDR-adjusted p-values.
#'   \item Define stereotyped recurrent variants (\code{is_recurrent}) as those
#'         that (i) are not \code{"no variant"}, (ii) are present in at least
#'         \code{min_samples} samples, (iii) are significant at FDR
#'         \code{alpha}, and (iv) have entropy greater than or equal to
#'         \code{entropy_thresh}.
#'   \item Remove \code{"no variant"} rows and recurrent stereotyped variants
#'         to obtain \code{filtered_counts}.
#'   \item Summarize per-sample metrics (total reads, fractions in kept,
#'         recurrent, and no-variant reads, unique barcodes, entropy before
#'         and after filtering).
#'   \item Compute per-variant collision probabilities on the filtered matrix
#'         using \code{estimate_barcode_collision_probs()} with
#'         \code{background_count} and \code{background_freq}.
#'   \item Construct and return a \code{Saber} object with all of these
#'         components.
#' }
#'
#' @return A \code{Saber} object (S7/R7 class defined elsewhere) containing:
#' \itemize{
#'   \item \code{description}: Analysis description.
#'   \item \code{filtered_counts}: Matrix of counts with stereotyped junk and
#'         no-variant rows removed.
#'   \item \code{recurrent_variants}: Data frame of per-variant recurrence and
#'         entropy statistics.
#'   \item \code{motif_blacklist}: Character vector of variant IDs treated as
#'         uninformative (recurrent or no-variant).
#'   \item \code{per_sample_uniques}: Data frame of unique barcodes per sample.
#'   \item \code{per_sample_qc}: Data frame of per-sample QC metrics.
#'   \item \code{per_variant_qc}: Tibble of per-variant collision probabilities
#'         based on the background model.
#'   \item \code{thresholds}: List of thresholds and analysis parameters.
#' }
#'
#' @seealso \code{CrispRVariants::variantCounts()}, \code{estimate_barcode_thresholds()},
#'   \code{estimate_barcode_collision_probs()}
#'
#' @export
#' @importFrom CrispRVariants variantCounts
#' @importFrom stats p.adjust
#' @importFrom tibble tibble
make.Saber <- function(crispr_set,
                       description = "No description provided.",
                       min_sample_depth   = 1000,
                       background_count   = 100,
                       background_freq    = 0.001,
                       depths             = NULL,
                       no_variant_pattern = "^no variant$",
                       thresholds         = NULL,
                       min_samples        = 2L,
                       alpha              = 0.05,
                       background         = c("singletons", "all"),
                       noise_quantile     = 0.99,
                       entropy_thresh     = 0.8,
                       B                  = 1000L,
                       seed               = NULL,
                       progress           = TRUE) {
  count_mat <- CrispRVariants::variantCounts(crispr_set)
  count_mat <- count_mat[, colSums(count_mat) > min_sample_depth]

  background <- match.arg(background)

  rn <- rownames(count_mat)
  if (is.null(rn)) {
    rn <- paste0("var", seq_len(nrow(count_mat)))
    rownames(count_mat) <- rn
  }

  # Mark no-variant rows
  no_var_rows <- grepl(no_variant_pattern, rn)

  n_var  <- nrow(count_mat)
  n_samp <- ncol(count_mat)

  # Depths
  if (is.null(depths)) {
    depths <- colSums(count_mat)
  } else if (length(depths) != n_samp) {
    stop("length(depths) must equal ncol(count_mat).")
  }

  # Thresholds from data if needed
  if (is.null(thresholds)) {
    thresholds <- estimate_barcode_thresholds(
      count_mat,
      depths             = depths,
      no_variant_pattern = no_variant_pattern,
      background         = background,
      noise_quantile     = noise_quantile
    )
  }

  # VAF matrix
  vaf_mat <- sweep(count_mat, 2, depths, "/")

  # Presence definition (depth-aware)
  present <- (count_mat >= thresholds$count_thresh) &
    (vaf_mat    >= thresholds$freq_thresh)

  # Never treat "no variant" as present
  present[no_var_rows, ] <- FALSE

  # Observed recurrence
  k_obs <- rowSums(present)

  # Variant-level entropy across samples
  entropy <- apply(count_mat, 1, row_entropy_normalized)

  # ---- Permutation-based recurrence test ----
  if (!is.null(seed))
    set.seed(seed)

  exceed_counts <- integer(n_var)  # # perms where k_perm >= k_obs
  sum_perm_k    <- numeric(n_var)  # to estimate expected k under null

  if (progress) {
    message("Running ",
            B,
            " permutations on ",
            n_var,
            " variants x ",
            n_samp,
            " samples...")
  }

  for (b in seq_len(B)) {
    # Permute rows within each column, preserving per-sample prevalence
    perm_present <- present
    for (j in seq_len(n_samp)) {
      perm_present[, j] <- perm_present[sample.int(n_var), j]
    }

    k_perm <- rowSums(perm_present)

    exceed_counts <- exceed_counts + (k_perm >= k_obs)
    sum_perm_k    <- sum_perm_k + k_perm

    if (progress && (b %% max(1L, B %/% 10L) == 0L)) {
      message("  ... permutation ", b, "/", B)
    }
  }

  # Empirical p-values (+1 smoothing)
  p_emp <- (exceed_counts + 1) / (B + 1)
  p_adj <- stats::p.adjust(p_emp, method = "BH")

  expected_k <- sum_perm_k / B   # expected # samples under permuted null

  # Recurrent stereotyped junk:
  is_recurrent <- (!no_var_rows) &
    (k_obs >= min_samples) &
    (p_adj < alpha) &
    (entropy >= entropy_thresh)

  recurrent_tbl <- data.frame(
    variant       = rn,
    n_samples     = k_obs,
    expected_n    = expected_k,
    total_reads   = rowSums(count_mat),
    entropy       = entropy,
    p_value_emp   = p_emp,
    p_adj_emp     = p_adj,
    is_recurrent  = is_recurrent,
    is_no_variant = no_var_rows,
    stringsAsFactors = FALSE
  )

  # Matrix after removing no-variant + recurrent junk
  keep_rows    <- !no_var_rows & !is_recurrent
  filtered_mat <- count_mat[keep_rows, , drop = FALSE]

  # Per-sample unique informative barcodes
  unique_per_sample <- colSums(filtered_mat > 0)
  unique_tbl <- data.frame(
    sample          = colnames(count_mat),
    unique_barcodes = unique_per_sample,
    stringsAsFactors = FALSE
  )

  # Sample-level QC
  total_reads_all  <- colSums(count_mat)
  total_reads_keep <- colSums(filtered_mat)
  total_reads_junk <- colSums(count_mat[is_recurrent, , drop = FALSE])
  total_reads_no_var <- colSums(count_mat[no_var_rows, , drop = FALSE])

  sample_entropy_all  <- apply(count_mat, 2, col_entropy_normalized)
  sample_entropy_keep <- apply(filtered_mat, 2, col_entropy_normalized)

  per_sample_qc <- tibble(
    sample                 = colnames(count_mat),
    total_reads_all        = total_reads_all,
    total_reads_kept       = total_reads_keep,
    total_reads_junk       = total_reads_junk,
    total_reads_no_var     = total_reads_no_var,
    frac_reads_kept        = ifelse(
      total_reads_all > 0,
      total_reads_keep / total_reads_all,
      NA_real_
    ),
    frac_reads_junk        = ifelse(
      total_reads_all > 0,
      total_reads_junk / total_reads_all,
      NA_real_
    ),
    frac_reads_no_var      = ifelse(
      total_reads_all > 0,
      total_reads_no_var / total_reads_all,
      NA_real_
    ),
    unique_barcodes_kept   = unique_per_sample,
    entropy_all            = sample_entropy_all,
    entropy_after_filter   = sample_entropy_keep
  )

  # Motif/variant blacklist for upstream filtering
  motif_blacklist <- recurrent_tbl$variant[recurrent_tbl$is_recurrent |
                                             recurrent_tbl$is_no_variant]

  attr(filtered_mat, "thresholds") <- thresholds

  saber_args <- list(
    description         = description,
    filtered_counts     = filtered_mat,
    recurrent_variants  = recurrent_tbl[order(-recurrent_tbl$n_samples, recurrent_tbl$variant), ],
    motif_blacklist     = motif_blacklist,
    per_sample_uniques  = unique_tbl,
    per_sample_qc       = per_sample_qc,
    per_variant_qc      = estimate_barcode_collision_probs(
      filtered_mat,
      count_thresh = background_count,
      freq_thresh = background_freq
    ),
    thresholds          = list(
      thresholds,
      min_sample_depth   = min_sample_depth,
      background_count   = background_count,
      background_freq    = background_freq,
      alpha          = alpha,
      min_samples    = min_samples,
      entropy_thresh = entropy_thresh,
      B              = B
    )
  )
  Saber(
    description = saber_args$description,
    filtered_counts = saber_args$filtered_counts,
    recurrent_variants = saber_args$recurrent_variants,
    motif_blacklist = saber_args$motif_blacklist,
    per_sample_uniques = saber_args$per_sample_uniques,
    per_sample_qc = saber_args$per_sample_qc,
    per_variant_qc = saber_args$per_variant_qc,
    thresholds = saber_args$thresholds
  )
}
