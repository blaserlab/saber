# Variant entropy: how evenly a variant's counts are spread across samples
row_entropy_normalized <- function(x) {
  x <- as.numeric(x)
  x <- x[x > 0]
  if (length(x) <= 1L)
    return(0)
  p <- x / sum(x)
  H <- -sum(p * log(p))          # natural log
  H / log(length(p))             # normalize to [0, 1]
}

# Sample entropy: how evenly a sample's reads are spread across variants
col_entropy_normalized <- function(x) {
  x <- as.numeric(x)
  x <- x[x > 0]
  if (length(x) <= 1L)
    return(0)
  p <- x / sum(x)
  H <- -sum(p * log(p))
  H / log(length(p))
}

# Estimate count & VAF thresholds from empirical noise
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


generate_present_mat <- function(count_mat, count_thresh, freq_thresh) {
  depths  <- colSums(count_mat)
  vaf_mat <- sweep(count_mat, 2, depths, "/")

  present_mat <- (count_mat >= count_thresh) &
    (vaf_mat >= freq_thresh)


}

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
    P_ge2_given_present = P_ge2_given_ge1,
    stringsAsFactors   = FALSE
  )
}

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

  per_sample_qc <- data.frame(
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
    entropy_after_filter   = sample_entropy_keep,
    stringsAsFactors = FALSE
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
