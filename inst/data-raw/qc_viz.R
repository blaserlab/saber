top_n_any_column_scaled <- function(m, n, na.rm = TRUE) {
  if (!is.matrix(m)) {
    stop("`m` must be a matrix.")
  }
  if (n <= 0) {
    stop("`n` must be a positive integer.")
  }

  n <- as.integer(n)

  ## 1) Find rows that are in the top n values of any column
  idx_cols_list <- lapply(seq_len(ncol(m)), function(j) {
    col <- m[, j]
    o <- order(col, decreasing = TRUE, na.last = NA)
    if (length(o) == 0) return(integer(0))
    o[seq_len(min(n, length(o)))]
  })

  rows_keep <- sort(unique(unlist(idx_cols_list)))

  # if nothing qualifies, return empty matrix with same columns
  if (length(rows_keep) == 0L) {
    return(m[0, , drop = FALSE])
  }

  m_sub <- m[rows_keep, , drop = FALSE]

  ## 2) Column-wise 0–1 scaling
  mins   <- apply(m_sub, 2, min, na.rm = na.rm)
  maxs   <- apply(m_sub, 2, max, na.rm = na.rm)
  ranges <- maxs - mins

  # avoid division by zero for constant / degenerate columns
  ranges[ranges == 0 | is.na(ranges)] <- 1

  m_scaled <- sweep(m_sub, 2, mins, FUN = "-")
  m_scaled <- sweep(m_scaled, 2, ranges, FUN = "/")

  m_scaled
}

order <- colSums(res$filtered_counts) |> enframe() |> arrange(desc(value)) |> pull(name)
filtered_scaled_mat <- top_n_any_column_scaled(res$filtered_counts, n = 10)
filtered_scaled_sh <- SummarizedHeatmap(filtered_scaled_mat, colOrder = order)
bb_plot_heatmap_main(filtered_scaled_sh) + theme(axis.text.y = element_blank()) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
colSums(res$filtered_counts)

res$per_sample_qc |> View()


