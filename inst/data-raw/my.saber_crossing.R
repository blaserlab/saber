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

  tibble::tibble(sample_name = col_names,
                 as.data.frame(res_mat, stringsAsFactors = FALSE))
}

rec_name <- function(input) {
  group_order <- input |>
    pull(group) |>
    unique()
  id_table <- input |>
    dplyr::filter(group == group_order[1]) |>
    dplyr::mutate(barcode_cat = paste0("_", top_1, "_", top_2, "_", top_3, "_", top_4, "_", top_5, "_")) |>
    dplyr::select(first_name = sample_name, barcode_cat)

  result <- dplyr::cross_join(input, id_table) |>
    dplyr::mutate(
      clone1_score = stringr::str_detect(barcode_cat, paste0("_", top_1, "_")),
      clone2_score = stringr::str_detect(barcode_cat, paste0("_", top_2, "_")),
      clone3_score = stringr::str_detect(barcode_cat, paste0("_", top_3, "_")),
      clone4_score = stringr::str_detect(barcode_cat, paste0("_", top_4, "_")),
      clone5_score = stringr::str_detect(barcode_cat, paste0("_", top_5, "_")),
    ) |>
    dplyr::mutate(match_score = clone1_score + clone2_score + clone3_score + clone4_score + clone5_score) |>
    dplyr::group_by(sample_name) |>
    dplyr::filter(match_score == max(match_score)) |>
    dplyr::filter(match_score >= 2) # this is your hard-coded threshold setting for matching. 2 or more clones matching to the concatenated barcode will call a match.
  newinput <- anti_join(input, result, by = "sample_name")

  if (nrow(newinput) == 0) {
    # exit case:  this means that everything has been matched.
    # Anything that matches to only 1 timepoint, including the last timepoint is a potential outlier.
    # best to filter the result outside the function, but it could be filtered below.
    result
  } else {
    # recursive case:  when not everything has been matched yet.
    result <- bind_rows(result, rec_name(input = newinput))
  }

}

cross_sabers <- function(saber_list) {
  input <- map2(.x = saber_list,
                .y = names(saber_list),
                .f = \(x, y) {
                  top_n_rownames_by_column(barcodes(x), n = 5) |>
                    mutate(group = y) |>
                    relocate(group)
                }) |> bind_rows()

  rec_name(input) |>
    select(group,
           sample_name,
           true_id = first_name,
           barcode_cat,
           starts_with("top"),
           match_score)


}
