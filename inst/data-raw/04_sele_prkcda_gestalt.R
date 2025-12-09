source("00_packages_functions.R")

# data to load
gestalt_dataload <- read_csv("data_in/gestalt_dataload.csv")

# ingest the data
gestalt_data <- pmap(
  .l = list(
    filepath = gestalt_dataload$filepath,
    class = gestalt_dataload$class,
    exp = gestalt_dataload$exp,
    timepoint = gestalt_dataload$timepoint,
    dob = gestalt_dataload$dob,
    mpf = gestalt_dataload$mpf
  ),
  .f = function(filepath, class, exp, timepoint, dob, mpf) {
    data <- read_csv(filepath) %>%
      mutate(class = class,
             exp = exp,
             timepoint = timepoint,
             dob = dob,
             mpf = mpf)
  }
) %>%
  bind_rows() %>%
  filter(dob == 20180900)

# make a function to recursively rename the fish based on matching a set of barcodes at each timepoin (say 2 or more out of 5) to the earliest occurrence of that set of barcodes

rec_name_fish <- function(input) {
  # base case
  # generate a table of concatenated barcodes, the month you are working with and the name of the fish
  # the name of the fish here will become the first name for all of the subsequent samples
  fish_id_table <- input %>%
    filter(mpf == min(mpf)) %>%
    select(value,
           sample_name,
           starts_with("clone"),
           class,
           exp,
           timepoint,
           mpf) %>%
    mutate(barcode_cat = paste(clone1, clone2, clone3, clone4, clone5, sep = "_")) %>%
    select(barcode_cat, class, first_name = sample_name, mpf)

  # now take the same input data and calculate a score for barcode matching.
  # You get 1 point for matching each component of barcode_cat for a maximum of 5, since we have 5 clones we are matching for
  result <- input %>%
    select(value,
           sample_name,
           starts_with("clone"),
           class,
           exp,
           timepoint) %>%
    left_join(., fish_id_table) %>%
    mutate(
      clone1_score = str_detect(barcode_cat, clone1),
      clone2_score = str_detect(barcode_cat, clone2),
      clone3_score = str_detect(barcode_cat, clone3),
      clone4_score = str_detect(barcode_cat, clone4),
      clone5_score = str_detect(barcode_cat, clone5)
    ) %>%
    mutate(match_score = clone1_score + clone2_score + clone3_score + clone4_score + clone5_score) %>%
    group_by(sample_name) %>%
    filter(match_score == max(match_score)) %>%
    filter(match_score >= 2) # this is your hard-coded threshold setting for matching. 2 or more clones matching to the concatenated barcode will call a match.

  # now take the input data and get rid of all of the samples that successfully matched
  # by definition of the function, everything from the earliest timepoint gets a 5/5 score so will be removed
  # also by the last step of the function that generates result, everything with a score bigger than or equal to 2 will be removed
  newinput <- anti_join(input, result, by = "sample_name")

  if (nrow(newinput)==0) {
    # exit case:  this means that everything has been matched.
    # Anything that matches to only 1 timepoint, including the last timepoint is a potential outlier.
    # best to filter the result outside the function, but it could be filtered below.
    result
  } else {
    # recursive case:  when not everything has been matched yet.
    result <- bind_rows(result, rec_name_fish(input = newinput))
  }
}


gestalt_20180900_named <- rec_name_fish(input = gestalt_data)

gestalt_20180900_named_1 <- gestalt_20180900_named %>%
  ungroup() %>%
  select(first_name, timepoint, class, value) %>%
  mutate(timepoint = factor(timepoint, levels = c("3 mpf", "6 mpf", "9 mpf", "12 mpf", "22 mpf"))) %>%
  mutate(class = factor(class, levels = c("mcs", "prkcda")))


save.pigz(gestalt_20180900_named, file = "data/gestalt_20180900_named.rda")
save.pigz(gestalt_20180900_named_1, file = "data/gestalt_20180900_named_1.rda")

ggqqplot(gestalt_20180900_named_1, "value") + facet_grid(timepoint ~ class) # looks good

gestalt_20180900_named_1 %>%
  group_by(timepoint) %>%
  levene_test(value ~ class)

box_m(gestalt_20180900_named_1[, "value", drop = FALSE], gestalt_20180900_named_1$class)
bwtrim(value ~ timepoint*class, id = first_name, data = gestalt_20180900_named_1)
gestalt_20180900_named_1 %>%
  group_by(timepoint) %>%
  t_test(value ~ class, p.adjust.method = "bonferroni",detailed = T,var.equal = F)
gestalt_20180900_named_1 %>% group_by(timepoint, class) %>% get_summary_stats()
