source("inst/data-raw/dependencies.R")
library(S7)
devtools::load_all()

load("~/network/X/labs/Blaser/staff/ngs_archive/gestalt_DP/analysis_20251202155739785562/crispr_set.rda")

res_DP <- make.Saber(
  crispr_set         = crispr_set,
  no_variant_pattern = "^no variant$",
  background         = "singletons",
  noise_quantile     = 0.99,
  min_samples        = 2,
  alpha              = 0.05,
  entropy_thresh     = 0.2,
  B                  = 1000,   # bump up if you want smoother p-values
  seed               = 1,
  progress           = TRUE
)

load("~/network/X/labs/Blaser/staff/ngs_archive/gestalt_DZ/gestalt_DZ_analysis_20251126221935711982/crispr_set.rda")

res_DZ <- make.Saber(
  crispr_set         = crispr_set,
  no_variant_pattern = "^no variant$",
  background         = "singletons",
  noise_quantile     = 0.99,
  min_samples        = 2,
  alpha              = 0.05,
  entropy_thresh     = 0.2,
  B                  = 1000,   # bump up if you want smoother p-values
  seed               = 1,
  progress           = TRUE
)

devtools::load_all()

cross_result <- cross_sabers(list(DP = res_DP, DZ = res_DZ), n_top = 10, min_match_score = 4)
View(cross_result)
cross_result |> count(group, sample_name) |> arrange(desc(n))
cross_result |> group_by(sample_name) |> mutate(n = n()) |> filter(n > 1) |> View()

sample_data(res_DP) |> glimpse()
variant_data(res_DP) |> glimpse()
ggplot(sample_data(res_DP), aes(x = sample, y = entropy_after_filter)) +
  geom_bar(stat = ("identity"))


devtools::load_all()

# Variant-level diagnostics
plot_saber_variant_entropy(res_DP)

plot_saber_sample_qc(res_DP)

# Variant-level collision risk, flagging >1% or >5%
plot_saber_variant_collision(res_DP, p_ge2_thresh = 0.05) |> View()
plot_saber_variant_collision(res_DZ, p_ge2_thresh = 0.05)

variant_data(res_DZ) |> View()

bind_rows(
sample_data(res_DP) |> select(sample, entropy_after_filter),
sample_data(res_DZ) |> select(sample, entropy_after_filter))

cross_result |> View()


cross_result |> count(group, true_id) |> arrange(desc(n))


cross_result |> select(group, true_id, sample_name) |> View()
cross_result |> count(true_id) |> arrange(desc(n)) |> View()

fish_ids <- cross_result |>
  select(group, true_id, sample_name) |>
  group_by(group, true_id) |>
  mutate(n = n()) |>
  filter(n == 1) |>
  ungroup() |>
  select(-n) |>
  pivot_wider(names_from = group, values_from = sample_name)
View(fish_ids)

entropy_DP <- sample_data(res_DP) |>
  select(sample_name = sample, entropy_DP = entropy_after_filter)

entropy_DZ <- sample_data(res_DZ) |>
  select(sample_name = sample, entropy_DZ = entropy_after_filter)

multiplex_barcode_key <- read_csv("~/network/X/Labs/Blaser/Brad/barcodes.csv", col_names = TRUE)

genotypes <- read_csv("~/network/X/Labs/Blaser/Brad/bleeds tracking.csv") |>
    select(DP = `bleed 1 (3mo) sequencing ID`, genotype, multiplex_barcode = `multiplex barcode`) |> left_join(multiplex_barcode_key, by = join_by(multiplex_barcode == sequence)) |> mutate(sample_name = paste0(barcode, ".", DP)) |>
  select(sample_name, genotype)
genotypes

entropy_data <- left_join(fish_ids, entropy_DP, by = join_by(DP == sample_name)) |>
  left_join(entropy_DZ, by = join_by(DZ == sample_name)) |>
  left_join(genotypes, by = join_by(true_id == sample_name)) |>
  filter(!is.na(genotype)) |>
  select(genotype, entropy_DP, entropy_DZ) |>
  pivot_longer(c(entropy_DP, entropy_DZ))
entropy_data

ggplot(entropy_data, aes(x = genotype, y = value)) +
  geom_jitter() +
  facet_wrap(~name) +
  ggpubr::stat_compare_means(method = "t.test")

