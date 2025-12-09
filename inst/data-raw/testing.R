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


cross_result <- cross_sabers(list(DP = res_DP, DZ = res_DZ), n_top = 10, min_match_score = 5)
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
plot_saber_variant_collision(res_DZ, p_ge2_thresh = 0.05) |> View()

bind_rows(
sample_data(res_DP) |> select(sample, entropy_after_filter),
sample_data(res_DZ) |> select(sample, entropy_after_filter))

cross_result |> View()

cross_result |> select(true_id, sample_name, group) |> pivot_wider(names_from = group, values_from = sample_name, values_fill = NA) |> View()
