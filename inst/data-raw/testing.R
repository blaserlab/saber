source("inst/data-raw/dependencies.R")
library(S7)
conflicted::conflicts_prefer(base::intersect)
conflicted::conflicts_prefer(dplyr::rename)
conflicted::conflicts_prefer(dplyr::count)
conflicted::conflicts_prefer(dplyr::filter)
devtools::load_all()

load("~/network/X/labs/Blaser/staff/ngs_archive/gestalt_DP/analysis_20251202155739785562/crispr_set.rda")

res_DP <- make.Saber(
  crispr_set         = crispr_set,
  no_variant_pattern = "^no variant$",
  background         = "singletons",
  noise_quantile     = 0.99,
  min_samples        = 2,
  alpha              = 0.05,
  entropy_thresh     = 0.8,
  B                  = 1000,   # bump up if you want smoother p-values
  seed               = 1,
  progress           = TRUE
)

res_DP@motif_blacklist
res_DP@per_variant_qc |> arrange(desc(P_ge2_given_present)) |>
  mutate(blacklist = case_when(variant %in% res_DP@motif_blacklist ~ TRUE, .default = FALSE))


load("~/network/X/labs/Blaser/staff/ngs_archive/gestalt_DZ/gestalt_DZ_analysis_20251126221935711982/crispr_set.rda")

res_DZ <- make.Saber(
  crispr_set         = crispr_set,
  no_variant_pattern = "^no variant$",
  background         = "singletons",
  noise_quantile     = 0.99,
  min_samples        = 2,
  alpha              = 0.05,
  entropy_thresh     = 0.8,
  B                  = 1000,   # bump up if you want smoother p-values
  seed               = 1,
  progress           = TRUE
)

new_rownames <- intersect(
rownames(barcodes(res_DP)),
rownames(barcodes(res_DZ)))
new_rownames
dim(barcodes(res_DP)[new_rownames,])
dim(barcodes(res_DZ)[new_rownames,])
dp_dz_cor0 <- cor(cbind(barcodes(res_DP)[new_rownames,],
          (barcodes(res_DZ)[new_rownames,])))
dp_dz_sh0 <- SummarizedHeatmap(dp_dz_cor0)
bb_plot_heatmap_main(dp_dz_sh0) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# dp_dz_dist <- dist(t(cbind(barcodes(res_DP)[new_rownames,],
#            (barcodes(res_DZ)[new_rownames,]))))
dp_dz_dist <- dist(t(scale(cbind(barcodes(res_DP)[new_rownames,],
                                 (barcodes(res_DZ)[new_rownames,])))))
fish_ids <- dp_dz_dist |>
  as.matrix() |>
  as_tibble(rownames = "sample") |>
  pivot_longer(-sample) |>
  filter(str_detect(sample, "DP")) |>
  filter(str_detect(name, "DZ")) |>
  group_by(name) |>
  slice_min(order_by = value, n = 1) |>
  group_by(sample) |>
  mutate(n = n()) |>
  filter(n == 1) |>
  rename(DP = sample, DZ = name) |>
  select(-value, -n) |>
  mutate(true_id = DP) |>
  relocate(true_id) |>
  ungroup()
fish_ids

list() |> length()

str(res_DP)

res_list <- list(t3mo = res_DP, t6mo = res_DZ)
res_list

tibble( = colnames(barcodes(res_DP)))



res <- cross_sabers(res_list)
res
res |> filter(!is.na(t6mo))

dp_dz_cor1 <- cor(cbind(barcodes(res_DP)[new_rownames,fish_ids$DP],
          (barcodes(res_DZ)[new_rownames,fish_ids$DZ])))

dp_dz_sh1 <- SummarizedHeatmap(dp_dz_cor1)
bb_plot_heatmap_main(dp_dz_sh1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

entropy_DP <- sample_data(res_DP) |>
  select(sample_name = sample, entropy_DP = entropy_after_filter, frk_DP = frac_reads_kept)

entropy_DZ <- sample_data(res_DZ) |>
  select(sample_name = sample, entropy_DZ = entropy_after_filter, frk_DZ = frac_reads_kept)

multiplex_barcode_key <- read_csv("~/network/X/Labs/Blaser/Brad/barcodes.csv", col_names = TRUE)

genotypes <- read_csv("~/network/X/Labs/Blaser/staff/ngs_archive/gestalt_DP/sample_metadata.csv") |> left_join(multiplex_barcode_key, by = join_by(multiplex_barcode == sequence)) |>
  mutate(sample_name = paste0(barcode, ".", sequencing_ID)) |>
  select(sample_name, hgfa_genotype)

entropy_data <- left_join(fish_ids, entropy_DP, by = join_by(DP == sample_name)) |>
  left_join(entropy_DZ, by = join_by(DZ == sample_name)) |>
  left_join(genotypes, by = join_by(DP == sample_name)) |>
  rename(DP_genotype = hgfa_genotype) |>
  left_join(genotypes, by = join_by(DZ == sample_name)) |>
  rename(DZ_genotype = hgfa_genotype) |>
  filter(DP_genotype == DZ_genotype) |>
  pivot_longer(c(entropy_DP, entropy_DZ), names_to = "timepoint", values_to = "entropy") |>
  filter(frk_DP > 0.25) |>
  filter(frk_DZ > 0.25) |>
  group_by(true_id) |>
  mutate(n = n()) |>
  filter(n == 2) |>
  ungroup()

ggplot(entropy_data |> filter(timepoint == "entropy_DP"), aes(x = frk_DP, y = entropy)) +
  geom_point()
ggplot(entropy_data |> filter(timepoint == "entropy_DZ"), aes(x = frk_DZ, y = entropy)) +
  geom_point()


ggplot(entropy_data, aes(x = DP_genotype, y = entropy)) +
  geom_jitter() +
  facet_wrap(~timepoint) +
  ggpubr::stat_compare_means(method = "t.test")

# Model: genotype effect, time effect, and their interaction


m <- lmerTest::lmer(entropy ~ DZ_genotype * timepoint + (1 | true_id), data = entropy_data)
summary(m)
anova(m)  # tests for genotype, time, interaction


library(emmeans)

emmeans(m, ~ DZ_genotype | timepoint)          # genotype at each timepoint
pvals <- contrast(emmeans(m, ~ DZ_genotype | timepoint), "pairwise") |> as_tibble()
pvals

devtools::load_all()
# x is a named list of Saber objects, ordered by timepoint
mapping <- cross_sabers(res_list)
mapping
# Simple presence/absence heatmap
plot_cross_sabers_heatmap(mapping)

# Same, but with sample IDs printed in tiles
plot_cross_sabers_heatmap(mapping, show_labels = TRUE)

devtools::load_all()
test <- saber_cor(res_list, mapping = mapping |> filter(!is.na(t6mo)))
test_sh <- SummarizedHeatmap(test)
bb_plot_heatmap_main(test_sh)

glimpse(test)
dim(test)

remove_na(c(mapping$t3mo, mapping$t6mo))
replace_na(c(mapping$t3mo, mapping$t6mo), replace = character(0))
