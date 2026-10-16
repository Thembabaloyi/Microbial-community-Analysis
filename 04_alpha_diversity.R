# =============================================================================
# 04_alpha_diversity.R
# Observed richness, Shannon entropy and Gini-Simpson (1 - D) from rounded rPM
# (phyloseq::estimate_richness). Kruskal-Wallis by plastic type (Figure 4),
# Friedman tests with reactor x time blocks, and Kruskal-Wallis by exposure time.
# =============================================================================
source("R/00_setup.R")
o    <- load_ps()
meta <- o$meta

alpha <- estimate_richness(o$ps_int, measures = c("Observed", "Shannon", "Simpson"))
rownames(alpha) <- sub("^X", "", rownames(alpha))
alpha <- data.frame(meta[rownames(alpha), c("sample_id", "Plastic_type", "rpm",
                                            "Time", "Block")], alpha)
save_table(alpha, "Table_alpha_diversity_per_sample.csv")

metrics <- c("Observed", "Shannon", "Simpson")

# Mean +/- SD per plastic type
alpha_summary <- alpha |> group_by(Plastic_type) |>
  summarise(across(all_of(metrics), list(mean = base::mean, sd = stats::sd)), n = n())
print(alpha_summary)
save_table(alpha_summary, "Table_alpha_diversity_by_plastic.csv")

# Tests
alpha_tests <- lapply(metrics, function(v) data.frame(
  metric = v,
  KW_plastic_p  = kruskal.test(alpha[[v]] ~ alpha$Plastic_type)$p.value,
  Friedman_plastic_blocks_p = friedman.test(alpha[[v]], alpha$Plastic_type, alpha$Block)$p.value,
  KW_time_p     = kruskal.test(alpha[[v]] ~ alpha$Time)$p.value,
  KW_rpm_p      = kruskal.test(alpha[[v]] ~ factor(alpha$rpm))$p.value)) |> bind_rows()
print(alpha_tests)
save_table(alpha_tests, "Table_alpha_diversity_tests.csv")

# Eight samples with the lowest richness
save_table(head(alpha[order(alpha$Observed), c("sample_id", "Plastic_type", "rpm",
                                               "Time", "Observed")], 8),
           "Table_alpha_lowest_richness_samples.csv")

# ---- Figure 4 ------------------------------------------------------------------
labs_y <- c(Observed = "Observed richness", Shannon = "Shannon entropy",
            Simpson = "Simpson's index")
panel <- function(v, tag) {
  p <- alpha_tests$KW_plastic_p[alpha_tests$metric == v]
  ggplot(alpha, aes(Plastic_type, .data[[v]], fill = Plastic_type)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.75) +
    scale_fill_manual(values = COL_PLASTIC, guide = "none") +
    scale_y_continuous(labels = label_comma()) +
    labs(x = "Plastic type", y = labs_y[[v]], tag = tag,
         subtitle = sprintf("Kruskal-Wallis, p = %.2f", p)) +
    theme_pub() + theme(plot.subtitle = element_text(hjust = 0.5),
                        plot.tag = element_text(face = "bold"))
}
fig4 <- panel("Observed", "A") + panel("Shannon", "B") + panel("Simpson", "C")
save_fig(fig4, "Figure4_alpha_diversity_plastic", width = 11, height = 4)
cat("Figure 4 and alpha-diversity tables written.\n")
