# =============================================================================
# 02_taxonomic_composition.R
# Phylum- and genus-level composition (Results: Taxonomic composition)
#   Figure 1 - phylum composition per sample
#   Figure 2 - top-30 genera per sample
#   Tables   - mean relative abundances, per-plastic means, KW and Friedman tests
# Relative abundance = % of each sample's filtered rPM; reads not assigned at
# the rank are reported as "Unclassified".
# =============================================================================
source("R/00_setup.R")
o    <- load_ps()
meta <- o$meta

# ---- Relative abundance (%) per rank, by taxon name -------------------------
rel_by_rank <- function(rank) {
  m <- by_name(glom_matrix(o$ps_rel, rank)) * 100
  rbind(m, Unclassified = 100 - colSums(m))
}
phy <- rel_by_rank("Phylum")
gen <- rel_by_rank("Genus")

mean_tab <- function(m, rank) {
  data.frame(rank = rank, taxon = rownames(m),
             mean_rel_abund_pct = rowMeans(m), sd = apply(m, 1, stats::sd),
             row.names = NULL) |>
    arrange(desc(mean_rel_abund_pct))
}
phy_mean <- mean_tab(phy, "Phylum")
gen_mean <- mean_tab(gen, "Genus")
save_table(phy_mean, "Table_phylum_mean_relative_abundance.csv")
save_table(gen_mean, "Table_genus_mean_relative_abundance.csv")
save_table(data.frame(taxon = rownames(phy), phy, check.names = FALSE),
           "Table_phylum_relative_abundance_per_sample.csv")
save_table(data.frame(taxon = rownames(gen), gen, check.names = FALSE),
           "Table_genus_relative_abundance_per_sample.csv")

top5_phyla <- head(phy_mean$taxon[phy_mean$taxon != "Unclassified"], 5)
cat("Top 5 phyla together:", round(sum(phy_mean$mean_rel_abund_pct[
  match(top5_phyla, phy_mean$taxon)]), 1), "%\n")

# Planctomycetota at 150 RPM by day (reported in the text)
plancto <- data.frame(meta[, c("sample_id", "Plastic_type", "rpm", "Time")],
                      Planctomycetota = phy["Planctomycetota", rownames(meta)])
save_table(plancto |> group_by(rpm, Time) |>
             summarise(mean_pct = mean(Planctomycetota), .groups = "drop"),
           "Table_Planctomycetota_by_rpm_time.csv")

# ---- Genera compared among plastics -----------------------------------------
genus_tests <- lapply(c("Pseudomonas", "Flavobacterium", "Bradyrhizobium"), function(g) {
  v <- gen[g, rownames(meta)]
  data.frame(genus = g,
             PVC = mean(v[meta$Plastic_type == "PVC"]),
             HDPE = mean(v[meta$Plastic_type == "HDPE"]),
             PP = mean(v[meta$Plastic_type == "PP"]),
             kruskal_p = kruskal.test(v ~ meta$Plastic_type)$p.value,
             friedman_p = friedman.test(v, meta$Plastic_type, meta$Block)$p.value)
}) |> bind_rows()
genus_tests$friedman_p_BH <- p.adjust(genus_tests$friedman_p, "BH")
print(genus_tests)
save_table(genus_tests, "Table_genus_by_plastic_tests.csv")

# Co-variation of Bradyrhizobium with Planctomycetota
save_table(data.frame(
  comparison = "Bradyrhizobium vs Planctomycetota",
  spearman_rho = cor(gen["Bradyrhizobium", ], phy["Planctomycetota", ], method = "spearman")),
  "Table_Bradyrhizobium_Planctomycetota_spearman.csv")

# ---- Stacked bar helper (Figures 1 and 2) -----------------------------------
# Taxon names in italics, "Other" and "Unclassified" in plain text
italic_labels <- function(x) {
  parse(text = ifelse(x %in% c("Other", "Unclassified"),
                      paste0('"', x, '"'), paste0('italic("', x, '")')))
}
stack_plot <- function(m, keep, pal, legend_title, ncol_legend) {
  other <- colSums(m[setdiff(rownames(m), c(keep, "Unclassified")), , drop = FALSE])
  d <- rbind(m[keep, , drop = FALSE], Other = other,
             Unclassified = m["Unclassified", ])
  df <- as.data.frame(as.table(d)); names(df) <- c("Taxon", "sample", "pct")
  df <- left_join(df, data.frame(sample = rownames(meta), meta), by = "sample")
  lev <- c(keep, "Other", "Unclassified")
  df$Taxon <- factor(df$Taxon, levels = rev(lev))
  ggplot(df, aes(RPM, pct, fill = Taxon)) +
    geom_col(width = 0.8, colour = NA) +
    facet_grid(Time ~ Plastic_type) +
    scale_fill_manual(values = pal, breaks = lev, labels = italic_labels(lev)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
    labs(x = "Agitation speed (RPM)", y = "Relative abundance (%)",
         fill = legend_title) +
    guides(fill = guide_legend(ncol = ncol_legend, title.position = "top",
                               title.hjust = 0.5)) +
    theme_pub() +
    theme(legend.position = "bottom", panel.spacing = unit(1, "lines"),
          legend.text.align = 0)
}

# ---- Figure 1: phylum composition (10 most abundant phyla) ------------------
top_phyla <- head(phy_mean$taxon[phy_mean$taxon != "Unclassified"], 10)
top_phyla <- sort(top_phyla)
pal_phy <- c(Actinomycetota = "#E41A1C", Arthropoda = "#377EB8", Bacillota = "#4DAF4A",
             Bacteroidota = "#984EA3", Cyanobacteriota = "#FF7F00", Myxococcota = "#A65628",
             Planctomycetota = "#F781BF", Pseudomonadota = "#66C2A5",
             Streptophyta = "#FC8D62", Uroviricota = "#8DA0CB")
missing_col <- setdiff(top_phyla, names(pal_phy))
if (length(missing_col)) pal_phy[missing_col] <- hue_pal()(length(missing_col))
pal_phy <- c(pal_phy[top_phyla], Other = "#CCCCCC", Unclassified = "#808080")
fig1 <- stack_plot(phy, top_phyla, pal_phy, "Phylum", 4)
save_fig(fig1, "Figure1_phylum_composition", width = 9, height = 9.5)

# ---- Figure 2: top-30 genera -------------------------------------------------
top_genera <- head(gen_mean$taxon[gen_mean$taxon != "Unclassified"], TOP_N)
pal30 <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A",
           "#D62728", "#FF9896", "#9467BD", "#C5B0D5", "#8C564B", "#C49C94",
           "#E377C2", "#F7B6D2", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5",
           "#393B79", "#6B6ECF", "#637939", "#B5CF6B", "#8C6D31", "#E7BA52",
           "#843C39", "#D6616B", "#7B4173", "#CE6DBD", "#FDD835", "#E6C594")
pal_gen <- c(setNames(pal30[seq_along(top_genera)], top_genera),
             Other = "#CCCCCC", Unclassified = "#808080")
fig2 <- stack_plot(gen, top_genera, pal_gen, "Genus", 4)
save_fig(fig2, "Figure2_genus_composition_top30", width = 9, height = 11)

cat("Figures 1-2 and composition tables written.\n")
