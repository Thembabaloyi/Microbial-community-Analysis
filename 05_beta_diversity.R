# =============================================================================
# 05_beta_diversity.R
# Bray-Curtis dissimilarity on filtered rPM, PCoA (Figure 5), PERMANOVA
# (vegan::adonis2), pairwise PERMANOVA (BH), PERMDISP (vegan::betadisper) and
# the sensitivity analyses (without the 150 RPM reactor; without the inoculated
# genera Escherichia and Listeria; on condition means).
# =============================================================================
source("R/00_setup.R")
o    <- load_ps()
meta <- o$meta

otu <- t(as(otu_table(o$ps_filt), "matrix"))          # samples x taxa
otu <- otu[rownames(meta), ]
bc  <- vegdist(otu, method = "bray")

perm_free  <- function() how(nperm = N_PERM)
perm_block <- function(b) how(blocks = b, nperm = N_PERM)
run_adonis <- function(formula, data, permutations = perm_free(), by = "terms") {
  set.seed(SEED)
  adonis2(formula, data = data, permutations = permutations, by = by)
}

res <- list()
# ---- 1. One-way PERMANOVA (all 36 samples) ---------------------------------
for (f in c("Plastic_type", "Time", "RPM")) {
  a <- run_adonis(as.formula(paste("bc ~", f)), meta)
  res[[length(res) + 1]] <- adonis_row(a, f, paste("One-way, all samples:", f))
}
# ---- 2. Plastic type with reactor x time block in the model -----------------
a <- run_adonis(bc ~ Block + Plastic_type, meta, perm_block(meta$Block))
res[[length(res) + 1]] <- adonis_row(a, "Plastic_type",
                                     "Plastic type | block term, permutations within blocks")
# ---- 3. Three-factor model, marginal tests ------------------------------------
a <- run_adonis(bc ~ Time + RPM + Plastic_type, meta, by = "margin")
for (f in c("Time", "RPM", "Plastic_type"))
  res[[length(res) + 1]] <- adonis_row(a, f, paste("Marginal, Time + RPM + Plastic:", f))

# ---- 4. Sensitivity: run 1 only (0, 90, 125 RPM; n = 27) ----------------------
k   <- meta$rpm != 150
m1  <- droplevels(meta[k, ])
bc1 <- vegdist(otu[k, ], method = "bray")
for (f in c("Time", "RPM", "Plastic_type")) {
  a <- run_adonis(as.formula(paste("bc1 ~", f)), m1)
  res[[length(res) + 1]] <- adonis_row(a, f, paste("Run 1 (no 150 RPM), one-way:", f))
}
# Sequential model: RPM after Time, Time x RPM interaction, then plastic type
a <- run_adonis(bc1 ~ Time + RPM + Plastic_type, m1, by = "margin")
res[[length(res) + 1]] <- adonis_row(a, "RPM", "Run 1: RPM after Time (marginal, Time + RPM + Plastic)")
a <- run_adonis(bc1 ~ Time * RPM + Plastic_type, m1)
res[[length(res) + 1]] <- adonis_row(a, "Time:RPM", "Run 1: Time x RPM interaction (Time * RPM + Plastic)")
a <- run_adonis(bc1 ~ Block + Plastic_type, m1, perm_block(m1$Block))
res[[length(res) + 1]] <- adonis_row(a, "Plastic_type",
                                     "Run 1: Plastic type | block term, permutations within blocks")
# Condition means (mean rPM of the three plastics per reactor x time; n = 9)
cond   <- m1$Block
otu_cm <- rowsum(otu[k, ], cond) / as.vector(table(cond)[rownames(rowsum(otu[k, ], cond))])
m_cm   <- m1[match(rownames(otu_cm), m1$Block), c("Time", "RPM")]
a <- run_adonis(vegdist(otu_cm, "bray") ~ Time + RPM, m_cm, by = "margin")
for (f in c("Time", "RPM"))
  res[[length(res) + 1]] <- adonis_row(a, f, paste("Run 1 condition means (n = 9), marginal:", f))

# ---- 5. Sensitivity: without inoculated genera Escherichia and Listeria -------
genus <- as.character(tax_table(o$ps_filt)[colnames(otu), "Genus"])
keep_tax <- !(genus %in% c("Escherichia", "Listeria"))
bc_ns <- vegdist(otu[, keep_tax], method = "bray")
for (f in c("Plastic_type", "Time", "RPM")) {
  a <- run_adonis(as.formula(paste("bc_ns ~", f)), meta)
  res[[length(res) + 1]] <- adonis_row(a, f, paste("Without Escherichia/Listeria:", f))
}

permanova <- bind_rows(res)
permanova$permutations <- N_PERM
print(permanova, digits = 3)
save_table(permanova, "Table_PERMANOVA.csv")

# ---- 6. Pairwise PERMANOVA (BH-adjusted) ---------------------------------------
pairwise_adonis <- function(d, groups, label) {
  dm <- as.matrix(d); lv <- levels(droplevels(groups))
  out <- combn(lv, 2, simplify = FALSE) |> lapply(function(pr) {
    s  <- groups %in% pr
    g  <- droplevels(groups[s])
    set.seed(SEED)
    a  <- adonis2(as.dist(dm[s, s]) ~ g, permutations = N_PERM)
    data.frame(factor = label, group1 = pr[1], group2 = pr[2],
               R2 = a$R2[1], F = a$F[1], p = a$`Pr(>F)`[1])
  }) |> bind_rows()
  out$p_BH <- p.adjust(out$p, "BH")
  out
}
pairwise <- bind_rows(pairwise_adonis(bc, meta$Time, "Time"),
                      pairwise_adonis(bc, meta$RPM, "RPM"),
                      pairwise_adonis(bc, meta$Plastic_type, "Plastic_type"))
print(pairwise, digits = 3)
save_table(pairwise, "Table_PERMANOVA_pairwise.csv")

# ---- 7. PERMDISP -----------------------------------------------------------------
disp_test <- function(d, groups, label) {
  bd <- betadisper(d, groups, type = "median", bias.adjust = TRUE)
  set.seed(SEED)
  pt <- permutest(bd, permutations = N_PERM)
  means <- tapply(bd$distances, groups, mean)
  data.frame(data = label, factor = deparse(substitute(groups)),
             F = pt$tab$F[1], p = pt$tab$`Pr(>F)`[1],
             mean_distance_to_median = paste(names(means), round(means, 3),
                                             sep = ": ", collapse = "; "))
}
permdisp <- bind_rows(
  disp_test(bc, meta$Plastic_type, "All samples"),
  disp_test(bc, meta$Time, "All samples"),
  disp_test(bc, meta$RPM, "All samples"),
  disp_test(bc1, m1$Plastic_type, "Run 1 (no 150 RPM)"),
  disp_test(bc1, m1$Time, "Run 1 (no 150 RPM)"),
  disp_test(bc1, m1$RPM, "Run 1 (no 150 RPM)"))
permdisp$factor <- sub("^.*\\$", "", permdisp$factor)
print(permdisp, digits = 3)
save_table(permdisp, "Table_PERMDISP.csv")

# ---- 8. Bray-Curtis within vs between reactor x time blocks ---------------------
dm <- as.matrix(bc); same <- outer(meta$Block, meta$Block, "==")
lt <- lower.tri(dm)
save_table(data.frame(comparison = c("Within reactor x time block", "Between blocks"),
                      mean_bray_curtis = c(mean(dm[lt & same]), mean(dm[lt & !same]))),
           "Table_BrayCurtis_within_between_blocks.csv")
write.csv(dm, file.path(paths$tables, "Table_BrayCurtis_matrix.csv"))

# ---- 9. PCoA and Figure 5 -------------------------------------------------------
pc  <- ordinate(o$ps_filt, method = "PCoA", distance = "bray")
ev  <- pc$values$Relative_eig * 100
pts <- data.frame(meta, Axis.1 = pc$vectors[rownames(meta), 1],
                  Axis.2 = pc$vectors[rownames(meta), 2])
save_table(data.frame(axis = 1:5, variance_explained_pct = ev[1:5]),
           "Table_PCoA_variance_explained.csv")
save_table(pts[, c("sample_id", "Plastic_type", "rpm", "Time", "Axis.1", "Axis.2")],
           "Table_PCoA_coordinates.csv")

p_lab <- function(f) {
  r <- permanova[permanova$test == paste("One-way, all samples:", f), ]
  p <- ifelse(r$p < 0.001, '"< 0.001"', sprintf('"= %.3f"', r$p))
  sprintf('R^2 == "%.3f" * "," ~ p ~ %s', r$R2, p)
}
pcoa_panel <- function(f, pal, title, tag) {
  ggplot(pts, aes(Axis.1, Axis.2, colour = .data[[f]])) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    stat_ellipse(type = "t", level = 0.95, linetype = "dashed", linewidth = 0.6,
                 show.legend = FALSE) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = pal) +
    annotate("text", x = Inf, y = -Inf, label = "PERMANOVA", hjust = 1.05,
             vjust = -2.2, size = 3, colour = "grey25") +
    annotate("text", x = Inf, y = -Inf, label = p_lab(f), parse = TRUE,
             hjust = 1.05, vjust = -0.4, size = 3, colour = "grey25") +
    labs(x = sprintf("PCoA1 (%.1f%%)", ev[1]), y = sprintf("PCoA2 (%.1f%%)", ev[2]),
         colour = NULL, title = title, tag = tag) +
    theme_pub() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5),
          plot.tag = element_text(face = "bold"))
}
fig5 <- pcoa_panel("Plastic_type", COL_PLASTIC, "Plastic type", "a") +
        pcoa_panel("Time", COL_TIME, "Exposure time", "b") +
        pcoa_panel("RPM", COL_RPM_PCOA, "Agitation speed (RPM)", "c")
save_fig(fig5, "Figure5_PCoA_BrayCurtis", width = 13, height = 5)
cat("Figure 5 and beta-diversity tables written.\n")
