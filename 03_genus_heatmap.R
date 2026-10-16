# =============================================================================
# 03_genus_heatmap.R
# Figure 3 - heatmap of the 30 most abundant genera (log10 relative abundance,
# row-scaled to z-scores; rows clustered, samples ordered by plastic, time, RPM).
# Also: mean relative abundance of selected genera by RPM x time (text, Fig. 3).
# =============================================================================
source("R/00_setup.R")
o    <- load_ps()
meta <- o$meta

gen <- by_name(glom_matrix(o$ps_rel, "Genus")) * 100
top <- names(sort(rowMeans(gen), decreasing = TRUE))[seq_len(TOP_N)]

ord  <- with(meta, order(Plastic_type, Time, as.numeric(as.character(RPM))))
smp  <- rownames(meta)[ord]
PSEUDO <- 1e-3                         # pseudocount (%) before log10
mat  <- log10(gen[top, smp] + PSEUDO)

ann <- data.frame(`Plastic type` = meta[smp, "Plastic_type"],
                  `Exposure time` = meta[smp, "Time"],
                  `Agitation (RPM)` = meta[smp, "RPM"],
                  row.names = smp, check.names = FALSE)
ann_col <- list(`Plastic type` = COL_PLASTIC, `Exposure time` = COL_TIME,
                `Agitation (RPM)` = COL_RPM)
gaps <- cumsum(table(meta[smp, "Plastic_type"]))[1:2]

draw <- function(file = NA) {
  pheatmap(mat, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
           clustering_distance_rows = "euclidean", clustering_method = "complete",
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
           breaks = seq(-4, 4, length.out = 102),
           annotation_col = ann, annotation_colors = ann_col,
           gaps_col = gaps, fontsize = 10, fontsize_row = 10,
           labels_row = parse(text = paste0("italic('", rownames(mat), "')")),
           border_color = NA, angle_col = 90, filename = file,
           width = 10, height = 9)
}
for (ext in c("pdf", "png", "tiff")) {
  draw(file.path(paths$figures, paste0("Figure3_genus_heatmap_top30.", ext)))
}

# Selected genera by RPM x time (described in the text)
sel <- c("Cupriavidus", "Duganella", "Massilia", "Novosphingobium",
         "Bradyrhizobium", "Afipia", "Mesorhizobium", "Gemmata", "Frigoriglobus",
         "Escherichia", "Listeria")
sel_tab <- data.frame(meta[, c("rpm", "time")], t(gen[sel, rownames(meta)])) |>
  group_by(rpm, time) |> summarise(across(all_of(sel), mean), .groups = "drop")
save_table(sel_tab, "Table_selected_genera_by_rpm_time.csv")
cat("Figure 3 written.\n")
