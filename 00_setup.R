# =============================================================================
# 00_setup.R - packages, paths, shared settings and helper functions
# Sourced by every other script. Run scripts from the repository root.
# =============================================================================

required_pkgs <- c("biomformat", "phyloseq", "vegan", "permute", "ggplot2",
                   "pheatmap", "ggpubr", "dplyr", "tidyr", "patchwork",
                   "RColorBrewer", "scales")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                 logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall them first (see README, 'Installation').")
}

suppressPackageStartupMessages({
  library(biomformat); library(phyloseq); library(vegan); library(permute)
  library(ggplot2); library(pheatmap); library(ggpubr); library(dplyr)
  library(tidyr); library(patchwork); library(RColorBrewer); library(scales)
})

# ---- Paths (relative to the repository root) --------------------------------
paths <- list(
  biom      = "data/czid_combined_NT_rPM.biom",
  metadata  = "data/sample_metadata.csv",
  reports   = "data/czid_sample_taxon_reports",   # 36 CZ ID taxon reports
  processed = "results/processed",
  tables    = "results/tables",
  figures   = "results/figures"
)
for (d in c(paths$processed, paths$tables, paths$figures)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Analysis settings -------------------------------------------------------
N_PERM      <- 9999   # permutations for PERMANOVA / PERMDISP
SEED        <- 123    # random seed used before every permutation test
MIN_SAMPLES <- 2      # taxon must be present in >= 2 samples ...
MIN_TOTAL   <- 10     # ... and have a total rPM >= 10 across samples
TOP_N       <- 30     # genera shown in Figure 2 and Figure 3

# ---- Factor levels and colours (shared by all figures) ----------------------
PLASTIC_LEVELS <- c("PVC", "HDPE", "PP")
TIME_LEVELS    <- c("7 days", "14 days", "21 days")
RPM_LEVELS     <- c("0", "90", "125", "150")

COL_PLASTIC <- c(PVC = "#E41A1C", HDPE = "#377EB8", PP = "#4DAF4A")
COL_TIME    <- c("7 days" = "#FDBF6F", "14 days" = "#FF7F00", "21 days" = "#B15928")
COL_RPM     <- c("0" = "#EDF8E9", "90" = "#BAE4B3", "125" = "#74C476", "150" = "#238B45")
COL_RPM_PCOA <- c("0" = "#9ECAE1", "90" = "#4292C6", "125" = "#08306B", "150" = "#AE017E")

# ---- Helpers -----------------------------------------------------------------
# Write a data frame as CSV into results/tables
save_table <- function(df, name) {
  write.csv(df, file.path(paths$tables, name), row.names = FALSE)
  invisible(df)
}

# Save a ggplot as PDF, PNG and TIFF
save_fig <- function(plot, name, width, height, dpi = 600) {
  base <- file.path(paths$figures, name)
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(base, ".tiff"), plot, width = width, height = height, dpi = dpi,
         compression = "lzw")
}

# Tidy one row of an adonis2 table
adonis_row <- function(a, term, label = term) {
  data.frame(test = label, term = term, Df = a[term, "Df"],
             R2 = a[term, "R2"], F = a[term, "F"], p = a[term, "Pr(>F)"],
             row.names = NULL)
}

# Load the phyloseq objects created by 01_import_and_filter.R
load_ps <- function() {
  f <- file.path(paths$processed, "phyloseq_objects.rds")
  if (!file.exists(f)) stop("Run R/01_import_and_filter.R first.")
  readRDS(f)
}

# Agglomerate a phyloseq object to a taxonomic rank. Same result as
# phyloseq::tax_glom(ps, rank, NArm = TRUE) - taxa are merged when their full
# lineage down to `rank` is identical and taxa unassigned at `rank` are dropped -
# but much faster for ~12,000 taxa. Returns a matrix (lineages x samples) with
# the rank name in the attribute "taxon" and the lineage in rownames.
glom_matrix <- function(ps, rank) {
  tt    <- as(tax_table(ps), "matrix")
  m     <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) m <- t(m)
  ranks <- colnames(tt)[seq_len(match(rank, colnames(tt)))]
  keep  <- !is.na(tt[, rank])
  lin   <- apply(tt[keep, ranks, drop = FALSE], 1, paste, collapse = "|")
  agg   <- rowsum(m[keep, , drop = FALSE], lin)
  attr(agg, "taxon") <- sub("^.*\\|", "", rownames(agg))
  agg
}

# Sum glommed rows by taxon name (lineages with the same name are combined)
by_name <- function(gm) rowsum(gm, attr(gm, "taxon"))

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(colour = "black"),
          strip.background = element_rect(fill = "grey90", colour = NA),
          strip.text = element_text(face = "bold"),
          legend.text = element_text(size = base_size - 1))
}
