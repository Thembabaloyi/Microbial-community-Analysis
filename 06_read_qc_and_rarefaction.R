# =============================================================================
# 06_read_qc_and_rarefaction.R
# Per-sample read statistics from the 36 CZ ID sample taxon reports
# (Table S3) and observed richness after rarefying species-level read counts
# to the smallest sample (Table S4).
#   - reads analysed      = nt_count / nt_rpm * 1e6 (CZ ID's rPM denominator)
#   - classified reads    = sum of species-level (tax_level 1) NT read counts
#   - thresholded species = >= 5 reads, NT rPM >= 1, alignment length >= 50 bp
# =============================================================================
source("R/00_setup.R")
o    <- load_ps()
meta <- o$meta

files <- list.files(paths$reports, pattern = "_taxon_report\\.csv(\\.gz)?$", full.names = TRUE)
if (length(files) != 36) stop("Expected 36 CZ ID taxon reports in ", paths$reports)
czid_id <- sub("^(\\d+(_\\d)?)_\\d+_taxon_report.*$", "\\1", basename(files))

read_report <- function(f) {
  d <- read.csv(f)
  num <- c("nt_count", "nt_rpm", "nt_alignment_length")
  d[num] <- lapply(d[num], function(x) { x <- suppressWarnings(as.numeric(x)); x[is.na(x)] <- 0; x })
  d
}
reports <- setNames(lapply(files, read_report), czid_id)

qc <- lapply(names(reports), function(id) {
  d   <- reports[[id]]
  sp  <- d[d$tax_level == 1 & d$nt_count > 0, ]
  ge  <- d[d$tax_level == 2, ]
  ref <- d[d$nt_count > 0 & d$nt_rpm > 0, ][1, ]
  total <- ref$nt_count / ref$nt_rpm * 1e6
  thr <- sp$nt_count >= 5 & sp$nt_rpm >= 1 & sp$nt_alignment_length >= 50
  cat_reads <- tapply(sp$nt_count, sp$category, sum)
  data.frame(czid_sample_id = id,
             reads_analysed = total,
             classified_reads_species = sum(sp$nt_count),
             pct_classified_species = 100 * sum(sp$nt_count) / total,
             n_species = nrow(sp), n_species_thresholded = sum(thr),
             n_genera = sum(ge$nt_count > 0),
             n_genera_thresholded = sum(ge$nt_count >= 5 & ge$nt_rpm >= 1),
             pct_bacteria  = 100 * sum(cat_reads["bacteria"], na.rm = TRUE) / sum(sp$nt_count),
             pct_eukaryota = 100 * sum(cat_reads["eukaryota"], na.rm = TRUE) / sum(sp$nt_count),
             pct_viruses   = 100 * sum(cat_reads["viruses"], na.rm = TRUE) / sum(sp$nt_count),
             pct_archaea   = 100 * sum(cat_reads["archaea"], na.rm = TRUE) / sum(sp$nt_count))
}) |> bind_rows()
qc <- left_join(meta[, c("czid_sample_id", "sample_id", "plastic_type", "rpm", "time")],
                qc, by = "czid_sample_id") |> arrange(sample_id)
save_table(qc, "TableS3_CZID_read_statistics.csv")

qc_summary <- data.frame(
  statistic = c("Reads analysed, min", "Reads analysed, max", "Reads analysed, CV (%)",
                "Species-level classified (%), mean", "min", "max",
                "Species per sample, min", "max",
                "Species per sample after thresholds, min", "max",
                "Genera after thresholds, min", "max",
                "Classified reads: bacteria (%)", "eukaryota (%)", "viruses (%)", "archaea (%)"),
  value = c(min(qc$reads_analysed), max(qc$reads_analysed),
            100 * sd(qc$reads_analysed) / mean(qc$reads_analysed),
            mean(qc$pct_classified_species), min(qc$pct_classified_species),
            max(qc$pct_classified_species), min(qc$n_species), max(qc$n_species),
            min(qc$n_species_thresholded), max(qc$n_species_thresholded),
            min(qc$n_genera_thresholded), max(qc$n_genera_thresholded),
            mean(qc$pct_bacteria), mean(qc$pct_eukaryota), mean(qc$pct_viruses),
            mean(qc$pct_archaea)))
print(qc_summary)
save_table(qc_summary, "TableS3_summary.csv")

# ---- Rarefied richness (binomial thinning to the smallest sample) -------------
# Rarefaction is random, so it is repeated N_RARE times; Table S4 reports the
# first draw per sample and the tests are summarised over all draws.
N_RARE <- 100
depth  <- min(qc$reads_analysed)
sp_counts <- lapply(qc$czid_sample_id, function(id) {
  d <- reports[[id]]; as.integer(d$nt_count[d$tax_level == 1 & d$nt_count > 0])
})
frac <- depth / qc$reads_analysed
set.seed(SEED)
rare_mat <- replicate(N_RARE, vapply(seq_along(sp_counts), function(i)
  sum(rbinom(length(sp_counts[[i]]), sp_counts[[i]], frac[i]) > 0), numeric(1)))

alpha <- read.csv(file.path(paths$tables, "Table_alpha_diversity_per_sample.csv"))
alpha <- alpha[match(qc$sample_id, alpha$sample_id), ]
rare  <- data.frame(sample_id = qc$sample_id, Plastic_type = alpha$Plastic_type,
                    Block = alpha$Block, reads_analysed = qc$reads_analysed,
                    Observed = alpha$Observed,
                    rarefied_richness = rare_mat[, 1],
                    rarefied_richness_mean = rowMeans(rare_mat))
save_table(rare, "TableS4_rarefied_richness.csv")

per_draw <- apply(rare_mat, 2, function(r) c(
  rho = cor(r, rare$Observed, method = "spearman"),
  KW  = kruskal.test(r ~ factor(rare$Plastic_type))$p.value,
  Fr  = friedman.test(r, factor(rare$Plastic_type), factor(rare$Block))$p.value))
summ <- function(v) sprintf("%.3f (range %.3f-%.3f)", median(v), min(v), max(v))
rare_tests <- data.frame(
  rarefaction_depth_reads = depth, n_rarefactions = N_RARE,
  spearman_rarefied_vs_unrarefied_median = summ(per_draw["rho", ]),
  KW_plastic_p_median = summ(per_draw["KW", ]),
  Friedman_plastic_blocks_p_median = summ(per_draw["Fr", ]),
  spearman_unrarefied_richness_vs_reads = cor(rare$Observed, rare$reads_analysed, method = "spearman"),
  spearman_unrarefied_richness_vs_pct_classified =
    cor(rare$Observed, qc$pct_classified_species, method = "spearman"))
print(t(rare_tests))
save_table(rare_tests, "TableS4_rarefaction_tests.csv")
cat("Tables S3 and S4 written.\n")
