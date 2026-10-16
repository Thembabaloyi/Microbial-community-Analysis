# =============================================================================
# 01_import_and_filter.R
# Import the CZ ID BIOM file (NT rPM), rename samples, attach metadata and
# taxonomy, and filter low-abundance taxa (Methods: Statistical and Downstream
# Analysis). Output: results/processed/phyloseq_objects.rds
# =============================================================================
source("R/00_setup.R")

# ---- 1. BIOM file ----------------------------------------------------------
biom    <- read_biom(paths$biom)
otu_mat <- as(biom_data(biom), "matrix")

# CZ ID column ids look like "1_2:887036" -> CZ ID sample name before ":"
czid_names <- sub(":.*$", "", colnames(otu_mat))

# ---- 2. Sample metadata and renaming to 1-36 --------------------------------
meta <- read.csv(paths$metadata, colClasses = c(czid_sample_id = "character"))
stopifnot(setequal(czid_names, meta$czid_sample_id))
colnames(otu_mat) <- as.character(meta$sample_id[match(czid_names, meta$czid_sample_id)])
otu_mat <- otu_mat[, as.character(1:36)]

meta <- meta[order(meta$sample_id), ]
rownames(meta)    <- as.character(meta$sample_id)
meta$Plastic_type <- factor(meta$plastic_type, levels = PLASTIC_LEVELS)
meta$Time         <- factor(meta$time, levels = TIME_LEVELS)
meta$RPM          <- factor(as.character(meta$rpm), levels = RPM_LEVELS)
meta$Block        <- factor(paste0("R", meta$reactor, "_", meta$time))  # reactor x time

# ---- 3. Taxonomy ------------------------------------------------------------
# Row names hold the full CZ ID lineage ("Bacteria;;Pseudomonadota;...");
# observation metadata columns taxonomy1-7 = kingdom ... species.
# Empty strings mean the rank is not assigned (treated as "Unclassified").
tax_raw <- observation_metadata(biom)
tax_mat <- cbind(Domain = sub(";.*$", "", rownames(otu_mat)),
                 as.matrix(tax_raw[, paste0("taxonomy", 1:7)]))
colnames(tax_mat) <- c("Domain", "Kingdom", "Phylum", "Class", "Order",
                       "Family", "Genus", "Species")
tax_mat[tax_mat == ""] <- NA
rownames(tax_mat) <- rownames(otu_mat)

samp <- sample_data(meta)
sample_names(samp) <- rownames(meta)   # keep "1"-"36" (phyloseq would rename them)
ps_raw <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE),
                   tax_table(tax_mat), samp)

# ---- 4. Filtering -------------------------------------------------------------
ps      <- prune_taxa(taxa_sums(ps_raw) > 0, ps_raw)
ps_filt <- filter_taxa(ps, function(x) sum(x > 0) >= MIN_SAMPLES &
                                         sum(x) >= MIN_TOTAL, prune = TRUE)
# Relative abundance (proportion of each sample's filtered rPM)
ps_rel  <- transform_sample_counts(ps_filt, function(x) x / sum(x))
# Integer table for alpha diversity (rPM rounded to the nearest integer)
ps_int  <- ps_filt
otu_table(ps_int) <- otu_table(round(as(otu_table(ps_filt), "matrix")),
                               taxa_are_rows = TRUE)

# Number of genera / phyla = number of distinct lineages after agglomeration
# (as counted by phyloseq::tax_glom)
genera <- rownames(glom_matrix(ps_filt, "Genus"))
phyla  <- rownames(glom_matrix(ps_filt, "Phylum"))

summary_tab <- data.frame(
  step  = c("Taxa in BIOM file", "Taxa with non-zero abundance",
            "Taxa after filtering (>= 2 samples, total rPM >= 10)",
            "Unique genera after filtering", "Unique phyla after filtering"),
  value = c(ntaxa(ps_raw), ntaxa(ps), ntaxa(ps_filt), length(genera), length(phyla))
)
print(summary_tab)
save_table(summary_tab, "Table_filtering_summary.csv")

saveRDS(list(ps_raw = ps_raw, ps_filt = ps_filt, ps_rel = ps_rel,
             ps_int = ps_int, meta = meta),
        file.path(paths$processed, "phyloseq_objects.rds"))
cat("Saved results/processed/phyloseq_objects.rds\n")
