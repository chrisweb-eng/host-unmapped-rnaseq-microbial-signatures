
# ==========================================================
# MICROBIAL COMMUNITY ANALYSIS
# Rebuilt top-five Bracken abundance table
# ==========================================================

library(vegan)
library(ggplot2)

# ----------------------------------------------------------
# 1. CREATE OUTPUT DIRECTORY
# ----------------------------------------------------------

outdir <- "results/analysis_top5_rebuilt"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------
# 2. IMPORT ABUNDANCE TABLE
# ----------------------------------------------------------

abund <- read.csv(
  "results/abundance/top5_rebuilt/abundance_rebuilt.csv",
  row.names = 1,
  check.names = FALSE
)

cat("\nOriginal abundance table dimensions:\n")
print(dim(abund))

# Species are rows; samples are columns.

abund[] <- lapply(abund, as.numeric)

if (anyNA(abund)) {
  stop("Missing or non-numeric abundance values detected.")
}

if (any(as.matrix(abund) < 0)) {
  stop("Negative abundance values detected.")
}

# ----------------------------------------------------------
# 3. TRANSPOSE ABUNDANCE TABLE
# ----------------------------------------------------------

samples <- t(as.matrix(abund))

expected_samples <- c(
  paste0("SKIN_", 1:5),
  paste0("LUNG_", 1:5),
  paste0("COLON_", 1:5)
)

stopifnot(
  setequal(rownames(samples), expected_samples),
  nrow(samples) == 15
)

samples <- samples[
  expected_samples,
  ,
  drop = FALSE
]

if (any(rowSums(samples) == 0)) {
  stop("One or more samples have zero total abundance.")
}

cat("\nTransposed sample matrix dimensions:\n")
print(dim(samples))

cat("\nTotal abundance per sample:\n")
print(rowSums(samples))

write.csv(
  samples,
  file.path(outdir, "transposed_abundance.csv")
)

# ----------------------------------------------------------
# 4. BRAY-CURTIS DISSIMILARITY
# ----------------------------------------------------------

bray <- vegdist(
  samples,
  method = "bray"
)

bray_matrix <- as.matrix(bray)

write.csv(
  bray_matrix,
  file.path(outdir, "bray_curtis_matrix.csv")
)

# ----------------------------------------------------------
# 5. PRINCIPAL COORDINATES ANALYSIS
# ----------------------------------------------------------

pcoa <- cmdscale(
  bray,
  k = 2,
  eig = TRUE
)

pcoa_points <- as.data.frame(pcoa$points)

colnames(pcoa_points) <- c(
  "PCoA1",
  "PCoA2"
)

pcoa_points$Sample <- rownames(pcoa_points)

pcoa_points$Tissue <- sub(
  "_.*",
  "",
  pcoa_points$Sample
)

pcoa_points$Tissue <- factor(
  pcoa_points$Tissue,
  levels = c("COLON", "LUNG", "SKIN")
)

positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]

axis1 <- round(
  100 * pcoa$eig[1] / sum(positive_eigenvalues),
  2
)

axis2 <- round(
  100 * pcoa$eig[2] / sum(positive_eigenvalues),
  2
)

cat("\nPCoA variance explained:\n")
cat("PCoA1:", axis1, "%\n")
cat("PCoA2:", axis2, "%\n")

# ----------------------------------------------------------
# 6. GENERATE PCoA PLOT
# ----------------------------------------------------------

p <- ggplot(
  pcoa_points,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = Tissue
  )
) +
  geom_point(size = 4) +
  geom_text(
    aes(label = Sample),
    vjust = -0.7,
    size = 3
  ) +
  scale_colour_manual(
    values = c(
      "COLON" = "#E76F73",
      "LUNG" = "#28A745",
      "SKIN" = "#4285F4"
    )
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "PCoA of Microbial Communities",
    x = paste0("PCoA1 (", axis1, "%)"),
    y = paste0("PCoA2 (", axis2, "%)"),
    colour = "Tissue"
  )

ggsave(
  filename = file.path(outdir, "PCoA_rebuilt.png"),
  plot = p,
  width = 10,
  height = 7,
  dpi = 300
)

write.csv(
  pcoa_points,
  file.path(outdir, "PCoA_coordinates.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# 7. METADATA
# ----------------------------------------------------------

metadata <- data.frame(
  Tissue = factor(
    sub("_.*", "", rownames(samples)),
    levels = c("COLON", "LUNG", "SKIN")
  ),
  row.names = rownames(samples)
)

# ----------------------------------------------------------
# 8. OVERALL PERMANOVA
# ----------------------------------------------------------

set.seed(123)

permanova <- adonis2(
  bray ~ Tissue,
  data = metadata,
  permutations = 999
)

cat("\nOVERALL PERMANOVA:\n")
print(permanova)

capture.output(
  permanova,
  file = file.path(outdir, "PERMANOVA_overall.txt")
)

# ----------------------------------------------------------
# 9. PAIRWISE PERMANOVA
# ----------------------------------------------------------

pairwise_test <- function(group1, group2) {

  idx <- metadata$Tissue %in% c(group1, group2)

  bray_sub <- as.dist(
    as.matrix(bray)[idx, idx]
  )

  tissue_sub <- droplevels(
    metadata$Tissue[idx]
  )

  set.seed(123)

  adonis2(
    bray_sub ~ tissue_sub,
    permutations = 999
  )
}

skin_lung <- pairwise_test("SKIN", "LUNG")
skin_colon <- pairwise_test("SKIN", "COLON")
lung_colon <- pairwise_test("LUNG", "COLON")

cat("\nSKIN vs LUNG:\n")
print(skin_lung)

cat("\nSKIN vs COLON:\n")
print(skin_colon)

cat("\nLUNG vs COLON:\n")
print(lung_colon)

capture.output(
  {
    cat("SKIN vs LUNG\n")
    print(skin_lung)

    cat("\nSKIN vs COLON\n")
    print(skin_colon)

    cat("\nLUNG vs COLON\n")
    print(lung_colon)
  },
  file = file.path(outdir, "PERMANOVA_pairwise.txt")
)

# ----------------------------------------------------------
# 10. PERMDISP
# ----------------------------------------------------------

disp <- betadisper(
  bray,
  metadata$Tissue
)

disp_anova <- anova(disp)

set.seed(123)

disp_perm <- permutest(
  disp,
  permutations = 999
)

cat("\nPERMDISP ANOVA:\n")
print(disp_anova)

cat("\nPERMDISP PERMUTATION TEST:\n")
print(disp_perm)

capture.output(
  {
    print(disp_anova)
    print(disp_perm)
  },
  file = file.path(outdir, "PERMDISP_results.txt")
)

# ----------------------------------------------------------
# 11. SAVE ANALYSIS OBJECTS
# ----------------------------------------------------------

save(
  abund,
  samples,
  bray,
  pcoa,
  pcoa_points,
  metadata,
  permanova,
  skin_lung,
  skin_colon,
  lung_colon,
  disp,
  disp_anova,
  disp_perm,
  file = file.path(outdir, "microbial_analysis.RData")
)

# ----------------------------------------------------------
# 12. COMPLETION MESSAGE
# ----------------------------------------------------------

cat("\n========================================\n")
cat("MICROBIAL COMMUNITY ANALYSIS COMPLETED\n")
cat("========================================\n")

cat("\nSpecies analysed:", ncol(samples))
cat("\nBiological samples analysed:", nrow(samples))

cat("\nPCoA1:", axis1, "%")
cat("\nPCoA2:", axis2, "%")

cat("\n\nResults saved in:", outdir, "\n")

