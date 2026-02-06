## Multi-panel heatmap across SRPs (same gene order; per-SRP row Z-score)
## - panels: SRP301469, SRP345877, SRP369109, SRP484543, lab_data
## - gene order: up genes first, then down genes (no clustering)
## - columns: Untreated/Cnt left, Apoptotic cells/AC right (no SRR labels)
## - scaling done separately within each SRP (row-wise Z over samples in that SRP)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
})

srp_order <- c("SRP301469", "SRP345877", "SRP369109", "SRP484543", "lab_data")

# gene order (keep your tibble order)
genes_use <- c(unique(up_genes$gene), unique(down_genes$gene))

# colors for top annotation (adjust if your labels differ)
treat_cols <- c(
  "Untreated" = "#2C7FB8",      # blue
  "Cnt" = "#2C7FB8",
  "Apoptotic cells" = "#D7301F",# red
  "AC" = "#D7301F"
)

# z-score color scale (fixed)
col_fun <- circlize::colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))

make_one_heatmap <- function(se, genes_use, title, show_rows = FALSE) {
  # genes present
  genes_present <- genes_use[genes_use %in% rownames(se)]
  if (length(genes_present) == 0) stop(title, ": none of the requested genes found in rownames(se).")
  
  # log2(CPM+1)
  mat <- assay(se, "cpm")[genes_present, , drop = FALSE]
  mat <- log2(mat + 1)
  
  # Treatment column and ordering: control left, AC right
  cd <- as.data.frame(colData(se))
  if (!"Treatment" %in% colnames(cd)) stop(title, ": colData is missing 'Treatment'.")
  
  # map possible control labels to "control first"
  ctrl_first <- cd$Treatment %in% c("Untreated", "Cnt", "Control")
  ord_cols <- order(!ctrl_first, as.character(cd$Treatment))
  mat <- mat[, ord_cols, drop = FALSE]
  cd <- cd[ord_cols, , drop = FALSE]
  
  # per-SRP row-wise Z-score (separate scaling per SRP)
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0  # genes with 0 variance across samples
  
  # top annotation (Treatment only)
  ha <- HeatmapAnnotation(
    Treatment = cd$Treatment,
    col = list(Treatment = treat_cols),
    show_annotation_name = FALSE,
    annotation_legend_param = list(title = "Treatment")
  )
  
  Heatmap(
    z,
    name = "Z-score",
    col = col_fun,
    top_annotation = ha,
    column_title = title,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = show_rows,
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 5),
    heatmap_legend_param = list(at = c(-4, -2, 0, 2, 4))
  )
}

# build heatmaps, only last panel shows rownames
ht_list <- NULL
for (i in seq_along(srp_order)) {
  nm <- srp_order[i]
  se <- se_list[[nm]]
  if (is.null(se)) stop("Missing se_list[['", nm, "']]")
  
  ht <- make_one_heatmap(
    se = se,
    genes_use = genes_use,
    title = nm,
    show_rows = (i == length(srp_order))
  )
  
  ht_list <- if (is.null(ht_list)) ht else (ht_list + ht)
}

draw(
  ht_list,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
