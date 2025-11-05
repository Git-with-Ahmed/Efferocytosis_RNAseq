library(UpSetR)

alpha <- 0.05

# Build sets: upregulated (Apoptotic cells > Untreated)
sig_up_sets <- lapply(de_results, function(df) df$gene[df$FDR < alpha & df$logFC > 0])
sig_up_sets <- sig_up_sets[lengths(sig_up_sets) > 0]  # drop empty sets

# Build sets: downregulated (Apoptotic cells < Untreated)
sig_down_sets <- lapply(de_results, function(df) df$gene[df$FDR < alpha & df$logFC < 0])
sig_down_sets <- sig_down_sets[lengths(sig_down_sets) > 0]

# Plot 1: Upregulated only
if (length(sig_up_sets)) {
  inc_up <- fromList(sig_up_sets)
  upset(inc_up,
        nsets = length(sig_up_sets),
        nintersects = 30,
        order.by = "freq",
        keep.order = TRUE,
        mainbar.y.label = "Intersection size (Up)",
        sets.x.label = "Genes per dataset (Up)")
} else {
  message("No upregulated sets passed the filter.")
}

# Plot 2: Downregulated only
if (length(sig_down_sets)) {
  inc_down <- fromList(sig_down_sets)
  upset(inc_down,
        nsets = length(sig_down_sets),
        nintersects = 30,
        order.by = "freq",
        keep.order = TRUE,
        mainbar.y.label = "Intersection size (Down)",
        sets.x.label = "Genes per dataset (Down)")
} else {
  message("No downregulated sets passed the filter.")
}
