library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tibble)
library(purrr)
library(matrixStats)

# ---- Helper: compute PCA for one SE and return a data.frame + label ----
pca_df_one <- function(se, id, top_n = 1000L) {
  mat <- assay(se, "counts")
  x   <- log1p(mat)
  
  rv   <- matrixStats::rowVars(as.matrix(x))
  keep <- which(rv > 0)
  stopifnot(length(keep) > 1)          # need ≥2 genes for PCA
  ord  <- order(rv[keep], decreasing = TRUE)
  k    <- min(top_n, length(keep))
  x    <- x[keep[ord[seq_len(k)]], , drop = FALSE]
  
  pc  <- prcomp(t(x), center = TRUE, scale. = FALSE)
  var <- (pc$sdev^2) / sum(pc$sdev^2)
  
  df <- as_tibble(pc$x[, 1:2], rownames = "sample")
  colnames(df)[2:3] <- c("PC1", "PC2")
  
  ann <- rownames_to_column(as.data.frame(colData(se)), "sample") %>%
    select(sample, Treatment)
  
  df <- left_join(df, ann, by = "sample") %>%
    mutate(id = id)
  
  list(
    df = df,
    label = sprintf("%s (PC1 %.1f%%, PC2 %.1f%%)", id, 100 * var[1], 100 * var[2])
  )
}

# ---- Build one big DF for all experiments + facet labels ----
res_list <- imap(se_list, ~ pca_df_one(.x, .y, top_n = 1000L))
panel_df <- bind_rows(map(res_list, "df"))
facet_labels <- setNames(map_chr(res_list, "label"), map_chr(res_list, ~ unique(.x$df$id)))

# ---- Plot: one panel with facet per experiment ----
p_facets <- ggplot(panel_df, aes(PC1, PC2, color = Treatment)) +
  geom_point(size = 2.5) +
  facet_wrap(~ id, scales = "free", labeller = as_labeller(facet_labels)) +
  labs(
    title = "PCA (log1p counts) by experiment",
    subtitle = "Color = Treatment; each facet is an SRP or lab_data",
    x = "PC1", y = "PC2", color = "Treatment"
  ) +
  theme_light(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )
