library(edgeR)
library(SummarizedExperiment)
library(purrr)
library(tibble)

run_edgeR_two_group <- function(se) {
  counts <- assay(se, "counts")
  # normalize Treatment labels and set desired order
  group <- factor(gsub("_", " ", as.character(colData(se)$Treatment)),
                  levels = c("Untreated", "Apoptotic cells"))
  stopifnot(nlevels(group) == 2)
  
  # DGEList, filter, TMM
  y <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(y, group = group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  
  # Design, dispersion, fit
  design <- model.matrix(~ 0 + group)
  colnames(design) <- make.names(levels(group))  # "Untreated", "Apoptotic.cells"
  y    <- estimateDisp(y, design)
  fit  <- glmQLFit(y, design)
  
  # Contrast: Apoptotic cells - Untreated
  contr <- makeContrasts(Apoptotic.cells - Untreated, levels = design)
  qlf   <- glmQLFTest(fit, contrast = contr)
  
  # Full DE table
  topTags(qlf, n = Inf)$table |>
    rownames_to_column("gene")
}

# ---- run across all datasets ----
de_results <- purrr::imap(se_list, ~ run_edgeR_two_group(.x))

# Example: peek
# head(de_results$SRP301469)
