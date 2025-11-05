# ---- Packages ----
library(tidyverse)
library(SummarizedExperiment)
library(S4Vectors)

setwd("/scratch/fa453/all_effero/Efferocytosis_RNAseq")
# ---- Read meta (for colData) ----
meta_all <- readr::read_csv(file.path("meta_all.csv"), show_col_types = FALSE) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(
    SRP       = as.factor(SRP),
    Treatment = as.factor(Treatment),
    Timepoint = as.factor(Timepoint)
  )

# ---- Files to include (exclude SRP311.., SRP4674.., SRP166533) ----
srp_files <- list.files(pattern = "^SRP\\d+_counts\\.csv$", full.names = TRUE)

# Add lab_data explicitly
lab_file <- file.path("lab_data_counts.csv")

# ---- Helper: build SE from a counts csv + id (SRP or 'lab_data') ----
make_se <- function(csv_path, meta_tbl, id = NULL) {
  if (is.null(id)) {
    id <- stringr::str_extract(basename(csv_path), "^SRP\\d+")
  }
  stopifnot(!is.na(id) || id == "lab_data")
  
  counts_df <- readr::read_csv(csv_path, show_col_types = FALSE)
  stopifnot("gene" %in% names(counts_df))
  
  mat <- counts_df %>%
    mutate(across(-gene, as.numeric)) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  cd <- meta_tbl %>%
    filter(SRP == id) %>%
    select(sample, SRP, Treatment, Timepoint)
  
  common <- intersect(colnames(mat), cd$sample)
  if (length(common) == 0L) stop(sprintf("No sample overlap for %s.", id))
  
  mat <- mat[, common, drop = FALSE]
  cd  <- cd %>% filter(sample %in% common) %>% arrange(match(sample, common))
  rownames(cd) <- cd$sample
  cd <- cd %>% select(-sample)
  
  se <- SummarizedExperiment(
    assays  = list(counts = mat),
    colData = S4Vectors::DataFrame(cd)
  )
  metadata(se) <- list(
    SRP = id,
    n_genes = nrow(mat),
    n_samples = ncol(mat),
    source_file = basename(csv_path)
  )
  se
}

# ---- Build SEs for SRPs ----
se_list <- map(srp_files, ~ make_se(.x, meta_tbl = meta_all))
names(se_list) <- str_extract(basename(srp_files), "^SRP\\d+")

# ---- Add lab_data SE ----
se_list[["lab_data"]] <- make_se(lab_file, meta_tbl = meta_all, id = "lab_data")