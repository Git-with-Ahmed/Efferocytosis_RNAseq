library("tidyr")

add_cpm <- function(se) {
  y <- DGEList(counts = assay(se, "counts"))
  y <- calcNormFactors(y)
  assay(se, "cpm") <- cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
  se
}

se_list <- map(se_list, add_cpm) 


# keep only apoptotic-cell samples from your metadata
apo_meta <- meta_all %>%
  dplyr::filter(Treatment == "Apoptotic cells") %>%
  dplyr::mutate(SRP = as.character(SRP))

# build CPM long table only for those samples
cpm_long <- purrr::imap_dfr(se_list, ~ {
  study <- as.character(.y)
  stopifnot("cpm" %in% assayNames(.x))
  keep_samps <- apo_meta %>% dplyr::filter(SRP == study) %>% dplyr::pull(sample)
  if (length(keep_samps) == 0) return(tibble())                  # nothing to plot for this study
  se_sub <- .x[, colnames(.x) %in% keep_samps, drop = FALSE]
  if (ncol(se_sub) == 0) return(tibble())
  as_tibble(assay(se_sub, "cpm"), rownames = "gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "CPM") %>%
    mutate(study = study)
}) %>%
  mutate(logCPM1 = log10(CPM + 1))

# plot panels by study
ggplot(cpm_long, aes(logCPM1)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  facet_wrap(~ study, scales = "free_y") +
  labs(title = "Apoptotic cells only: log10(CPM+1)",
       x = "log10(CPM+1)", y = "Frequency") +
  theme_classic(base_size = 12)



# genes with median log10(CPM+1) ≥ 1 per study (Apoptotic cells only, from cpm_long)
thr <- 0.5
genes_by_study <- cpm_long |>
  dplyr::group_by(study, gene) |>
  dplyr::summarise(med_logCPM1 = median(logCPM1), .groups = "drop_last") |>
  dplyr::filter(med_logCPM1 >= thr) |>
  dplyr::arrange(study, desc(med_logCPM1))

# named list
genes_lists <- genes_by_study |>
  dplyr::group_split(study, .keep = FALSE) |>
  purrr::set_names(unique(genes_by_study$study)) |>
  purrr::map(~ .x$gene)

# union and ≥k-studies consensus
genes_union <- sort(unique(unlist(genes_lists)))
k <- 5
genes_kplus <- names(which(table(unlist(
  purrr::imap(genes_lists, ~ setNames(.x, rep(.y, length(.x))))
)) >= k))

