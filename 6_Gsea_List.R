library(dplyr)
library(purrr)
library(tidyr)

# parameters
p_cut <- 0.05
studies <- names(de_results)                   # includes "lab_data"
stopifnot("lab_data" %in% studies)
meta_studies <- studies                        # use all 5
n_req <- 3                                     # significant in >=4 studies

# long table of all studies
long <- imap_dfr(de_results, ~{
  df <- as.data.frame(.x)
  df %>% transmute(study = .y, gene, logFC, PValue)
})

# per-gene summaries across studies
gene_summ <- long %>%
  group_by(gene) %>%
  summarize(
    n_studies = n(),
    n_sig = sum(PValue < p_cut, na.rm = TRUE),
    all_up  = all(logFC > 0, na.rm = TRUE),
    all_down= all(logFC < 0, na.rm = TRUE),
    .groups = "drop"
  )

# keep genes that are directionally consistent AND significant in >=4 studies
keepers <- gene_summ %>%
  filter(n_sig >= n_req, (all_up | all_down)) %>%
  mutate(direction = ifelse(all_up, "up","down")) %>%
  select(gene, direction)

# pull lab_data logFC for ranking
lab_fc <- de_results$lab_data %>%
  select(gene, lab_logFC = logFC)

rank_tbl <- keepers %>%
  inner_join(lab_fc, by = "gene") %>%         # require gene present in lab_data
  mutate(rank_score = lab_logFC) %>%
  arrange(desc(rank_score))

# outputs
up_genes   <- rank_tbl %>% filter(direction=="up")   %>% arrange(desc(rank_score))
down_genes <- rank_tbl %>% filter(direction=="down") %>% arrange(rank_score)
