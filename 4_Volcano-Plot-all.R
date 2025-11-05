suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(ggplot2); library(tibble)
})

# combine to long df
volc_df <- imap_dfr(de_results, ~{
  as_tibble(.x) %>%
    transmute(gene, logFC, PValue,
              study = .y,
              neglog10p = -log10(pmax(PValue, 1e-300)),
              sig = PValue < 0.05 & abs(logFC) > 2)
})

# plot
ggplot(volc_df, aes(x = logFC, y = neglog10p, color = sig)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
  facet_wrap(~ study, ncol = 3, scales = "free_y") +
  labs(x = "log2 Fold Change", y = expression(-log[10](P)),
       title = "Volcano plots (P < 0.05 and |logFC| > 2 highlighted)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
