# ============================================================
# Resampling analysis:
# Human-Human vs Bot-Human factor matching distributions
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(psych)
  library(ggplot2)
})

if (!requireNamespace("clue", quietly = TRUE)) {
  install.packages("clue")
}
library(clue)

# -----------------------------
# file paths
# -----------------------------
humans_csv <- "C:/Users/alonafa.WISMAIN/Weizmann Institute Dropbox/Alona Faktor/post_doc/validity_paper/code/validity-EFA-main/validity-EFA-main/chemistry_humans.csv"
bots_csv   <- "C:/Users/alonafa.WISMAIN/Weizmann Institute Dropbox/Alona Faktor/post_doc/validity_paper/code/validity-EFA-main/validity-EFA-main/chemistry_chatbots.csv"

#humans_csv <- "C:/Users/alonafa.WISMAIN/Weizmann Institute Dropbox/Alona Faktor/post_doc/validity_paper/code/validity-EFA-main/validity-EFA-main/psychometric_q1_humans.csv"
#bots_csv   <- "C:/Users/alonafa.WISMAIN/Weizmann Institute Dropbox/Alona Faktor/post_doc/validity_paper/code/validity-EFA-main/validity-EFA-main/psychometric_q1_chatbots.csv"

output_dir <- "C:/Users/alonafa.WISMAIN/Weizmann Institute Dropbox/Alona Faktor/post_doc/validity_paper/code/Results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# settings
# -----------------------------
n_factors <- 3
sample_size <- 120
n_iter <- 100
rotate_method <- "oblimin"
fm_method <- "minres"
master_seed <- 20260308

# -----------------------------
# helpers
# -----------------------------
read_binary_noheader <- function(path, prefix = "Q") {
  df <- read_csv(path, col_names = FALSE, show_col_types = FALSE)
  names(df) <- paste0(prefix, seq_len(ncol(df)))
  df
}

coerce_binary <- function(df) {
  df %>%
    mutate(across(everything(), ~{
      v <- .x
      if (is.character(v)) v <- trimws(v)
      v <- suppressWarnings(as.numeric(v))
      v[!is.na(v) & !(v %in% c(0, 1))] <- NA_real_
      v
    }))
}

drop_zero_var <- function(df) {
  keep <- sapply(df, function(x) {
    u <- unique(x[!is.na(x)])
    length(u) >= 2
  })
  df[, keep, drop = FALSE]
}

tetra_cor_smooth <- function(df) {
  R <- psych::tetrachoric(df)$rho
  R <- psych::cor.smooth(R)
  diag(R) <- 1
  R
}

fit_efa_loadings <- function(df, nfactors, rotate, fm) {
  R <- tetra_cor_smooth(df)
  fa_obj <- psych::fa(
    r = R,
    nfactors = nfactors,
    n.obs = nrow(df),
    rotate = rotate,
    fm = fm
  )
  L <- as.matrix(unclass(fa_obj$loadings))
  rownames(L) <- colnames(df)
  colnames(L) <- paste0("F", seq_len(ncol(L)))
  L
}

hungarian_match_from_congruence <- function(congruence_mat) {
  nr <- nrow(congruence_mat)
  nc <- ncol(congruence_mat)
  
  cost_mat <- 1 - congruence_mat
  
  if (nr != nc) {
    n <- max(nr, nc)
    pad_value <- max(cost_mat, na.rm = TRUE) + 1
    cost_pad <- matrix(pad_value, nrow = n, ncol = n)
    cost_pad[1:nr, 1:nc] <- cost_mat
    assignment <- clue::solve_LSAP(cost_pad)
    matched_cols <- as.vector(assignment)[1:nr]
  } else {
    assignment <- clue::solve_LSAP(cost_mat)
    matched_cols <- as.vector(assignment)
  }
  
  matched_scores <- congruence_mat[cbind(seq_len(nr), matched_cols)]
  
  data.frame(
    Factor_1 = rownames(congruence_mat),
    Matched_Factor_2 = colnames(congruence_mat)[matched_cols],
    Congruence = matched_scores
  )
}

prepare_common_items <- function(df1, df2) {
  df1_nz <- drop_zero_var(df1)
  df2_nz <- drop_zero_var(df2)
  keep_items <- intersect(names(df1_nz), names(df2_nz))
  
  df1_keep <- df1[, keep_items, drop = FALSE]
  df2_keep <- df2[, keep_items, drop = FALSE]
  
  list(df1 = df1_keep, df2 = df2_keep, items = keep_items)
}

sample_rows <- function(df, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  idx <- sample(seq_len(nrow(df)), size = n, replace = FALSE)
  df[idx, , drop = FALSE]
}

# -----------------------------
# load original data
# -----------------------------
hum_all <- read_binary_noheader(humans_csv, prefix = "Q") %>% coerce_binary()
bot_all <- read_binary_noheader(bots_csv,   prefix = "Q") %>% coerce_binary()

# align number of columns
min_n <- min(ncol(hum_all), ncol(bot_all))
hum_all <- hum_all[, seq_len(min_n), drop = FALSE]
bot_all <- bot_all[, seq_len(min_n), drop = FALSE]

names(hum_all) <- paste0("Q", seq_len(min_n))
names(bot_all) <- paste0("Q", seq_len(min_n))

cat("Total humans:", nrow(hum_all), "\n")
cat("Total bots:  ", nrow(bot_all), "\n")
cat("Total items: ", ncol(hum_all), "\n")

if (nrow(hum_all) < sample_size) stop("Not enough human rows for the requested sample size.")
if (nrow(bot_all) < sample_size) stop("Not enough bot rows for the requested sample size.")

# -----------------------------
# containers
# -----------------------------
results <- vector("list", n_iter)

set.seed(master_seed)

# -----------------------------
# main loop
# -----------------------------
for (i in seq_len(n_iter)) {
  cat("Iteration", i, "of", n_iter, "\n")
  
  # Human sample 1 and Human sample 2
  hum_sampled_1 <- sample_rows(hum_all, sample_size)
  hum_sampled_2 <- sample_rows(hum_all, sample_size)
  
  hum_pair <- prepare_common_items(hum_sampled_1, hum_sampled_2)
  hum1_use <- hum_pair$df1
  hum2_use <- hum_pair$df2
  
  if (ncol(hum1_use) < 3) {
    warning(sprintf("Iteration %d skipped for human-human: too few common usable items.", i))
    next
  }
  
  L_hum_1 <- fit_efa_loadings(hum1_use, n_factors, rotate_method, fm_method)
  L_hum_2 <- fit_efa_loadings(hum2_use, n_factors, rotate_method, fm_method)
  
  cong_hh <- psych::factor.congruence(L_hum_1, L_hum_2)
  match_hh <- hungarian_match_from_congruence(cong_hh)
  mean_hh <- mean(match_hh$Congruence, na.rm = TRUE)
  
  # Bot sample vs same human sample 2
  bot_sampled <- sample_rows(bot_all, sample_size)
  
  hb_pair <- prepare_common_items(bot_sampled, hum_sampled_2)
  bot_use <- hb_pair$df1
  hum2b_use <- hb_pair$df2
  
  if (ncol(bot_use) < 3) {
    warning(sprintf("Iteration %d skipped for bot-human: too few common usable items.", i))
    next
  }
  
  L_bot <- fit_efa_loadings(bot_use, n_factors, rotate_method, fm_method)
  L_hum_2b <- fit_efa_loadings(hum2b_use, n_factors, rotate_method, fm_method)
  
  cong_bh <- psych::factor.congruence(L_hum_2b, L_bot)
  match_bh <- hungarian_match_from_congruence(cong_bh)
  mean_bh <- mean(match_bh$Congruence, na.rm = TRUE)
  
  results[[i]] <- list(
    iteration = i,
    mean_human_human = mean_hh,
    mean_bot_human = mean_bh,
    congruence_hh = cong_hh,
    congruence_bh = cong_bh,
    matching_hh = match_hh,
    matching_bh = match_bh
  )
}

# -----------------------------
# collect summary table
# -----------------------------
summary_df <- bind_rows(lapply(results, function(x) {
  if (is.null(x)) return(NULL)
  data.frame(
    iteration = x$iteration,
    Human_Human = x$mean_human_human,
    Bot_Human = x$mean_bot_human
  )
}))

if (nrow(summary_df) == 0) {
  stop("No successful iterations were completed.")
}

# long format for plotting
plot_df <- bind_rows(
  data.frame(
    iteration = summary_df$iteration,
    comparison = "Human-Human",
    mean_congruence = summary_df$Human_Human
  ),
  data.frame(
    iteration = summary_df$iteration,
    comparison = "Bot-Human",
    mean_congruence = summary_df$Bot_Human
  )
)

# -----------------------------
# save summary tables
# -----------------------------
summary_file <- file.path(output_dir, "resampling_factor_matching_summary.csv")
plot_file_csv <- file.path(output_dir, "resampling_factor_matching_long.csv")
rds_file <- file.path(output_dir, "resampling_factor_matching_full_results.rds")

write_csv(summary_df, summary_file)
write_csv(plot_df, plot_file_csv)
saveRDS(results, rds_file)

cat("\nSaved summary table to:\n", summary_file, "\n")
cat("Saved long plot table to:\n", plot_file_csv, "\n")
cat("Saved full results to:\n", rds_file, "\n")

# -----------------------------
# summary stats
# -----------------------------
stats_df <- plot_df %>%
  group_by(comparison) %>%
  summarise(
    n = n(),
    mean = mean(mean_congruence, na.rm = TRUE),
    sd = sd(mean_congruence, na.rm = TRUE),
    median = median(mean_congruence, na.rm = TRUE),
    min = min(mean_congruence, na.rm = TRUE),
    max = max(mean_congruence, na.rm = TRUE),
    .groups = "drop"
  )

print(stats_df)

stats_file <- file.path(output_dir, "resampling_factor_matching_stats.csv")
write_csv(stats_df, stats_file)

# -----------------------------
# plot distributions
# -----------------------------
p <- ggplot(plot_df, aes(x = mean_congruence, fill = comparison, color = comparison)) +
  geom_density(alpha = 0.25, adjust = 1.1) +
  labs(
    title = "Distribution of factor matching scores",
    subtitle = paste0(
      "EFA with ", n_factors, " factors; n = ", sample_size,
      " per sample; ", nrow(summary_df), " successful iterations"
    ),
    x = "Mean matched factor congruence",
    y = "Density"
  ) +
  theme_bw()

print(p)

plot_path <- file.path(output_dir, "resampling_factor_matching_distributions.png")
ggsave(plot_path, p, width = 9, height = 6, dpi = 300)

cat("Saved plot to:\n", plot_path, "\n")