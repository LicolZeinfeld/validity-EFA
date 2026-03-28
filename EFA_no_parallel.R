# ============================================================
# EFA_Sunbok.R no parallel analysis
# ============================================================
# ============================================================
# EFA_Kaiser_Scree_Loadings.R
# Purpose:
#   1) Compute tetrachoric correlation matrices for humans and bots
#   2) Produce scree plots
#   3) Retain factors using Kaiser rule (eigenvalues > 1)
#   4) Fit EFA separately for humans and bots
#   5) Produce loading plots for both groups
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(psych)
})

# -----------------------------
# USER CONFIG
# -----------------------------
humans_csv <- "~/Desktop/CASEd/Reliability/Factor Analysis/chemistry_humans.csv"
bots_csv   <- "~/Desktop/CASEd/Reliability/Factor Analysis/chemistry_chatbots.csv"

rotate_method <- "oblimin"
fm_method     <- "minres"

# optional saving
save_plots <- TRUE
out_dir <- "~/Desktop/CASEd/Reliability/Factor Analysis/efa_kaiser_outputs"

# -----------------------------
# Helpers
# -----------------------------
read_binary_noheader <- function(path, prefix = "Q") {
  df <- read_csv(path, col_names = FALSE, show_col_types = FALSE)
  names(df) <- paste0(prefix, seq_len(ncol(df)))
  df <- as.data.frame(df)
  
  # force numeric 0/1 if possible
  df[] <- lapply(df, function(x) as.integer(as.character(x)))
  df
}

drop_zero_variance <- function(df) {
  keep <- sapply(df, function(x) var(x, na.rm = TRUE) > 0)
  df[, keep, drop = FALSE]
}

make_common_item_set <- function(df1, df2) {
  common <- intersect(names(df1), names(df2))
  list(
    df1 = df1[, common, drop = FALSE],
    df2 = df2[, common, drop = FALSE]
  )
}

get_tetra_R <- function(df) {
  psych::tetrachoric(df)$rho
}

get_kaiser_k <- function(R) {
  eig <- eigen(R)$values
  sum(eig > 1)
}

plot_scree_base <- function(R, title_text) {
  eig <- eigen(R)$values
  plot(
    seq_along(eig), eig,
    type = "b", pch = 19,
    xlab = "Component Number",
    ylab = "Eigenvalue",
    main = title_text
  )
  abline(h = 1, lty = 2)
}

plot_loadings_base <- function(fa_obj, title_text) {
  L <- as.matrix(unclass(fa_obj$loadings))
  nf <- ncol(L)
  
  matplot(
    L,
    type = "b",
    pch = 19,
    lty = 1,
    xlab = "Item Index",
    ylab = "Loading",
    main = title_text
  )
  abline(h = 0, lty = 2)
  legend(
    "topright",
    legend = colnames(L),
    col = seq_len(nf),
    lty = 1,
    pch = 19,
    bty = "n"
  )
}

# -----------------------------
# Load data
# -----------------------------
hum <- read_binary_noheader(humans_csv, prefix = "Q")
bot <- read_binary_noheader(bots_csv, prefix = "Q")

# manually drop selected items
items_to_drop <- c("Q8", "Q10", "Q12", "Q1", "Q17")

hum <- hum[, !(names(hum) %in% items_to_drop), drop = FALSE]
bot <- bot[, !(names(bot) %in% items_to_drop), drop = FALSE]

# remove zero-variance items separately first
hum <- drop_zero_variance(hum)
bot <- drop_zero_variance(bot)

# keep only common items across both groups
aligned <- make_common_item_set(hum, bot)
hum <- aligned$df1
bot <- aligned$df2

cat("Final #items used (common, nonzero variance):", ncol(hum), "\n\n")

# -----------------------------
# Tetrachoric correlations
# -----------------------------
R_h <- get_tetra_R(hum)
R_b <- get_tetra_R(bot)

# -----------------------------
# Kaiser factor retention
# -----------------------------
eig_h <- eigen(R_h)$values
eig_b <- eigen(R_b)$values

K_h <- sum(eig_h > 1)
K_b <- sum(eig_b > 1)

cat("Kaiser suggested K (humans):", K_h, "\n")
cat("Kaiser suggested K (bots):  ", K_b, "\n\n")

if (K_h < 1) K_h <- 1
if (K_b < 1) K_b <- 1

# -----------------------------
# Fit EFA separately
# -----------------------------
efa_h <- psych::fa(
  r = R_h,
  nfactors = K_h,
  n.obs = nrow(hum),
  rotate = rotate_method,
  fm = fm_method
)

efa_b <- psych::fa(
  r = R_b,
  nfactors = K_b,
  n.obs = nrow(bot),
  rotate = rotate_method,
  fm = fm_method
)

source("~/Desktop/CASEd/Reliability/Factor Analysis/Generate_Loading_Graph.R")
print(plot_comparison_loading_graph)

plot_comparison_loading_graph(
  efa_h = efa_h,
  efa_b = efa_b,
  title_text = "Psychometric Q1: Humans vs Bots Loading Graph",
  cutoff = 0.30,
  label_cutoff = 0.40,
  file = file.path(out_dir, "comparison_loading_graph.png")
)

# -----------------------------
# Print results
# -----------------------------
cat("===== HUMANS EFA =====\n")
print(efa_h)

cat("\n===== BOTS EFA =====\n")
print(efa_b)

# -----------------------------
# Output plots
# -----------------------------
if (save_plots) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Scree: humans
if (save_plots) {
  png(file.path(out_dir, "scree_humans.png"), width = 900, height = 700)
}
plot_scree_base(R_h, "Scree Plot - Humans")
if (save_plots) dev.off()

# Scree: bots
if (save_plots) {
  png(file.path(out_dir, "scree_bots.png"), width = 900, height = 700)
}
plot_scree_base(R_b, "Scree Plot - Bots")
if (save_plots) dev.off()

# Loading plot: humans
if (save_plots) {
  png(file.path(out_dir, "loadings_humans.png"), width = 1000, height = 700)
}
plot_loadings_base(efa_h, "Factor Loadings - Humans")
if (save_plots) dev.off()

# Loading plot: bots
if (save_plots) {
  png(file.path(out_dir, "loadings_bots.png"), width = 1000, height = 700)
}
plot_loadings_base(efa_b, "Factor Loadings - Bots")
if (save_plots) dev.off()

cat("\nPlots saved to:\n", out_dir, "\n")
