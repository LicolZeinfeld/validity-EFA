# ============================================================
# EFA_Sunbok.R
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(psych)
})

# -----------------------------
# file paths
# -----------------------------
humans_csv <- "~/Desktop/CASEd/Reliability/Factor Analysis/psychometric_q2_humans.csv"
bots_csv   <- "~/Desktop/CASEd/Reliability/Factor Analysis/psychometric_q2_chatbots.csv"


pairs <- list()

drop_items <- c()

rotate_method  <- "oblimin"
fm_method      <- "minres"
loading_cutoff <- 0.30

use_parallel_analysis <- TRUE
pa_n_iter <- 30

alpha <- 0.05

# Optional: save results to RDS
save_results_rds <- TRUE
results_rds_path <- "~/Desktop/CASEd/Reliability/Factor Analysis/efa_chemistry_results.rds"
save_plots <- TRUE
out_dir <- "~/Desktop/CASEd/Reliability/Factor Analysis/efa_kaiser_outputs"

# ============================================================
# Helpers
# ============================================================

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

combine_pairs_strict <- function(df, pairs) {
  if (length(pairs) == 0) return(df)
  for (p in pairs) {
    stopifnot(length(p) >= 2, all(p %in% names(df)))
    new_name <- paste0(paste(p, collapse = "_AND_"))
    df[[new_name]] <- as.integer(rowSums(df[, p, drop = FALSE], na.rm = FALSE) == length(p))
    df <- df[, setdiff(names(df), p), drop = FALSE]
  }
  df
}

tetra_cor_smooth <- function(df) {
  R <- psych::tetrachoric(df)$rho
  R <- psych::cor.smooth(R)
  diag(R) <- 1
  R <- R + diag(1e-8, nrow(R))
  diag(R) <- 1
  R
}

kaiser_nfactors <- function(R) {
  ev <- eigen(R, only.values = TRUE)$values
  nf <- sum(ev > 1)
  if (nf < 1) nf <- 1
  list(eigenvalues = ev, nfactors = nf)
}

run_fa <- function(R, n_obs, nfactors, rotate, fm) {
  psych::fa(r = R, nfactors = nfactors, n.obs = n_obs, rotate = rotate, fm = fm)
}

has_heywood <- function(fa_obj, tol = 1e-6) {
  h2 <- fa_obj$communality
  u2 <- fa_obj$uniquenesses
  any(is.na(h2)) || any(is.na(u2)) || any(h2 > 1 + tol) || any(u2 < -tol)
}

get_loadings_matrix <- function(fa_obj) {
  L <- unclass(fa_obj$loadings)
  if (is.null(dim(L))) L <- matrix(L, ncol = fa_obj$factors)
  L
}

# Step down until admissible; logs explicitly when it rejects a K
fit_group_solution <- function(R, n_obs, k_start, rotate, fm, label) {
  k_start <- max(1, as.integer(k_start))
  for (k in seq(k_start, 1, by = -1)) {
    fa_k <- run_fa(R, n_obs, k, rotate, fm)
    hey  <- has_heywood(fa_k)
    cat(sprintf("  [%s] Try K=%d | Heywood/NA: %s\n", label, k, hey))
    if (!hey) {
      if (k < k_start) {
        cat(sprintf("  [%s] NOTE: Requested K=%d rejected as inadmissible; using largest admissible K=%d.\n",
                    label, k_start, k))
      }
      return(list(K = k, fa = fa_k, rejected_from = k_start))
    }
  }
  stop(sprintf("[%s] All K produced Heywood/NA issues. Consider fewer items / more N / different settings.", label))
}

# ============================================================
# Load & Prepare (align items so groups match)
# ============================================================

hum_raw <- read_binary_noheader(humans_csv, prefix = "Q") %>% coerce_binary()
bot_raw <- read_binary_noheader(bots_csv,   prefix = "Q") %>% coerce_binary()

min_n <- min(ncol(hum_raw), ncol(bot_raw))
hum <- hum_raw[, seq_len(min_n), drop = FALSE]
bot <- bot_raw[, seq_len(min_n), drop = FALSE]
names(hum) <- paste0("Q", seq_len(min_n))
names(bot) <- paste0("Q", seq_len(min_n))

hum <- combine_pairs_strict(hum, pairs)
bot <- combine_pairs_strict(bot, pairs)

if (length(drop_items) > 0) {
  hum <- hum[, setdiff(names(hum), drop_items), drop = FALSE]
  bot <- bot[, setdiff(names(bot), drop_items), drop = FALSE]
}

hum_nz <- drop_zero_var(hum)
bot_nz <- drop_zero_var(bot)
keep_items <- intersect(names(hum_nz), names(bot_nz))

hum <- hum[, keep_items, drop = FALSE]
bot <- bot[, keep_items, drop = FALSE]

cat("Final #items used (common, nonzero variance):", ncol(hum), "\n")
if (ncol(hum) < 3) stop("Too few usable items after filtering.")

# ============================================================
# Factor count suggestion
# ============================================================

R_h <- tetra_cor_smooth(hum)
R_b <- tetra_cor_smooth(bot)

kai_h <- kaiser_nfactors(R_h)
kai_b <- kaiser_nfactors(R_b)

cat("Kaiser suggested K (humans):", kai_h$nfactors, "\n")
cat("Kaiser suggested K (bots):  ", kai_b$nfactors, "\n")

K_h_start <- kai_h$nfactors
K_b_start <- kai_b$nfactors

if (use_parallel_analysis) {
  pa_h <- psych::fa.parallel(
    hum,
    fm = fm_method,
    fa = "fa",
    cor = "tet",
    n.iter = pa_n_iter,
    plot = FALSE
  )
  
  pa_b <- psych::fa.parallel(
    bot,
    fm = fm_method,
    fa = "fa",
    cor = "tet",
    n.iter = pa_n_iter,
    plot = FALSE
  )
  
  if (!is.null(pa_h$nfact) && !is.na(pa_h$nfact) && pa_h$nfact >= 1) K_h_start <- pa_h$nfact
  if (!is.null(pa_b$nfact) && !is.na(pa_b$nfact) && pa_b$nfact >= 1) K_b_start <- pa_b$nfact
  
  cat("Parallel Analysis suggested K (humans):", K_h_start, "\n")
  cat("Parallel Analysis suggested K (bots):  ", K_b_start, "\n")
}

# ============================================================
# Fit EFA (separate groups; step-down only if inadmissible)
# ============================================================

# ============================================================
# Fit EFA (separate groups; step-down only if inadmissible)
# ============================================================

# ============================================================
# Fit EFA exactly at the selected K (NO step-down)
# ============================================================

cat("\n=== PRIMARY: Separate EFA fits (NO step-down) ===\n")
cat("Humans: fitting K =", K_h_start, "\n")
fa_h <- run_fa(R_h, nrow(hum), K_h_start, rotate_method, fm_method)

cat("Bots:   fitting K =", K_b_start, "\n")
fa_b <- run_fa(R_b, nrow(bot), K_b_start, rotate_method, fm_method)

K_h_final <- K_h_start
K_b_final <- K_b_start

cat("\n[FINAL separate K]\n")
cat(" - Humans K =", K_h_final, "\n")
cat(" - Bots   K =", K_b_final, "\n\n")

if (save_plots) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("~/Desktop/CASEd/Reliability/Factor Analysis/Generate_Loading_Graph.R")

plot_comparison_loading_graph(
  efa_h = fa_h,
  efa_b = fa_b,
  title_text = "Chemistry: Humans vs Bots Loading Graph",
  cutoff = loading_cutoff,
  label_cutoff = 0.40,
  file = file.path(out_dir, "comparison_loading_graph.png")
)

# Optional safety check: report but do NOT step down
if (has_heywood(fa_h)) {
  warning(sprintf("[HUMANS] K=%d produced Heywood/NA issues, but no step-down was applied.", K_h_final))
}

if (has_heywood(fa_b)) {
  warning(sprintf("[BOTS] K=%d produced Heywood/NA issues, but no step-down was applied.", K_b_final))
}

# ============================================================
# Return results
# ============================================================

results <- list(
  data = list(hum = hum, bot = bot, keep_items = keep_items),
  selection = list(
    humans = list(
      kaiser = kai_h,
      K_start = K_h_start,
      K_final = K_h_final,
      rejected_from = NULL
    ),
    bots = list(
      kaiser = kai_b,
      K_start = K_b_start,
      K_final = K_b_final,
      rejected_from = NULL
    ),
    parallel = if (use_parallel_analysis) list(humans = pa_h, bots = pa_b) else NULL
  ),
  efa = list(
    humans = list(
      fa = fa_h,
      loadings = get_loadings_matrix(fa_h)
    ),
    bots = list(
      fa = fa_b,
      loadings = get_loadings_matrix(fa_b)
    )
  ),
  meta = list(
    rotate_method = rotate_method,
    fm_method = fm_method,
    loading_cutoff = loading_cutoff,
    alpha = alpha,
    pa_n_iter = pa_n_iter
  )
)



if (save_results_rds) {
  saveRDS(results, results_rds_path)
  cat("Saved results RDS to:", results_rds_path, "\n")
}



cat("\n===== HUMANS EFA (final) =====\n")
print(results$efa$humans$fa)

cat("\n===== BOTS EFA (final) =====\n")
print(results$efa$bots$fa)
