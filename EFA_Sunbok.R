# ============================================================
# EFA_Sunbok.R


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(psych)
})

# -----------------------------
# file paths
# -----------------------------
humans_csv <- "~/Desktop/CASEd/Reliability/Factor Analysis/chemistry_humans.csv"
bots_csv   <- "~/Desktop/CASEd/Reliability/Factor Analysis/chemistry_chatbots.csv"

pairs <- list()
drop_items <- character(0)

rotate_method  <- "oblimin"
fm_method      <- "minres"
loading_cutoff <- 0.30

use_parallel_analysis <- TRUE
pa_n_iter <- 30

alpha <- 0.05
boot_B <- 300

# Determinism controls
master_seed <- 20260223   
set_rng_kind <- TRUE      



# Optional: save results to RDS
save_results_rds <- TRUE
results_rds_path <- "~/Desktop/CASEd/Reliability/Factor Analysis/efa_chemistry_results.rds"

# ============================================================
# Helpers
# ============================================================

if (set_rng_kind) {
  # helps reproducibility across machines/R versions (not perfect, but better)
  suppressWarnings(RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection"))
}
set.seed(master_seed)

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

# ------------------------------------------------------------
# Bootstrap significance for loadings (deterministic)
# ------------------------------------------------------------
bootstrap_loading_sig <- function(df, nfactors, rotate, fm, B, seed, alpha) {
  stopifnot(alpha > 0 && alpha < 1)
  set.seed(seed)
  
  R0 <- tetra_cor_smooth(df)
  fit0 <- run_fa(R0, nrow(df), nfactors, rotate, fm)
  L0 <- get_loadings_matrix(fit0)
  rownames(L0) <- colnames(df)
  colnames(L0) <- paste0("F", seq_len(ncol(L0)))
  
  boot_array <- array(NA_real_, dim = c(nrow(L0), ncol(L0), B),
                      dimnames = list(rownames(L0), colnames(L0), NULL))
  
  for (b in seq_len(B)) {
    idx <- sample.int(nrow(df), replace = TRUE)
    dfb <- df[idx, , drop = FALSE]
    
    # skip replicate if any item becomes constant (tetrachoric fails)
    ok <- TRUE
    for (j in seq_len(ncol(dfb))) {
      u <- unique(dfb[[j]][!is.na(dfb[[j]])])
      if (length(u) < 2) { ok <- FALSE; break }
    }
    if (!ok) next
    
    Rb <- try(tetra_cor_smooth(dfb), silent = TRUE)
    if (inherits(Rb, "try-error")) next
    
    fb <- try(run_fa(Rb, nrow(dfb), nfactors, rotate, fm), silent = TRUE)
    if (inherits(fb, "try-error")) next
    
    Lb <- get_loadings_matrix(fb)
    rownames(Lb) <- colnames(df)
    colnames(Lb) <- paste0("F", seq_len(ncol(Lb)))
    
    # Align permutation using congruence; ensure 1-to-1 mapping via greedy unique assignment
    cong <- psych::factor.congruence(Lb, L0)  # (boot factors) x (orig factors) in most versions
    
    # Make sure dimensions match expectation
    if (nrow(cong) != ncol(Lb) || ncol(cong) != ncol(L0)) next
    
    # Greedy unique assignment: for each orig factor, pick best unused boot factor
    perm <- integer(ncol(L0))
    used <- rep(FALSE, ncol(Lb))
    for (k in seq_len(ncol(L0))) {
      scores <- cong[, k]
      scores[used] <- -Inf
      j <- which.max(scores)
      if (!is.finite(scores[j])) { perm[k] <- NA_integer_; next }
      perm[k] <- j
      used[j] <- TRUE
    }
    if (any(is.na(perm))) next
    
    Lb_aligned <- Lb[, perm, drop = FALSE]
    colnames(Lb_aligned) <- colnames(L0)
    
    # Sign alignment
    for (k in seq_len(ncol(L0))) {
      if (cor(Lb_aligned[, k], L0[, k], use = "pairwise.complete.obs") < 0) {
        Lb_aligned[, k] <- -Lb_aligned[, k]
      }
    }
    
    boot_array[, , b] <- Lb_aligned
  }
  
  lo_p <- alpha / 2
  hi_p <- 1 - alpha / 2
  lo <- apply(boot_array, c(1, 2), quantile, probs = lo_p, na.rm = TRUE)
  hi <- apply(boot_array, c(1, 2), quantile, probs = hi_p, na.rm = TRUE)
  sig <- (lo > 0) | (hi < 0)
  
  list(
    fit = fit0,
    loadings = L0,
    ci_lo = lo,
    ci_hi = hi,
    sig = sig,
    alpha = alpha,
    B = B,
    seed = seed
  )
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
# Factor count suggestion (deterministic)
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
  set.seed(master_seed + 100)  # deterministic PA (humans)
  pa_h <- psych::fa.parallel(hum, fm = fm_method, fa = "fa", cor = "tet",
                             n.iter = pa_n_iter, plot = TRUE)
  
  set.seed(master_seed + 200)  # deterministic PA (bots)
  pa_b <- psych::fa.parallel(bot, fm = fm_method, fa = "fa", cor = "tet",
                             n.iter = pa_n_iter, plot = TRUE)
  
  if (!is.null(pa_h$nfact) && !is.na(pa_h$nfact) && pa_h$nfact >= 1) K_h_start <- pa_h$nfact
  if (!is.null(pa_b$nfact) && !is.na(pa_b$nfact) && pa_b$nfact >= 1) K_b_start <- pa_b$nfact
  
  cat("Parallel Analysis suggested K (humans):", K_h_start, "\n")
  cat("Parallel Analysis suggested K (bots):  ", K_b_start, "\n")
}

# ============================================================
# Fit EFA (separate groups; step-down only if inadmissible)
# ============================================================

cat("\n=== PRIMARY: Separate EFA fits ===\n")
cat("Humans: start K =", K_h_start, "\n")
fit_h <- fit_group_solution(R_h, nrow(hum), K_h_start, rotate_method, fm_method, "HUMANS")
fa_h <- fit_h$fa
K_h_final <- fit_h$K

cat("Bots:   start K =", K_b_start, "\n")
fit_b <- fit_group_solution(R_b, nrow(bot), K_b_start, rotate_method, fm_method, "BOTS")
fa_b <- fit_b$fa
K_b_final <- fit_b$K

cat("\n[FINAL separate K]\n")
cat(" - Humans K =", K_h_final, "\n")
cat(" - Bots   K =", K_b_final, "\n\n")


# ============================================================
# Bootstrap significance (deterministic, separate seeds)
# ============================================================

cat("=== Bootstrapping loading significance ===\n")
boot_h <- bootstrap_loading_sig(
  hum, nfactors = K_h_final, rotate = rotate_method, fm = fm_method,
  B = boot_B, seed = master_seed + 1000, alpha = alpha
)

boot_b <- bootstrap_loading_sig(
  bot, nfactors = K_b_final, rotate = rotate_method, fm = fm_method,
  B = boot_B, seed = master_seed + 2000, alpha = alpha
)

boot_b_forced <- NULL
if (compute_forced_models && !is.null(fa_b_forced) && !isTRUE(forced_rejected)) {
  boot_b_forced <- bootstrap_loading_sig(
    bot, nfactors = forced_K_bots, rotate = rotate_method, fm = fm_method,
    B = boot_B, seed = master_seed + 3000, alpha = alpha
  )
}

# ============================================================
# Return results
# ============================================================

results <- list(
  data = list(hum = hum, bot = bot, keep_items = keep_items),
  selection = list(
    humans = list(kaiser = kai_h, K_start = K_h_start, K_final = K_h_final, rejected_from = fit_h$rejected_from),
    bots   = list(kaiser = kai_b, K_start = K_b_start, K_final = K_b_final, rejected_from = fit_b$rejected_from),
    parallel = if (use_parallel_analysis) list(humans = pa_h, bots = pa_b) else NULL
  ),
  efa = list(
    humans = list(fa = fa_h, loadings = boot_h$loadings, ci_lo = boot_h$ci_lo, ci_hi = boot_h$ci_hi, sig = boot_h$sig),
    bots   = list(fa = fa_b, loadings = boot_b$loadings, ci_lo = boot_b$ci_lo, ci_hi = boot_b$ci_hi, sig = boot_b$sig),
    bots_forced = if (!is.null(boot_b_forced)) list(
      K_forced = forced_K_bots,
      fa = fa_b_forced,
      loadings = boot_b_forced$loadings,
      ci_lo = boot_b_forced$ci_lo,
      ci_hi = boot_b_forced$ci_hi,
      sig = boot_b_forced$sig
    ) else if (!is.null(fa_b_forced)) list(
      K_forced = forced_K_bots,
      fa = fa_b_forced,
      inadmissible = forced_rejected
    ) else NULL
  ),
  meta = list(
    rotate_method = rotate_method,
    fm_method = fm_method,
    loading_cutoff = loading_cutoff,
    alpha = alpha,
    boot_B = boot_B,
    master_seed = master_seed,
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

if (compute_forced_models && !is.null(results$efa$bots_forced)) {
  cat(sprintf("\n===== BOTS EFA (forced K=%d) =====\n", forced_K_bots))
  print(results$efa$bots_forced$fa)
  if (!is.null(results$efa$bots_forced$inadmissible)) {
    cat("Forced model inadmissible:", results$efa$bots_forced$inadmissible, "\n")
  }
}
