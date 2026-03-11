# ============================================================
# Generate Factor Loadings Graph
# Manual drawing version for exact edge anchoring
# ============================================================

plot_comparison_loading_graph <- function(efa_h,
                                          efa_b,
                                          title_text   = "Comparison Loading Graph",
                                          cutoff       = 0.30,
                                          label_cutoff = 0.50,
                                          file         = NULL) {
  # -----------------------------
  # Extract loading matrices
  # -----------------------------
  Lh <- as.matrix(unclass(efa_h$loadings))
  Lb <- as.matrix(unclass(efa_b$loadings))
  
  common_items <- intersect(rownames(Lh), rownames(Lb))
  if (length(common_items) == 0) {
    stop("No common items found between efa_h and efa_b loadings.")
  }
  
  Lh <- Lh[common_items, , drop = FALSE]
  Lb <- Lb[common_items, , drop = FALSE]
  
  human_factors <- paste0("H_", colnames(Lh))
  bot_factors   <- paste0("B_", colnames(Lb))
  items         <- common_items
  
  colnames(Lh) <- human_factors
  colnames(Lb) <- bot_factors
  
  # -----------------------------
  # Build edge data
  # -----------------------------
  eh <- which(abs(Lh) >= cutoff, arr.ind = TRUE)
  edges_h <- data.frame(
    from       = colnames(Lh)[eh[, 2]],
    to         = rownames(Lh)[eh[, 1]],
    weight     = Lh[eh],
    abs_weight = abs(Lh[eh]),
    group      = "human",
    stringsAsFactors = FALSE
  )
  
  eb <- which(abs(Lb) >= cutoff, arr.ind = TRUE)
  edges_b <- data.frame(
    from       = rownames(Lb)[eb[, 1]],
    to         = colnames(Lb)[eb[, 2]],
    weight     = Lb[eb],
    abs_weight = abs(Lb[eb]),
    group      = "bot",
    stringsAsFactors = FALSE
  )
  
  edge_df <- rbind(edges_h, edges_b)
  
  # -----------------------------
  # Layout
  # -----------------------------
  nH <- length(human_factors)
  nI <- length(items)
  nB <- length(bot_factors)
  
  # Wider item spacing
  x_h <- seq(4, 22, length.out = nH)
  x_i <- seq(1, 25, length.out = nI)
  x_b <- seq(4, 22, length.out = nB)
  
  y_h <- 10.0
  y_i <- 5.4
  y_b <- 1.0
  
  coords <- rbind(
    cbind(x_h, rep(y_h, nH)),
    cbind(x_i, rep(y_i, nI)),
    cbind(x_b, rep(y_b, nB))
  )
  rownames(coords) <- c(human_factors, items, bot_factors)
  
  # -----------------------------
  # Geometry settings
  # -----------------------------
  factor_r <- 1.25
  box_w    <- 0.72
  box_h    <- 0.46
  
  # -----------------------------
  # Output device
  # -----------------------------
  if (!is.null(file)) {
    png(file, width = 2800, height = 1600, res = 220)
  }
  
  op <- par(mar = c(1, 1, 3, 1), xpd = NA)
  on.exit({
    par(op)
    if (!is.null(file)) dev.off()
  }, add = TRUE)
  
  plot(
    NA,
    xlim = c(0, 26),
    ylim = c(0, 12),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = title_text
  )
  
  # -----------------------------
  # Helper: draw factor circles
  # -----------------------------
  draw_factor_node <- function(x, y, label, fill) {
    symbols(
      x, y,
      circles = factor_r,
      inches = FALSE,
      bg = fill,
      fg = "gray35",
      add = TRUE
    )
    text(
      x, y,
      labels = label,
      cex = 1.0,
      col = "navy",
      family = "sans"
    )
  }
  
  # -----------------------------
  # Helper: draw item box
  # -----------------------------
  draw_item_box <- function(x, y, label) {
    rect(
      x - box_w / 2, y - box_h / 2,
      x + box_w / 2, y + box_h / 2,
      col = "white",
      border = "gray35",
      lwd = 1.0
    )
    text(
      x, y,
      labels = label,
      cex = 0.88,
      col = "navy",
      family = "sans"
    )
  }
  
  # -----------------------------
  # Draw edges first
  # so nodes/boxes sit on top cleanly
  # -----------------------------
  
  # Human edges: factor bottom edge -> item top middle
  if (nrow(edges_h) > 0) {
    for (i in seq_len(nrow(edges_h))) {
      from <- edges_h$from[i]
      to   <- edges_h$to[i]
      w    <- edges_h$weight[i]
      
      x1 <- coords[from, 1]
      y1 <- coords[from, 2] - factor_r
      
      x2 <- coords[to, 1]
      y2 <- coords[to, 2] + box_h / 2
      
      segments(
        x1, y1, x2, y2,
        col = "steelblue4",
        lwd = 1.5,
        lty = ifelse(w < 0, 2, 1)
      )
    }
  }
  
  # Bot edges: item bottom middle -> factor top edge
  if (nrow(edges_b) > 0) {
    for (i in seq_len(nrow(edges_b))) {
      from <- edges_b$from[i]
      to   <- edges_b$to[i]
      w    <- edges_b$weight[i]
      
      x1 <- coords[from, 1]
      y1 <- coords[from, 2] - box_h / 2
      
      x2 <- coords[to, 1]
      y2 <- coords[to, 2] + factor_r
      
      segments(
        x1, y1, x2, y2,
        col = "firebrick3",
        lwd = 1.5,
        lty = ifelse(w < 0, 2, 1)
      )
    }
  }
  
  # -----------------------------
  # Draw factor nodes
  # -----------------------------
  for (hf in human_factors) {
    draw_factor_node(coords[hf, 1], coords[hf, 2], hf, "lightblue")
  }
  
  for (bf in bot_factors) {
    draw_factor_node(coords[bf, 1], coords[bf, 2], bf, "mistyrose")
  }
  
  # -----------------------------
  # Draw item boxes
  # -----------------------------
  for (it in items) {
    draw_item_box(coords[it, 1], coords[it, 2], it)
  }
  
  # -----------------------------
  # Helper: loading label with white bg
  # -----------------------------
  draw_edge_label <- function(x, y, label, cex = 0.74) {
    w <- strwidth(label, cex = cex)
    h <- strheight(label, cex = cex)
    
    rect(
      x - w * 0.62, y - h * 0.68,
      x + w * 0.62, y + h * 0.68,
      col = "white",
      border = NA
    )
    
    text(
      x, y,
      labels = label,
      cex = cex,
      col = "navy",
      family = "sans"
    )
  }
  
  # -----------------------------
  # Rank labels within each item/group
  # to spread them apart
  # -----------------------------
  if (nrow(edge_df) > 0) {
    edge_df$label_rank <- ave(
      edge_df$abs_weight,
      paste(edge_df$group,
            ifelse(edge_df$group == "human", edge_df$to, edge_df$from)),
      FUN = function(z) rank(-z, ties.method = "first")
    )
  }
  
  # -----------------------------
  # Draw labels for stronger edges
  # -----------------------------
  if (nrow(edge_df) > 0) {
    for (i in seq_len(nrow(edge_df))) {
      w <- edge_df$weight[i]
      if (abs(w) < label_cutoff) next
      
      from <- edge_df$from[i]
      to   <- edge_df$to[i]
      
      if (edge_df$group[i] == "human") {
        x1 <- coords[from, 1]
        y1 <- coords[from, 2] - factor_r
        x2 <- coords[to, 1]
        y2 <- coords[to, 2] + box_h / 2
        t  <- 0.72
        vertical_nudge <- 0.10
      } else {
        x1 <- coords[from, 1]
        y1 <- coords[from, 2] - box_h / 2
        x2 <- coords[to, 1]
        y2 <- coords[to, 2] + factor_r
        t  <- 0.28
        vertical_nudge <- -0.10
      }
      
      dx <- x2 - x1
      dy <- y2 - y1
      len <- sqrt(dx^2 + dy^2)
      if (len == 0) next
      
      px <- -dy / len
      py <-  dx / len
      
      lx <- x1 + t * dx
      ly <- y1 + t * dy
      
      k <- edge_df$label_rank[i]
      offset_seq <- c(0.00, 0.18, -0.18, 0.32, -0.32, 0.46, -0.46)
      off <- offset_seq[((k - 1) %% length(offset_seq)) + 1]
      
      lx <- lx + off * px
      ly <- ly + off * py + vertical_nudge
      
      draw_edge_label(lx, ly, sprintf("%.2f", w), cex = 0.74)
    }
  }
  
  # -----------------------------
  # Legend
  # -----------------------------
  legend(
    "topleft",
    legend = c("Human loadings", "Bot loadings", "Negative loading = dashed"),
    col    = c("steelblue4", "firebrick3", "black"),
    lty    = c(1, 1, 2),
    lwd    = c(2.5, 2.5, 2),
    bty    = "n",
    cex    = 1.1,
    text.col = "black"
  )
}