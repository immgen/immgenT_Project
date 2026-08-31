# ==============================================================================
# tail_dependence_functions.R
#
# All functions for asymmetric tail dependence analysis via empirical copula.
# Source this file at the top of any analysis script:
#
#   source("tail_dependence_functions.R")
#
# Functions exported:
#   to_copula()                -- empirical copula (rank) transform
#   tail_dependence_analysis() -- main analysis; returns tail_dependence_result
#   print.tail_dependence_result()
#   plot_kappa_curve()         -- kappa(u) curve with bootstrap CI
#   plot_conditional_prob()    -- P(Y high | X high) curve with bootstrap CI
#   plot_scatter_threshold()   -- annotated scatter at a chosen u threshold
#   plot_joint_tail()          -- joint survival surface comparison (3 panels)
# ==============================================================================


# ------------------------------------------------------------------------------
# Helper: empirical copula transform
#
# Applies rank-based probability integral transform within the supplied sample.
# Called separately for each bootstrap resample so that ranks reflect that
# resample's own marginal distributions, not the global ones.
# ------------------------------------------------------------------------------

to_copula <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x    <- x[keep]
  y    <- y[keep]
  n    <- length(x)
  U    <- rank(x, ties.method = "average") / (n + 1)
  V    <- rank(y, ties.method = "average") / (n + 1)
  list(U = U, V = V, n = n)
}


# ------------------------------------------------------------------------------
# Helper: compute kappa and test statistics for one (u, v_fixed) pair
#
# Returns a single-row data frame with counts, probabilities, kappa,
# conditional probability P(Y high | X high), and Fisher p-values for both
# directions (depletion and enrichment).
# ------------------------------------------------------------------------------

compute_kappa_row <- function(U, V, u, v_fixed, min_tail_count) {
  n       <- length(U)
  idx_u   <- U > u
  idx_v   <- V > v_fixed
  n_x     <- sum(idx_u)
  n_y     <- sum(idx_v)
  n_joint <- sum(idx_u & idx_v)

  enough  <- (n_x >= min_tail_count) && (n_y >= min_tail_count)

  p_x      <- n_x     / n
  p_y      <- n_y     / n
  p_joint  <- n_joint / n
  baseline <- (1 - u) * (1 - v_fixed)

  kappa    <- if (enough && baseline > 0) p_joint / baseline else NA_real_
  delta    <- if (enough)                 p_joint - baseline else NA_real_
  cond_y_x <- if (enough && n_x > 0)     n_joint / n_x      else NA_real_

  # Two-sided contingency table:
  #               Y high          Y low
  #   X high      n_joint         n_x - n_joint
  #   X low       n_y - n_joint   n - n_x - n_y + n_joint
  #
  # Depletion  (OR < 1): fewer co-occurring cells than independence predicts.
  # Enrichment (OR > 1): more co-occurring cells than independence predicts.
  p_depletion  <- NA_real_
  p_enrichment <- NA_real_
  odds_ratio   <- NA_real_

  if (enough) {
    tab <- matrix(
      c(n_joint,         n_x - n_joint,
        n_y - n_joint,   n - n_x - n_y + n_joint),
      nrow = 2, byrow = TRUE
    )
    ft_less    <- suppressWarnings(fisher.test(tab, alternative = "less"))
    ft_greater <- suppressWarnings(fisher.test(tab, alternative = "greater"))
    p_depletion  <- ft_less$p.value
    p_enrichment <- ft_greater$p.value
    odds_ratio   <- unname(ft_less$estimate)
  }

  data.frame(
    u            = u,
    v_fixed      = v_fixed,
    n            = n,
    n_x_tail     = n_x,
    n_y_tail     = n_y,
    n_joint      = n_joint,
    p_x_tail     = p_x,
    p_y_tail     = p_y,
    p_joint      = p_joint,
    baseline     = baseline,
    delta        = delta,
    cond_y_x     = cond_y_x,
    kappa        = kappa,
    odds_ratio   = odds_ratio,
    p_depletion  = p_depletion,
    p_enrichment = p_enrichment,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# Helper: scan across u_grid for a given (U, V) pair
# ------------------------------------------------------------------------------

scan_u_grid <- function(U, V, u_grid, v_fixed, min_tail_count) {
  rows <- lapply(
    u_grid, compute_kappa_row,
    U = U, V = V, v_fixed = v_fixed, min_tail_count = min_tail_count
  )
  do.call(rbind, rows)
}


# ------------------------------------------------------------------------------
# Main function: tail_dependence_analysis
#
# Arguments:
#   x             numeric vector, the "driver" program scores (X axis, scanned)
#   y             numeric vector, the "response" program scores (Y axis, fixed)
#   v_fixed       quantile threshold for Y tail; default 0.8 (top 20%)
#   u_grid        quantile thresholds to scan for X
#   test          "depletion", "enrichment", or "both"
#   min_tail_count minimum cells required in each tail for a valid test
#   B             bootstrap replicates (0 = skip)
#   conf_level    confidence level for bootstrap CIs
#   seed          RNG seed for reproducibility
# ------------------------------------------------------------------------------

tail_dependence_analysis <- function(
    x,
    y,
    v_fixed        = 0.8,
    u_grid         = seq(0.02, 0.98, by = 0.02),
    test           = c("depletion", "enrichment", "both"),
    min_tail_count = 5,
    B              = 0,
    conf_level     = 0.95,
    seed           = NULL
) {

  test <- match.arg(test)

  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
  stopifnot(v_fixed > 0, v_fixed < 1)
  stopifnot(all(u_grid > 0 & u_grid < 1))
  u_grid <- sort(unique(u_grid))

  if (!is.null(seed)) set.seed(seed)

  cop <- to_copula(x, y)
  U   <- cop$U
  V   <- cop$V
  n   <- cop$n

  if (n < 20) stop("Need at least 20 complete observations.")

  summary_df <- scan_u_grid(U, V, u_grid, v_fixed, min_tail_count)

  if (test %in% c("depletion", "both"))
    summary_df$p_depletion_adj  <- p.adjust(summary_df$p_depletion,  method = "BH")
  if (test %in% c("enrichment", "both"))
    summary_df$p_enrichment_adj <- p.adjust(summary_df$p_enrichment, method = "BH")

  if (test %in% c("depletion", "both"))
    summary_df$sig_depletion  <- !is.na(summary_df$kappa) &
      summary_df$kappa < 1 & summary_df$p_depletion_adj < 0.05
  if (test %in% c("enrichment", "both"))
    summary_df$sig_enrichment <- !is.na(summary_df$kappa) &
      summary_df$kappa > 1 & summary_df$p_enrichment_adj < 0.05

  boot_kappa   <- NULL
  boot_cond_yx <- NULL

  if (B > 0) {
    alpha        <- 1 - conf_level
    boot_kappa   <- matrix(NA_real_, nrow = B, ncol = length(u_grid))
    boot_cond_yx <- matrix(NA_real_, nrow = B, ncol = length(u_grid))

    for (b in seq_len(B)) {
      idx      <- sample.int(n, size = n, replace = TRUE)
      boot_cop <- to_copula(x[idx], y[idx])
      tmp      <- scan_u_grid(boot_cop$U, boot_cop$V, u_grid, v_fixed, min_tail_count)
      boot_kappa[b, ]   <- tmp$kappa
      boot_cond_yx[b, ] <- tmp$cond_y_x
    }

    summary_df$kappa_lower <- apply(boot_kappa,   2, quantile,
                                    probs = alpha / 2,     na.rm = TRUE)
    summary_df$kappa_upper <- apply(boot_kappa,   2, quantile,
                                    probs = 1 - alpha / 2, na.rm = TRUE)
    summary_df$cond_yx_lower <- apply(boot_cond_yx, 2, quantile,
                                      probs = alpha / 2,     na.rm = TRUE)
    summary_df$cond_yx_upper <- apply(boot_cond_yx, 2, quantile,
                                      probs = 1 - alpha / 2, na.rm = TRUE)
  }

  structure(
    list(
      summary      = summary_df,
      U            = U,
      V            = V,
      v_fixed      = v_fixed,
      n            = n,
      test         = test,
      conf_level   = conf_level,
      boot_kappa   = boot_kappa,
      boot_cond_yx = boot_cond_yx,
      call         = match.call()
    ),
    class = "tail_dependence_result"
  )
}


# ------------------------------------------------------------------------------
# Print method
# ------------------------------------------------------------------------------

print.tail_dependence_result <- function(x, ...) {
  df <- x$summary
  cat("Tail dependence analysis\n")
  cat(sprintf("  n             : %d\n", x$n))
  cat(sprintf("  v_fixed       : %.2f  (Y top %.0f%%)\n",
              x$v_fixed, (1 - x$v_fixed) * 100))
  cat(sprintf("  u_grid        : [%.2f, %.2f], %d steps\n",
              min(df$u), max(df$u), nrow(df)))
  cat(sprintf("  test          : %s\n", x$test))
  cat(sprintf("  bootstrap B   : %s\n",
              if (is.null(x$boot_kappa)) "none" else nrow(x$boot_kappa)))
  if (x$test %in% c("depletion", "both")) {
    u_dep <- df$u[df$sig_depletion]
    cat(sprintf("  Depleted u*   : %s\n",
                if (length(u_dep) == 0) "none"
                else paste(sprintf("%.2f", u_dep), collapse = ", ")))
  }
  if (x$test %in% c("enrichment", "both")) {
    u_enr <- df$u[df$sig_enrichment]
    cat(sprintf("  Enriched u*   : %s\n",
                if (length(u_enr) == 0) "none"
                else paste(sprintf("%.2f", u_enr), collapse = ", ")))
  }
  cat("  (BH-adjusted Fisher p < 0.05 and kappa in the expected direction)\n")
  invisible(x)
}


# ------------------------------------------------------------------------------
# Plot 1: kappa(u) curve with bootstrap CI and significance markers
# ------------------------------------------------------------------------------

plot_kappa_curve <- function(res,
                             xlab = "X quantile threshold u",
                             ylab = expression(hat(kappa)(u, v)),
                             main = NULL) {

  if (!inherits(res, "tail_dependence_result"))
    stop("res must be a tail_dependence_result object.")

  df     <- res$summary
  has_ci <- "kappa_lower" %in% names(df)

  if (is.null(main))
    main <- sprintf("Tail dependence  (v_fixed = %.2f, test = %s)",
                    res$v_fixed, res$test)

  pt_col <- rep("gray50", nrow(df))
  if ("sig_depletion"  %in% names(df)) pt_col[df$sig_depletion]  <- "firebrick"
  if ("sig_enrichment" %in% names(df)) pt_col[df$sig_enrichment] <- "forestgreen"

  ylim <- range(c(df$kappa,
                  if (has_ci) c(df$kappa_lower, df$kappa_upper), 1),
                na.rm = TRUE)

  plot(df$u, df$kappa, type = "b", pch = 16, cex = 0.9,
       col = pt_col, xlab = xlab, ylab = ylab, ylim = ylim, main = main)
  abline(h = 1, lty = 2, col = "gray40")

  if (has_ci)
    segments(df$u, df$kappa_lower, df$u, df$kappa_upper,
             col = "gray70", lwd = 0.8)

  if ("sig_depletion" %in% names(df) && any(df$sig_depletion, na.rm = TRUE))
    points(df$u[df$sig_depletion], df$kappa[df$sig_depletion],
           pch = 8, col = "firebrick", cex = 1.3)
  if ("sig_enrichment" %in% names(df) && any(df$sig_enrichment, na.rm = TRUE))
    points(df$u[df$sig_enrichment], df$kappa[df$sig_enrichment],
           pch = 8, col = "forestgreen", cex = 1.3)

  leg_lab <- "independence (kappa = 1)"
  leg_col <- "gray40"; leg_lty <- 2; leg_pch <- NA
  if (res$test %in% c("depletion", "both")) {
    leg_lab <- c(leg_lab, "depleted (BH p < 0.05)")
    leg_col <- c(leg_col, "firebrick"); leg_lty <- c(leg_lty, NA); leg_pch <- c(leg_pch, 8)
  }
  if (res$test %in% c("enrichment", "both")) {
    leg_lab <- c(leg_lab, "enriched (BH p < 0.05)")
    leg_col <- c(leg_col, "forestgreen"); leg_lty <- c(leg_lty, NA); leg_pch <- c(leg_pch, 8)
  }
  legend("topright", legend = leg_lab, col = leg_col,
         lty = leg_lty, pch = leg_pch, cex = 0.8, bty = "n")
}


# ------------------------------------------------------------------------------
# Plot 2: P(Y high | X high) with bootstrap CI ribbon
# ------------------------------------------------------------------------------

plot_conditional_prob <- function(res,
                                  xlab = "X quantile threshold u",
                                  ylab = expression(hat(P)(Y ~ high ~ "|" ~ X ~ high)),
                                  main = "Conditional probability P(Y high | X high)") {

  if (!inherits(res, "tail_dependence_result"))
    stop("res must be a tail_dependence_result object.")

  df      <- res$summary
  p_y_ref <- mean(df$p_y_tail, na.rm = TRUE)
  has_ci  <- "cond_yx_lower" %in% names(df)

  ylim <- range(c(df$cond_y_x, p_y_ref,
                  if (has_ci) c(df$cond_yx_lower, df$cond_yx_upper)),
                na.rm = TRUE)
  ylim <- ylim + c(-1, 1) * 0.02 * diff(ylim)

  plot(df$u, df$cond_y_x, type = "n",
       xlab = xlab, ylab = ylab, ylim = ylim, main = main)

  if (has_ci) {
    ok <- !is.na(df$cond_yx_lower) & !is.na(df$cond_yx_upper)
    if (any(ok))
      polygon(c(df$u[ok], rev(df$u[ok])),
              c(df$cond_yx_lower[ok], rev(df$cond_yx_upper[ok])),
              col = adjustcolor("steelblue", alpha.f = 0.20), border = NA)
  }

  lines(df$u, df$cond_y_x, type = "b", pch = 16, cex = 0.8, col = "steelblue")
  abline(h = p_y_ref, lty = 2, col = "firebrick")

  leg_lab  <- c("P(Y high | X > u)",
                sprintf("Marginal P(Y high) = %.3f", p_y_ref))
  leg_col  <- c("steelblue", "firebrick")
  leg_lty  <- c(1, 2); leg_pch <- c(16, NA); leg_fill <- c(NA, NA)

  if (has_ci) {
    leg_lab  <- c(leg_lab, sprintf("%d%% bootstrap CI", round(res$conf_level * 100)))
    leg_col  <- c(leg_col, "steelblue"); leg_lty <- c(leg_lty, NA)
    leg_pch  <- c(leg_pch, NA)
    leg_fill <- c(leg_fill, adjustcolor("steelblue", alpha.f = 0.20))
  }

  legend("topright", legend = leg_lab, col = leg_col,
         lty = leg_lty, pch = leg_pch, fill = leg_fill,
         border = NA, cex = 0.8, bty = "n")
}


# ------------------------------------------------------------------------------
# Plot 3: annotated scatter at a chosen u threshold
# ------------------------------------------------------------------------------

plot_scatter_threshold <- function(res, x_raw, y_raw,
                                   u_thresh = NULL,
                                   x_label  = "X Score",
                                   y_label  = "Y Score",
                                   main     = NULL,
                                   show_annotation = TRUE) {

  if (!inherits(res, "tail_dependence_result"))
    stop("res must be a tail_dependence_result object.")

  df <- res$summary

  if (is.null(u_thresh)) {
    sig_u    <- c(if ("sig_depletion"  %in% names(df)) df$u[df$sig_depletion],
                  if ("sig_enrichment" %in% names(df)) df$u[df$sig_enrichment])
    u_thresh <- if (length(sig_u) > 0) max(sig_u) else median(df$u)
  }

  row <- df[abs(df$u - u_thresh) < 1e-9, ]
  if (nrow(row) == 0) stop("u_thresh not found in u_grid.")

  x_cut <- quantile(x_raw, probs = u_thresh,   na.rm = TRUE)
  y_cut <- quantile(y_raw, probs = res$v_fixed, na.rm = TRUE)

  if (is.null(main)) {
    main <- sprintf("%s vs %s  |  u = %.2f, v = %.2f",
                    y_label, x_label, u_thresh, res$v_fixed)
  }

  plot(y_raw ~ x_raw, pch = 16, cex = 0.4, col = rgb(0, 0, 1, 0.1),
       xlab = x_label, ylab = y_label,
       main = main)
  abline(v = x_cut, col = "firebrick", lty = 2)
  abline(h = y_cut, col = "firebrick", lty = 2)

  if (!isTRUE(show_annotation)) return(invisible(NULL))

  usr    <- par("usr")
  labels <- c(
    sprintf("u = %.2f,  v = %.2f",      u_thresh,        res$v_fixed),
    sprintf("P(X > u)         = %.3f",  row$p_x_tail),
    sprintf("P(Y > v)         = %.3f",  row$p_y_tail),
    sprintf("P(joint)         = %.3f",  row$p_joint),
    sprintf("Expected (indep) = %.3f",  row$baseline),
    sprintf("kappa            = %.3f",  row$kappa),
    sprintf("P(Y>v | X>u)     = %.3f",  row$cond_y_x)
  )
  if ("p_depletion_adj"  %in% names(row))
    labels <- c(labels, sprintf("Fisher p dep  (BH) = %.4f", row$p_depletion_adj))
  if ("p_enrichment_adj" %in% names(row))
    labels <- c(labels, sprintf("Fisher p enr  (BH) = %.4f", row$p_enrichment_adj))

  text(x      = usr[1] + 0.03 * diff(usr[1:2]),
       y      = usr[4] - seq(0.05, by = 0.065, length.out = length(labels)) *
                  diff(usr[3:4]),
       labels = labels, adj = c(0, 1), cex = 0.80)
}

plot_kappa_heatmap_square <- function(
    x,
    y,
    x_label,
    y_label,
    n_grid     = 80,
    v_fixed    = 0.80,
    zlim_kappa = NULL,
    main       = ""
) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]; n <- length(x)
  stopifnot(n >= 20)

  # Use strictly increasing grids; duplicated quantiles can create blank bands in image().
  probs  <- seq(0.01, 0.99, length.out = n_grid)
  x_grid <- sort(unique(as.numeric(quantile(x, probs = probs, na.rm = TRUE))))
  y_grid <- sort(unique(as.numeric(quantile(y, probs = probs, na.rm = TRUE))))
  if (length(x_grid) < 2 || length(y_grid) < 2)
    stop("Not enough unique grid points to draw heatmap.")

  # Use >= to avoid zero-probability edges when a grid point equals the sample max.
  px <- sapply(x_grid, function(xi) mean(x >= xi))
  py <- sapply(y_grid, function(yi) mean(y >= yi))

  S_emp <- outer(x_grid, y_grid, FUN = Vectorize(function(a, b) mean(x >= a & y >= b)))
  S_ind <- outer(px, py)

  kappa_mat             <- S_emp / S_ind
  kappa_mat[S_ind == 0] <- NA_real_

  if (is.null(zlim_kappa)) {
    rng        <- range(kappa_mat, na.rm = TRUE)
    max_dev    <- max(abs(rng - 1))
    zlim_kappa <- c(1 - max_dev, 1 + max_dev)
  }

  kappa_lo <- zlim_kappa[1]; kappa_hi <- zlim_kappa[2]
  mid_frac <- max(0.01, min(0.99, (1 - kappa_lo) / (kappa_hi - kappa_lo)))
  n_lo     <- round(mid_frac * 256); n_hi <- 256 - n_lo
  pal_kappa <- c(
    colorRampPalette(c("#2166ac", "#4393c3", "#d1e5f0", "#f7f7f7"))(n_lo),
    colorRampPalette(c("#f7f7f7", "#fddbc7", "#d6604d", "#b2182b"))(n_hi)
  )

  # Do not use asp=1 with pty="s" here: when marginal score ranges differ a lot in
  # width (common for TH1 vs TFH/Tchron), R letterboxes the heatmap and large
  # white margins appear inside the panel. Stretch x/y independently so the
  # heatmap fills the square PDF page.
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(5, 5, 2, 2))

  image(
    x_grid, y_grid, kappa_mat,
    col  = pal_kappa,
    zlim = zlim_kappa,
    xlab = x_label,
    ylab = y_label,
    main = main
  )
  box()

  # Reference lines at v_fixed marginal quantile on both axes (raw score space)
  x_ref <- unname(quantile(x, probs = v_fixed, na.rm = TRUE))
  y_ref <- unname(quantile(y, probs = v_fixed, na.rm = TRUE))
  abline(v = x_ref, lty = 2, lwd = 1.4, col = "gray20")
  abline(h = y_ref, lty = 2, lwd = 1.4, col = "gray20")

  invisible(list(
    x_grid = x_grid, y_grid = y_grid,
    S_emp  = S_emp,  S_ind  = S_ind,
    kappa  = kappa_mat
  ))
}


# ------------------------------------------------------------------------------
# Plot 4: joint tail probability comparison (3 panels)
#
#   Left   -- -log10 P(X > x, Y > y)  [empirical]
#   Centre -- -log10 P(X > x) P(Y > y)  [independence]
#   Right  -- kappa(x, y) = empirical / independence  (linear scale, no log)
#
# Panels 1 and 2 share the same -log10 colour scale, so a grid cell that is
# darker in panel 1 than panel 2 signals depletion of the empirical joint tail
# relative to independence. White cells in panel 3 are where the denominator
# is exactly zero (no marginal mass) -- no NaN artefacts.
# ------------------------------------------------------------------------------

plot_joint_tail <- function(
    x,
    y,
    n_grid     = 80,
    x_label    = "X Score",
    y_label    = "Y Score",
    zlim_kappa = NULL
) {

  keep <- is.finite(x) & is.finite(y)
  x    <- x[keep]; y <- y[keep]; n <- length(x)
  stopifnot(n >= 20)

  x_grid <- quantile(x, probs = seq(0.01, 0.99, length.out = n_grid))
  y_grid <- quantile(y, probs = seq(0.01, 0.99, length.out = n_grid))

  px <- sapply(x_grid, function(xi) mean(x > xi))
  py <- sapply(y_grid, function(yi) mean(y > yi))

  S_emp <- outer(x_grid, y_grid, FUN = function(xi, yi)
    mapply(function(a, b) mean(x > a & y > b), xi, yi))
  S_ind <- outer(px, py)

  kappa_mat             <- S_emp / S_ind
  kappa_mat[S_ind == 0] <- NA_real_

  floor_emp     <- 0.5 / n
  floor_ind     <- (0.5 / n)^2
  neg_log10_emp <- -log10(pmax(S_emp, floor_emp))
  neg_log10_ind <- -log10(pmax(S_ind, floor_ind))
  zlim_log      <- range(c(neg_log10_emp, neg_log10_ind), na.rm = TRUE)

  if (is.null(zlim_kappa)) {
    rng        <- range(kappa_mat, na.rm = TRUE)
    max_dev    <- max(abs(rng - 1))
    zlim_kappa <- c(1 - max_dev, 1 + max_dev)
  }

  pal_surv <- colorRampPalette(
    c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(256)

  kappa_lo <- zlim_kappa[1]; kappa_hi <- zlim_kappa[2]
  mid_frac <- max(0.01, min(0.99, (1 - kappa_lo) / (kappa_hi - kappa_lo)))
  n_lo     <- round(mid_frac * 256); n_hi <- 256 - n_lo
  pal_kappa <- c(
    colorRampPalette(c("#2166ac", "#4393c3", "#d1e5f0", "#f7f7f7"))(n_lo),
    colorRampPalette(c("#f7f7f7", "#fddbc7", "#d6604d", "#b2182b"))(n_hi)
  )

  add_prob_contours <- function(x_g, y_g, neg_log_mat, prob_levels) {
    log_levels <- -log10(prob_levels)
    valid      <- log_levels >= zlim_log[1] & log_levels <= zlim_log[2]
    if (!any(valid)) return(invisible(NULL))
    contour(x_g, y_g, neg_log_mat,
            levels = log_levels[valid], labels = as.character(prob_levels[valid]),
            col = "white", lwd = 0.6, labcex = 0.55, add = TRUE)
  }

  prob_levels <- c(0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.001)

  old_par <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 2))
  on.exit(par(old_par))

  image(x_grid, y_grid, neg_log10_emp,
        col = pal_surv, zlim = zlim_log, xlab = x_label, ylab = y_label,
        main = expression("-" * log[10] ~ P(X > x * "," ~ Y > y) ~ "[empirical]"))
  add_prob_contours(x_grid, y_grid, neg_log10_emp, prob_levels)
  box()

  image(x_grid, y_grid, neg_log10_ind,
        col = pal_surv, zlim = zlim_log, xlab = x_label, ylab = y_label,
        main = expression("-" * log[10] ~ P(X>x) * P(Y>y) ~ "[independence]"))
  add_prob_contours(x_grid, y_grid, neg_log10_ind, prob_levels)
  box()

  image(x_grid, y_grid, kappa_mat,
        col = pal_kappa, zlim = zlim_kappa, xlab = x_label, ylab = y_label,
        main = expression(kappa(x * "," ~ y) * " = " *
                            frac(P(X>x*","~Y>y), P(X>x)*P(Y>y))))
  # No kappa = 1 contour per request
  box()

  invisible(list(x_grid = x_grid, y_grid = y_grid,
                 S_emp = S_emp, S_ind = S_ind, kappa = kappa_mat,
                 neg_log10_emp = neg_log10_emp, neg_log10_ind = neg_log10_ind))
}
