####################################################################################
# CONDITIONAL PROBABILITY OF RECESSION AND OF FINANCIAL STRESS - CHAPTER 1 (SECTIONS 1)
####################################################################################

# Camille Souffron - MASTER THESIS APE (PSE & ENS)

#NB: all packages in the DATA_LOADING file



###############################################################################
## Original R|R 
###############################################################################

# 0) Helpers
as_qtr_index <- function(date_vec) {
  # assume Date or yearqtr-like; just ensure strictly increasing integer index
  if (inherits(date_vec, "Date")) return(order(order(date_vec)))
  order(seq_along(date_vec)) 
}


# First peak" selector: max p-hat in [7,56], ties -> smallest k
# Works with any data frame containing columns {k, p_hat} and group vars.
first_peak_k <- function(df,
                         group_vars,      
                         k_col  = "k",
                         p_col  = "p_hat",
                         k_min  = 7L,
                         k_max  = 56L,
                         digits = 6) {
  stopifnot(all(c(group_vars, k_col, p_col) %in% names(df)))
  d <- df %>%
    dplyr::filter(is.finite(.data[[k_col]]), is.finite(.data[[p_col]])) %>%
    dplyr::filter(.data[[k_col]] >= k_min, .data[[k_col]] <= k_max) %>%
    dplyr::mutate(.p_adj = round(.data[[p_col]], digits))
  if (!nrow(d)) return(tibble::tibble()[0,])
  
  d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::mutate(.p_max = max(.p_adj)) %>%
    dplyr::filter(.p_adj == .p_max) %>%
    dplyr::slice_min(order_by = .data[[k_col]], n = 1, with_ties = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      dplyr::across(dplyr::all_of(group_vars)),
      first_peak_k    = .data[[k_col]],
      first_peak_prob = .data[[p_col]]
    )
}

# Tiny consistency check (use in snapshots)
assert_peaks_range <- function(peaks_df, k_min = 7L, k_max = 56L){
  bad <- with(peaks_df, is.finite(first_peak_k) &
                (first_peak_k < k_min | first_peak_k > k_max))
  if (any(bad)) stop("Peak k out of range [", k_min, ",", k_max, "].")
  invisible(TRUE)
}

make_window_indicator <- function(r, k, x){
  # y_t^{k,x} = 1 if ANY recession occurs in [t+k-x, t+k+x]
  Tn <- length(r)
  y  <- integer(Tn)
  for (t in seq_len(Tn)) {
    lo <- t + k - x
    hi <- t + k + x
    if (lo < 1L || hi > Tn) { y[t] <- NA_integer_; next }
    y[t] <- as.integer(any(r[lo:hi] == 1L))
  }
  y
}

nw_var_alpha_plus_beta <- function(fit, lag = NULL){
  # Var(β0 + β1) with HAC (Newey–West)
  S <- if (is.null(lag)) sandwich::NeweyWest(fit) else sandwich::NeweyWest(fit, lag = lag)
  e <- c(1,1)
  as.numeric(t(e) %*% S %*% e)
}

default_nw_lag <- function(Tn){
  # generic plug-in lag length for quarterly data
  floor(4*(Tn/100)^(2/9))
}

# 1) Main estimator for a given x and k-grid 
recession_probability_curve <- function(df, k_seq = 12:90, x_window = 5,
                                        min_gap_from_current = 12,
                                        restrict_to = NULL,
                                        nw_lag = NULL) {
  stopifnot(all(c("date","recession") %in% names(df)))
  df <- df %>% arrange(date)
  if (!is.null(restrict_to)) {
    df <- df %>% filter(date >= restrict_to[1], date <= restrict_to[2])
  }
  r <- as.integer(df$recession)
  Tn <- length(r)
  if (is.null(nw_lag)) nw_lag <- default_nw_lag(Tn)
  
  out <- lapply(k_seq, function(k){
    if (k < min_gap_from_current) return(tibble(k=k, p_hat=NA_real_, se=NA_real_, n=NA_integer_))
    y <- make_window_indicator(r, k, x_window)
    keep <- is.finite(y)
    if (sum(keep) < 30) return(tibble(k=k, p_hat=NA_real_, se=NA_real_, n=sum(keep)))
    fit <- lm(y[keep] ~ r[keep])
    p_hat <- sum(coef(fit)[1:2])
    var_p <- nw_var_alpha_plus_beta(fit, lag = nw_lag)
    se    <- sqrt(max(var_p, 0))
    tibble(k=k, p_hat=p_hat, se=se, n=sum(keep))
  })
  bind_rows(out) %>% mutate(x = x_window)
}

# 2) Pretty plotting wrapper for multiple x 
plot_recession_prob <- function(df, k_seq = 12:90, x_values = c(3,4,5),
                                restrict_to = NULL, level_bands = c(.66,.80,.90),
                                nw_lag = NULL, title = "Conditional probability of recession") {
  curves <- lapply(x_values, function(x)
    recession_probability_curve(df, k_seq, x_window = x,
                                restrict_to = restrict_to, nw_lag = nw_lag)
  )
  dfc <- dplyr::bind_rows(curves)
  
  bands <- lapply(level_bands, function(lev){
    z <- qnorm((1+lev)/2)
    dfc %>%
      dplyr::transmute(k, x, lev = paste0(round(100*lev),"%"),
                       lo = pmax(0, p_hat - z*se),
                       hi = pmin(1, p_hat + z*se))
  }) %>% dplyr::bind_rows()
  
  # choose the largest x for shading
  x_max <- max(dfc$x, na.rm = TRUE)
  bands_max <- bands[bands$x == x_max, , drop = FALSE]
  
  # Pick a distinct qualitative palette
  col_map <- c("3" = "#1b9e77",  # green
               "4" = "#d95f02",  # orange
               "5" = "#7570b3")  # purple
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = bands_max,
      ggplot2::aes(k, ymin = lo, ymax = hi, fill = lev),
      alpha = 0.25
    ) +
    ggplot2::scale_fill_brewer(palette = "Greys", direction = -1,
                               name = paste0("Bands (x = ", x_max, ")")) +
    ggplot2::geom_line(
      data = dfc,
      ggplot2::aes(k, p_hat, color = factor(x)),
      linewidth = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = col_map,
      breaks = names(col_map),
      labels = paste0("x = ", names(col_map)),
      name = "Window width ±x"
    ) +
    ggplot2::coord_cartesian(xlim = range(k_seq), ylim = c(0,1)) +
    ggplot2::labs(x = "k (quarters since recession)", y = "Probability",
                  title = title,
                  subtitle = "Solid: point estimates; shaded: HAC CIs (largest x)") +
    ggplot2::theme_minimal()
  
  
  list(data = dfc, bands = bands, plot = p)
}



# 3) Joint X² tests 
# 3a) Peak=>peak bootstrap covariance 
find_peaks_from_recessions <- function(r){
  # peaks: start of recessions (0 -> 1). If series starts in recession, include t=1.
  up <- c(0L, diff(r))
  idx <- which(up == 1L)
  if (r[1] == 1L) idx <- c(1L, idx)
  if (!length(idx)) stop("No peaks found (no recession starts).")
  idx
}

resample_peak_to_peak <- function(r, T_target){
  peaks <- find_peaks_from_recessions(r)
  # Build list of peak→peak segments (last segment ends at T)
  segs <- lapply(seq_along(peaks), function(i){
    a <- peaks[i]
    b <- if (i < length(peaks)) peaks[i+1]-1L else length(r)
    r[a:b]
  })
  # sample segments with replacement until reaching T_target
  out <- integer(0)
  while (length(out) < T_target) {
    s <- segs[[sample.int(length(segs), 1)]]
    out <- c(out, s)
  }
  out[seq_len(T_target)]
}

compute_p_vec <- function(r, K, x){
  # returns p_hat across selected K using OLS+HAC se only for consistency;
  # here we just need point estimates (means conditional on r_t=1)
  sapply(K, function(k){
    y <- make_window_indicator(r, k, x)
    keep <- is.finite(y) & (r == 1L)
    if (!any(keep)) return(NA_real_)
    mean(y[keep])
  })
}

joint_test_bootstrap <- function(df, K, x = 5, B = 10000, nw_lag = NULL){
  r <- as.integer(df %>% arrange(date) %>% pull(recession))
  Tn <- length(r)
  if (is.null(nw_lag)) nw_lag <- default_nw_lag(Tn)
  
  p_obs <- compute_p_vec(r, K, x)
  if (any(!is.finite(p_obs))) stop("Not enough data for chosen K window(s).")
  
  # bootstrap covariance Σ̂_boot
  Pmat <- matrix(NA_real_, nrow = B, ncol = length(K))
  set.seed(123)
  for (b in 1:B) {
    rb <- resample_peak_to_peak(r, Tn)
    Pmat[b, ] <- compute_p_vec(rb, K, x)
  }
  Sigma_hat <- stats::cov(Pmat, use = "pairwise.complete.obs")
  
  # X² stat
  pbar <- mean(p_obs)
  evec <- rep(1, length(K))
  # robust inverse sqrt via eigen
  ev <- eigen(Sigma_hat, symmetric=TRUE)
  # regularize tiny eigenvalues
  lam <- pmax(ev$values, .Machine$double.eps)
  Sminushalf <- ev$vectors %*% diag(1/sqrt(lam)) %*% t(ev$vectors)
  eta <- as.numeric(Sminushalf %*% (p_obs - pbar*evec))
  stat <- sum(eta^2)
  pval <- 1 - pchisq(stat, df = length(K) - 1L)
  list(stat = stat, df = length(K) - 1L, pval = pval,
       Sigma_hat = Sigma_hat, p_obs = p_obs, K = K, x = x)
}

# 3b) HAC-asymptotic covariance across K (quick-and-dirty)
# Treat each p_hat_k as a time-avg of z_t(k) = r_t * y_t^{k,x} / E[r_t]; estimate
# cross-covariances with multivariate Newey–West on stacked demeaned z_t(k).
joint_test_hac <- function(df, K, x = 5, L = NULL){
  df <- df %>% arrange(date)
  r  <- as.integer(df$recession)
  Tn <- length(r)
  
  # build matrix Z_t (T x m) with columns z_t(k) = r_t * y_t^{k,x}
  Z <- sapply(K, function(k) make_window_indicator(r, k, x)) # T x m (with NA near ends)
  keep <- apply(is.finite(Z), 1, all)
  Z <- Z[keep, , drop=FALSE]; r_keep <- r[keep]
  if (!all(r_keep %in% c(0L,1L))) stop("recession must be 0/1.")
  denom <- mean(r_keep)
  Z <- (r_keep * Z) / denom  # so that E[Z(·)] = p_{k,x}
  p_obs <- colMeans(Z)
  
  # multivariate NW: Σ̂ = Γ_0 + sum_{l=1..L} w_l (Γ_l + Γ_lᵀ), Γ_l = T^{-1} sum ( (Z_t-μ)(Z_{t-l}-μ)ᵀ )
  Zc <- scale(Z, center = TRUE, scale = FALSE); T0 <- nrow(Zc)
  T0 <- nrow(Zc)
  if (is.null(L)) L <- default_nw_lag(T0)  # <<< T0, not raw Tn
  
  Gamma0 <- crossprod(Zc)/T0
  w <- 1 - (1:L)/(L+1)
  Gsum <- matrix(0, ncol(Zc), ncol(Zc))
  for (l in 1:L) {
    A <- t(Zc[(l+1):T0, , drop=FALSE]) %*% Zc[1:(T0-l), , drop=FALSE] / T0
    Gsum <- Gsum + w[l] * (A + t(A))
  }
Sigma_hat <- (Gamma0 + Gsum) / T0
  
  pbar <- mean(p_obs); evec <- rep(1, length(K))
  ev <- eigen(Sigma_hat, symmetric=TRUE)
  lam <- pmax(ev$values, .Machine$double.eps)
  Sminushalf <- ev$vectors %*% diag(1/sqrt(lam)) %*% t(ev$vectors)
  eta <- as.numeric(Sminushalf %*% (p_obs - pbar*evec))
  stat <- sum(eta^2)
  pval <- 1 - pchisq(stat, df = length(K) - 1L)
  list(stat = stat, df = length(K) - 1L, pval = pval,
       Sigma_hat = Sigma_hat, p_obs = p_obs, K = K, x = x)
}



# 1) NBER recession intervals (monthly, official)
nber_month_ranges <- tribble(
  ~start,      ~end,
  "1948-11-01","1949-10-01",
  "1953-07-01","1954-05-01",
  "1957-08-01","1958-04-01",
  "1960-04-01","1961-02-01",
  "1969-12-01","1970-11-01",
  "1973-11-01","1975-03-01",
  "1980-01-01","1980-07-01",
  "1981-07-01","1982-11-01",
  "1990-07-01","1991-03-01",
  "2001-03-01","2001-11-01",
  "2007-12-01","2009-06-01",
) %>%
  mutate(start = as.Date(start), end = as.Date(end))


start_q <- lubridate::floor_date(min(nber_month_ranges$start), "quarter")
end_q   <- lubridate::floor_date(max(nber_month_ranges$end),   "quarter")
q_dates <- seq(from = start_q, to = end_q, by = "quarter")


# 3) Flag quarters with any recession month
#    (a quarter is recession if ANY day in that quarter is inside an NBER interval)
is_recession_quarter <- function(q_date){
  # months inside this quarter
  q_start <- q_date
  q_end   <- q_date %m+% months(3) - days(1)
  any(apply(nber_month_ranges, 1, function(r){
    r_start <- r[[1]]; r_end <- r[[2]]
    # overlap test between [q_start, q_end] and [r_start, r_end]
    (q_start <= r_end) && (r_start <= q_end)
  }))
}

recession_flag <- vapply(q_dates, is_recession_quarter, logical(1))

df_nber <- tibble(
  date = q_dates,
  recession = as.integer(recession_flag)
)

# 4) Now call the recession-probability plot 
#    If want the exact same span as h/s analysis:
# restrict_rng <- if (exists("df_e")) range(df_e$date) else range(df_nber$date)

crp_real <- plot_recession_prob(
  df_nber,
  k_seq   = 12:90,
  x_values= c(3,4,5)
)

print(crp_real$plot)


first_peak_real <- first_peak_k(crp_real$data, group_vars = "x")
# if need "±x" label:
first_peak_real <- first_peak_real %>%
  dplyr::mutate(x_band = paste0("±", x)) %>%
  dplyr::select(x_band, first_peak_k, first_peak_prob)

print(first_peak_real)

# Convenience wrapper to return the *row* at the first peak (≤ k_cap)
peak_rows <- function(df, group_vars, k_col = "k", p_col = "p_hat",
                      k_cap = 56L, k_min = 7L) {
  stopifnot(all(c(group_vars, k_col, p_col) %in% names(df)))
  df %>%
    dplyr::filter(is.finite(.data[[k_col]]),
                  is.finite(.data[[p_col]]),
                  .data[[k_col]] >= k_min, .data[[k_col]] <= k_cap) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::slice_max(order_by = .data[[p_col]], n = 1, with_ties = TRUE) %>%
    dplyr::slice_min(order_by = .data[[k_col]], n = 1, with_ties = TRUE) %>%
    dplyr::ungroup()
}


# Joint tests (Appendix B cases):
# Case 1:  K = {38, 56}, x = 5
jt_boot_real <- joint_test_bootstrap(df_nber, K = c(38,56), x = 5, B = 10000)
jt_hac_real  <- joint_test_hac(df_nber,  K = c(38,56), x = 5)

cat("Case 1 (bootstrap) p-value:", sprintf("%.4f", jt_boot_real$pval), "\n")
cat("Case 1 (HAC)       p-value:", sprintf("%.4f", jt_hac_real$pval),  "\n")

# Case 2:  K = {36,37,38,39,55,56,57,58}, x = 5
jt2_real_boot <- joint_test_bootstrap(df_nber, K = c(36:39,55:58), x=5, B=10000)
jt2_real_hac  <- joint_test_hac(df_nber,  K = c(36:39,55:58), x=5)

cat("Case 2 (bootstrap) p-value:", sprintf("%.4f", jt2_real_boot$pval), "\n")
cat("Case 2 (HAC)       p-value:", sprintf("%.4f", jt2_real_hac$pval),  "\n")













###############################################################################
## Rolling R|R — robust windowing 
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(lubridate)
  library(ggplot2); library(sandwich)
})

# ───────────────────────── Helpers (baseline-consistent) ──────────────────────
make_window_indicator <- function(r, k, x){
  Tn <- length(r); y <- integer(Tn)
  for (t in seq_len(Tn)) {
    lo <- t + k - x; hi <- t + k + x
    if (lo < 1L || hi > Tn) { y[t] <- NA_integer_; next }
    y[t] <- as.integer(any(r[lo:hi] == 1L))
  }
  y
}
default_nw_lag <- function(Tn) floor(4 * (Tn/100)^(2/9))

# HAC var(α+β) but return NA on singular designs
nw_var_alpha_plus_beta_safe <- function(fit, lag = NULL){
  V <- tryCatch(
    if (is.null(lag)) sandwich::NeweyWest(fit) else sandwich::NeweyWest(fit, lag = lag),
    error = function(e) matrix(NA_real_, nrow = length(coef(fit)), ncol = length(coef(fit)))
  )
  e <- c(1,1)
  as.numeric(t(e) %*% V %*% e)
}

# estimator made robust for thin/degenerate subsamples
recession_probability_curve_safe <- function(df, k_seq = 12:90, x_window = 5,
                                             min_gap_from_current = 12L,
                                             restrict_to = NULL, nw_lag = NULL,
                                             min_keep = 30L, min_class_each = 5L){
  stopifnot(all(c("date","recession") %in% names(df)))
  z <- df %>%
    dplyr::arrange(.data$date) %>%
    dplyr::mutate(date = as.Date(.data$date), recession = as.integer(.data$recession))
  
  # robust restriction to [start, end]
  if (!is.null(restrict_to)) {
    rt <- as.Date(as.vector(restrict_to))
    if (length(rt) != 2L || any(!is.finite(rt))) stop("restrict_to must be length-2 Date vector.")
    z <- z %>% dplyr::filter(.data$date >= rt[1] & .data$date <= rt[2])
  }
  
  if (nrow(z) == 0L) {
    return(tibble(k = k_seq, p_hat = NA_real_, se = NA_real_, n = 0L, x = x_window))
  }
  
  r <- as.integer(z$recession); Tn <- length(r)
  if (is.null(nw_lag)) nw_lag <- default_nw_lag(Tn)
  
  out <- lapply(k_seq, function(k){
    if (k < min_gap_from_current)
      return(tibble(k=k, p_hat=NA_real_, se=NA_real_, n=NA_integer_))
    
    y <- make_window_indicator(r, k, x_window)
    keep <- is.finite(y); n_keep <- sum(keep)
    if (n_keep < min_keep)
      return(tibble(k=k, p_hat=NA_real_, se=NA_real_, n=n_keep))
    
    r_keep <- r[keep]; y_keep <- y[keep]
    if (min(sum(r_keep==1L), sum(r_keep==0L)) < min_class_each)
      return(tibble(k=k, p_hat=NA_real_, se=NA_real_, n=n_keep))
    
    fit <- lm(y_keep ~ r_keep)
    cf  <- coef(fit)
    if (!all(is.finite(cf[1:2])))
      return(tibble(k=k, p_hat=NA_real_, se=NA_real_, n=n_keep))
    
    p_hat <- sum(cf[1:2])
    var_p <- nw_var_alpha_plus_beta_safe(fit, lag = nw_lag)
    se    <- if (is.finite(var_p)) sqrt(max(var_p, 0)) else NA_real_
    tibble(k=k, p_hat=p_hat, se=se, n=n_keep)
  })
  
  dplyr::bind_rows(out) %>% dplyr::mutate(x = x_window)
}

# Schemas & binders
empty_curves_schema <- tibble(
  window_start = as.Date(character()),
  window_end   = as.Date(character()),
  window_label = character(),
  k = integer(), p_hat = double(), se = double(), n = integer(), x = integer()
)

empty_bands_schema <- tibble(
  window_start = as.Date(character()),
  window_end   = as.Date(character()),
  window_label = character(),
  k = integer(), x = integer(), lo = double(), hi = double()
)

empty_peaks_schema <- tibble(
  window_start = as.Date(character()),
  window_end   = as.Date(character()),
  window_label = character(),
  x = integer(), first_peak_k = integer(), first_peak_prob = double()
)

safe_bind <- function(lst, schema){
  if (!length(lst)) return(schema)
  out <- tryCatch(dplyr::bind_rows(lst), error = function(e) schema)
  miss <- setdiff(names(schema), names(out))
  if (length(miss)) for (nm in miss) out[[nm]] <- schema[[nm]]
  dplyr::select(out, dplyr::all_of(names(schema)))
}

# Date math
add_quarters <- function(d, n){
  d <- as.Date(d)
  y <- year(d); m <- month(d)
  q <- ((m-1) %/% 3) + 1
  q_new <- q + n
  y_new <- y + (q_new - 1) %/% 4
  q_new <- ((q_new - 1) %% 4) + 1
  m_new <- (q_new - 1) * 3 + 1
  as.Date(sprintf("%04d-%02d-01", y_new, m_new))
}
fmt_q <- function(d) sprintf("%dQ%d", year(as.Date(d)), quarter(as.Date(d)))

# Rolling R|R main runner 
run_rr_rolling_windows <- function(
    df_nber,                 # {date (Date), recession ∈ {0,1}}
    window_years   = 40,     # e.g., 30/40/50
    step_quarters  = 4,      # slide yearly
    k_seq          = 12:90,
    x_values       = c(3,4,5),
    ci_lev         = 0.80
){
  stopifnot(all(c("date","recession") %in% names(df_nber)))
  z <- df_nber %>%
    arrange(date) %>%
    mutate(date = as.Date(date), recession = as.integer(recession))
  stopifnot(inherits(z$date, "Date"))
  
  first_q <- floor_date(min(z$date), "quarter")
  last_q  <- floor_date(max(z$date), "quarter")
  
  starts_all <- seq.Date(from = first_q, to = last_q, by = "quarter")
  starts     <- starts_all[ seq(1, length(starts_all), by = step_quarters) ]
  
  # keep only starts that allow a full window: end = start + (4*Y - 1) quarters
  can_finish <- vapply(starts, function(s){
    e <- add_quarters(s, 4L*window_years - 1L)
    e <= last_q
  }, logical(1))
  starts <- starts[can_finish]
  if (!length(starts)) stop("No valid rolling windows with given window_years and sample.")
  
  all_labels <- vapply(starts, function(s){
    e <- add_quarters(s, 4L*window_years - 1L)
    paste0(fmt_q(s), "–", fmt_q(e))
  }, character(1))
  
  
  zval <- qnorm((1 + ci_lev)/2)
  curves_all <- list(); bands_all <- list(); peaks_all <- list()
  
  for (s in starts){
    e   <- add_quarters(s, 4L*window_years - 1L)     # inclusive last quarter (quarter-start)
    lab <- paste0(fmt_q(s), "–", fmt_q(e))
    
    # curves for each x
    cur_lst <- lapply(x_values, function(xw){
      recession_probability_curve_safe(
        df = z, k_seq = k_seq, x_window = xw, restrict_to = c(s, e)
      ) %>%
        dplyr::mutate(window_start = s, window_end = e, window_label = lab, .before = 1)
    })
    cur_df <- dplyr::bind_rows(cur_lst)
    
    if (!nrow(cur_df) || !any(is.finite(cur_df$p_hat))) {
      # create a placeholder so the facet still appears
      cur_df <- tidyr::expand_grid(k = k_seq, x = x_values) %>%
        dplyr::mutate(p_hat = NA_real_, se = NA_real_, n = NA_integer_) %>%
        dplyr::mutate(window_start = s, window_end = e, window_label = lab, .before = 1) %>%
        dplyr::select(window_start, window_end, window_label, k, p_hat, se, n, x)
    }
    
    # enforce minimal, ordered schema before storing
    cur_df <- cur_df %>%
      dplyr::select(window_start, window_end, window_label, k, p_hat, se, n, x)
    curves_all[[lab]] <- cur_df
    
    # bands at largest x
    x_max <- max(x_values)
    b_df <- cur_df %>%
      dplyr::filter(x == x_max, is.finite(p_hat), is.finite(se)) %>%
      dplyr::transmute(window_start, window_end, window_label,
                       k, x, lo = pmax(0, p_hat - zval*se), hi = pmin(1, p_hat + zval*se))
    bands_all[[lab]] <- b_df
    
    # first peak ≤ 56
    pk <- first_peak_k(cur_df, group_vars = "x") %>%
      dplyr::mutate(window_start = s, window_end = e, window_label = lab,
                    .before = 1)
    peaks_all[[lab]] <- pk
  }
  
  curves <- safe_bind(curves_all, empty_curves_schema)
  bands  <- safe_bind(bands_all,  empty_bands_schema)
  peaks  <- safe_bind(peaks_all,  empty_peaks_schema)
  
  # Ensure facets know about every window, even if all values are NA
  if (nrow(curves)) curves$window_label <- factor(curves$window_label, levels = all_labels)
  if (nrow(bands))  bands$window_label  <- factor(bands$window_label,  levels = all_labels)
  if (nrow(peaks))  peaks$window_label  <- factor(peaks$window_label,  levels = all_labels)
  
  
  # Plots
  humps <- data.frame(xmin = c(36,55), xmax = c(40,58))
  
  p_by_x <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = humps, inherit.aes = FALSE,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "red", alpha = 0.10
    )
  
  if (nrow(bands)) {
    p_by_x <- p_by_x +
      ggplot2::geom_ribbon(
        data = bands,
        ggplot2::aes(k, ymin = lo, ymax = hi, group = window_label),
        fill = "grey85"
      )
  }
  if (nrow(curves)) {
    sub_txt <- paste0(
      "Lines: x ∈ {", paste(sort(unique(curves$x)), collapse = ", "),
      "}; shaded: ", 100*ci_lev, "% HAC CI at largest x"
    )
    p_by_x <- p_by_x +
      ggplot2::geom_line(
        data = curves,
        ggplot2::aes(k, p_hat, color = factor(x)),
        linewidth = 0.9
      ) +
     ggplot2::facet_wrap(~ window_label, ncol = 2, drop = FALSE) +
      ggplot2::scale_color_manual(
        values = setNames(c("#1b9e77","#d95f02","#7570b3")[seq_along(sort(unique(curves$x)))],
                          sort(unique(curves$x))),
        name = "Window ±x"
      ) +
      ggplot2::coord_cartesian(ylim = c(0,1)) +
      ggplot2::labs(
        title = "Rolling R|R: conditional probability curves by time window",
        subtitle = sub_txt,
        x = "k (quarters since recession)", y = "Probability"
      ) +
      ggplot2::theme_minimal(base_size = 12)
  }
  
  p_bigx <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = humps, inherit.aes = FALSE,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "red", alpha = 0.10
    )
  
  if (nrow(bands)) {
    p_bigx <- p_bigx +
      ggplot2::geom_ribbon(
        data = bands,
        ggplot2::aes(k, ymin = lo, ymax = hi, group = window_label),
        fill = "grey85"
      ) +
      ggplot2::geom_line(
        data = bands,
        ggplot2::aes(k, (lo + hi)/2, color = window_label),
        linewidth = 0.9
      ) +
      ggplot2::facet_wrap(~ window_label, ncol = 2, drop = FALSE) +
      ggplot2::coord_cartesian(ylim = c(0,1)) +
      ggplot2::labs(
        title = "Rolling R|R at largest x (HAC CI shaded)",
        x = "k (quarters since recession)", y = "Probability", color = NULL
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(legend.position = "none")
  }
  
  list(curves = curves, bands = bands, peaks = peaks,
       plot_by_x = p_by_x, plot_bigx = p_bigx,
       meta = list(window_years = window_years, step_quarters = step_quarters,
                   k_seq = k_seq, x_values = x_values, ci_lev = ci_lev))
}

### CALL
#for disjoint: 
# window_years <- 40
#and put
# step_quarters = 4*window_years
out_rr_roll <- run_rr_rolling_windows(
  df_nber,
  window_years  = 40,
  step_quarters = 4,
  k_seq         = 12:90,
  x_values      = c(3,4,5),
  ci_lev        = 0.80
)
print(out_rr_roll$plot_by_x)
print(out_rr_roll$plot_bigx)
out_rr_roll$peaks




###############################################################################
## Helpers for diagnostics
###############################################################################

# counts for a given window/horizon
counts_for_window_k <- function(df_nber, start_date, end_date, k, x){
  z <- df_nber %>%
    dplyr::filter(date >= as.Date(start_date), date <= as.Date(end_date)) %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(recession = as.integer(recession))
  r <- z$recession
  if (!length(r)) return(list(n_keep = 0L, n_R1 = 0L))
  y <- make_window_indicator(r, k, x)
  keep <- is.finite(y)
  n_keep <- sum(keep)
  n_R1   <- if (n_keep > 0) sum(r[keep] == 1L) else 0L
  list(n_keep = n_keep, n_R1 = n_R1)
}

# main snapshot
snapshot_info_diagnostics <- function(out_roll,
                                      df_nber,
                                      k_focus = c(38L,56L),
                                      min_keep = 30L,
                                      x_focus  = max(out_roll$meta$x_values)){
  stopifnot(all(c("curves","peaks","meta") %in% names(out_roll)))
  cv   <- out_roll$curves
  pkdf <- out_roll$peaks %>% dplyr::filter(x == x_focus)
  stopifnot(nrow(cv) > 0)
  
  # distinct windows
  wins <- cv %>%
    dplyr::distinct(window_start, window_end, window_label) %>%
    dplyr::arrange(window_start)
  
  # counts and curve values at k_focus for x_focus
  per_win <- purrr::pmap_dfr(wins, function(window_start, window_end, window_label){
    # p-hat and n_keep from curves object (already computed) at x_focus
    kk <- cv %>% dplyr::filter(window_label == !!window_label, x == !!x_focus,
                               k %in% k_focus)
    # ensure both ks present (fill with NA otherwise)
    kk <- tidyr::complete(kk, k = k_focus, fill = list(p_hat = NA_real_, n = NA_integer_))
    
    # compute n_R1 by recomputing counts (not stored in curves)
    cnts <- purrr::map(k_focus, ~counts_for_window_k(df_nber, window_start, window_end, .x, x_focus))
    nR1  <- purrr::map_int(cnts, "n_R1")
    
    # feasible k_max given min_keep (based on curves' n)
    k_grid <- cv %>% dplyr::filter(window_label == !!window_label, x == !!x_focus) %>%
      dplyr::arrange(k) %>% dplyr::select(k, n)
    feasible_kmax <- suppressWarnings(max(k_grid$k[is.finite(k_grid$n) & k_grid$n >= min_keep], na.rm = TRUE))
    if (!is.finite(feasible_kmax)) feasible_kmax <- NA_integer_
    
    # actual designed k_max from meta
    designed_kmax <- max(out_roll$meta$k_seq)
    cap_recommended <- is.finite(feasible_kmax) && feasible_kmax < designed_kmax
    
    # first-peak info (x = x_focus)
    peak_row <- pkdf %>% dplyr::filter(window_label == !!window_label) %>% dplyr::slice_tail(n = 1)
    k_star <- if (nrow(peak_row)) peak_row$first_peak_k else NA_integer_
    p_star <- if (nrow(peak_row)) peak_row$first_peak_prob else NA_real_
    
    # p at 38 and 56; Δ = p38 - p56
    p38 <- kk$p_hat[match(38L, kk$k)]
    p56 <- kk$p_hat[match(56L, kk$k)]
    delta <- if (is.finite(p38) && is.finite(p56)) p38 - p56 else NA_real_
    
    tibble::tibble(
      window_label = window_label,
      window_start = window_start,
      window_end   = window_end,
      x_focus      = x_focus,
      k_star       = k_star,
      p_star       = p_star,
      p38          = p38,
      p56          = p56,
      delta_p38_56 = delta,
      n_keep_38    = kk$n[match(38L, kk$k)],
      n_R1_38      = nR1[which(k_focus == 38L)],
      n_keep_56    = kk$n[match(56L, kk$k)],
      n_R1_56      = nR1[which(k_focus == 56L)],
      feasible_kmax   = feasible_kmax,
      designed_kmax   = designed_kmax,
      cap_recommended = cap_recommended
    )
  })
  
  # Nice ordering
  per_win <- per_win %>% dplyr::arrange(window_start)
  
  # Simple, compact plot: kept observations at 38/56 and peak–trough contrast
  p_counts <- per_win %>%
    dplyr::select(window_label, n_keep_38, n_keep_56) %>%
    tidyr::pivot_longer(cols = c(n_keep_38, n_keep_56),
                        names_to = "which_k", values_to = "n_keep") %>%
    dplyr::mutate(which_k = dplyr::recode(which_k,
                                          n_keep_38 = "n_keep @ k=38",
                                          n_keep_56 = "n_keep @ k=56")) %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(window_label, n_keep),
                                 y = n_keep, fill = which_k)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Information per rolling window",
                  x = "Window", y = "Kept observations",
                  fill = NULL) +
    ggplot2::theme_minimal(base_size = 11)
  
  p_delta <- per_win %>%
    ggplot2::ggplot(ggplot2::aes(x = window_start, y = delta_p38_56)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_line() +
    ggplot2::labs(title = expression(Delta~"= p"[38] - " p"[56]),
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11)
  
  list(table = per_win, plot_counts = p_counts, plot_delta = p_delta)
}

### CALL
diag <- snapshot_info_diagnostics(out_rr_roll, df_nber,
                                  k_focus = c(38L,56L),
                                  min_keep = 30L,
                                  x_focus  = max(out_rr_roll$meta$x_values))


print(diag$table)
print(diag$table, n=21)
print(diag$plot_delta)











## PREPARATION - Financial series to 2018Q4 (inclusive)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(lubridate); library(tibble); library(purrr)
})

# IO helper
read_stress_local <- function(file, value_col = NULL) {
  x  <- readr::read_csv(file, show_col_types = FALSE)
  nm <- tolower(names(x)); names(x) <- nm
  date_col <- if ("date" %in% nm) "date" else nm[1]
  x <- x %>% mutate(!!date_col := as.Date(.data[[date_col]]))
  
  if (is.null(value_col)) {
    cand <- setdiff(names(x), date_col)
    cand_num <- cand[sapply(x[cand], is.numeric)]
    if (!length(cand_num)) stop("No numeric value column in: ", file)
    value_col <- cand_num[1]
  }
  
  x %>%
    transmute(date = as.Date(.data[[date_col]]),
              value = suppressWarnings(as.numeric(.data[[value_col]]))) %>%
    arrange(.data$date) %>%
    dplyr::filter(is.finite(.data$value))
}


# Native=>quarterly binary S_t via threshold q and min_run at native freq
stress_to_quarterly_indicator <- function(df_native, threshold = NULL, q = 0.95,
                                          min_run = 2L, q_grid = NULL) {
  stopifnot(all(c("date","value") %in% names(df_native)))
  z <- df_native %>% arrange(.data$date)
  thr  <- if (is.null(threshold)) stats::quantile(z$value, probs = q, na.rm = TRUE) else threshold
  flag <- as.integer(z$value >= thr)
  
  r <- rle(flag); vals <- r$values; lens <- r$lengths
  for (i in seq_along(vals)) if (vals[i] == 1L && lens[i] < min_run) vals[i] <- 0L
  flag2 <- inverse.rle(list(values = vals, lengths = lens))
  
  if (is.null(q_grid)) {
    q_grid <- seq(floor_date(min(z$date), "quarter"),
                  floor_date(max(z$date), "quarter"), by = "quarter")
  }
  q_grid <- as.Date(q_grid)
  
  stress_q <- vapply(q_grid, function(qd) {
    q_end <- qd %m+% months(3) - days(1)
    any(flag2[z$date >= qd & z$date <= q_end] == 1L)
  }, logical(1))
  
  tibble(date = q_grid, stress_high = as.integer(stress_q))
}


# files
stress_files <- tibble::tribble(
  ~file,                     ~label,
  "NFCI_weekly.csv",         "NFCI",
  "ANFCI_weekly.csv",        "ANFCI",
  "NFCIRISK_weekly.csv",     "NFCI RISK",
  "NFCICREDIT_weekly.csv",   "NFCI CREDIT",
  "NFCILEVERAGE_weekly.csv", "NFCI LEVERAGE",
)

# A) Load NATIVE frequency once 
financial_native <- purrr::map_dfr(seq_len(nrow(stress_files)), function(i){
  read_stress_local(stress_files$file[i]) %>%
    mutate(source = stress_files$label[i], freq = "native")
}) %>%
  arrange(source, date)

# B) Quarterly means for PLOTTING ONLY (not for thresholds)
financial_q <- financial_native %>%
  mutate(qdate = lubridate::floor_date(date, "quarter")) %>%
  group_by(source, qdate) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  transmute(source, date = qdate, value)

# C) Build S_t with NATIVE threshold + min_run (weeks), 
#              then aggregate to quarterly date grid 
stopifnot(exists("df_nber"), all(c("date","recession") %in% names(df_nber)))
q_grid <- as.Date(df_nber$date)  # quarterly grid used everywhere downstream

qS_all_native <- financial_native %>%
  group_split(source, .keep = TRUE) %>%
  map_dfr(function(dfw){
    lab <- unique(dfw$source)
    stress_to_quarterly_indicator(
      df_native = dfw %>% dplyr::select(date, value),
      q        = 0.95,   # same quantile
      min_run  = 2L,     # interpreted in NATIVE units (weeks)
      q_grid   = q_grid
    ) %>% mutate(source = lab)
  }) %>%
  transmute(source = as.character(source),
            date   = as.Date(date),
            stress_high = as.integer(stress_high))

#  D) Convenience alias for downstream code
qS_all <- qS_all_native




################################################################################
## R|S (NBER recession | high financial stress)
###############################################################################


stopifnot(exists("df_nber"), all(c("date","recession") %in% names(df_nber)))
df_nber <- df_nber %>%
  arrange(.data$date) %>%
  mutate(date = as.Date(.data$date), recession = as.integer(.data$recession))

# 1) IO + PREP HELPERS
read_stress_local <- function(file, value_col = NULL) {
  x  <- readr::read_csv(file, show_col_types = FALSE)
  nm <- tolower(names(x)); names(x) <- nm
  date_col <- if ("date" %in% nm) "date" else nm[1]
  x <- x %>% mutate(!!date_col := as.Date(.data[[date_col]]))
  
  if (is.null(value_col)) {
    cand <- setdiff(names(x), date_col)
    cand_num <- cand[sapply(x[cand], is.numeric)]
    if (!length(cand_num)) stop("No numeric value column in: ", file)
    value_col <- cand_num[1]
  }
  
  x %>%
    transmute(date = as.Date(.data[[date_col]]),
              value = suppressWarnings(as.numeric(.data[[value_col]]))) %>%
    arrange(.data$date) %>%
    dplyr::filter(is.finite(.data$value))
}

# 2) FORWARD-WINDOW + VARIANCE + LAGS
fw_indicator <- function(event01, k, x){
  Tn <- length(event01); y <- integer(Tn)
  for (t in seq_len(Tn)) {
    lo <- t + k - x; hi <- t + k + x
    if (lo < 1L || hi > Tn) { y[t] <- NA_integer_; next }
    y[t] <- as.integer(any(event01[lo:hi] == 1L))
  }
  y
}

default_nw_lag <- function(Tn) floor(4 * (Tn/100)^(2/9))

var_alpha_plus_beta_hac <- function(fit, lag = NULL){
  V <- if (is.null(lag)) sandwich::NeweyWest(fit, prewhite = FALSE, adjust = TRUE)
  else               sandwich::NeweyWest(fit, lag = lag, prewhite = FALSE, adjust = TRUE)
  e <- rep(1, length(coef(fit)))
  as.numeric(t(e) %*% V %*% e)
}

# Single-source curve: p_hat(k,x) = α+β from lm(y ~ S) with HAC se(α+β)
rs_curve_one <- function(df_event, df_cond, k_seq = 12:90, x = 5,
                         min_gap = 12L, min_keep = 30L, nw_lag = NULL) {
  a <- df_event %>% arrange(.data$date) %>% mutate(date = as.Date(.data$date))
  b <- df_cond  %>% arrange(.data$date) %>% mutate(date = as.Date(.data$date))
  inter <- intersect(a$date, b$date)
  if (length(inter) < 2L) return(tibble(k=integer(), x=integer(), p_hat=double(), se=double(), n=integer()))
  rng <- range(inter)
  a <- dplyr::filter(a, .data$date >= rng[1] & .data$date <= rng[2])
  b <- dplyr::filter(b, .data$date >= rng[1] & .data$date <= rng[2])
  stopifnot(nrow(a) == nrow(b), all(a$date == b$date))
  
  R <- as.integer(a$event); S <- as.integer(b$cond); Tn <- length(R)
  if (is.null(nw_lag)) nw_lag <- default_nw_lag(Tn)
  
  map_dfr(k_seq, function(k){
    if (k < min_gap) return(tibble(k=k, x=x, p_hat=NA_real_, se=NA_real_, n=NA_integer_))
    y <- fw_indicator(R, k, x); keep <- is.finite(y)
    n_keep <- sum(keep)
    if (n_keep < min_keep) return(tibble(k=k, x=x, p_hat=NA_real_, se=NA_real_, n=n_keep))
    d <- tibble(y = y[keep], S = S[keep])
    if (min(sum(d$S==1L), sum(d$S==0L)) < 5) return(tibble(k=k, x=x, p_hat=NA_real_, se=NA_real_, n=n_keep))
    fit   <- lm(y ~ S, data = d)
    p_hat <- sum(coef(fit), na.rm = TRUE)
    se    <- sqrt(pmax(var_alpha_plus_beta_hac(fit, nw_lag), 0))
    tibble(k=k, x=x, p_hat=p_hat, se=se, n=n_keep)
  })
}

# 3) JOINT TESTS (unchanged methods): HAC and episode bootstrap on S
find_S_starts <- function(S){
  up <- c(0L, diff(S))
  idx <- which(up == 1L)
  if (length(idx) == 0L && S[1] == 1L) idx <- 1L
  idx
}

resample_S_episodes <- function(S, T_target){
  idx <- find_S_starts(S); if (!length(idx)) stop("No S=1 episodes.")
  segs <- lapply(seq_along(idx), function(i){
    a <- idx[i]; b <- if (i < length(idx)) idx[i+1]-1L else length(S)
    S[a:b]
  })
  out <- integer(0)
  while (length(out) < T_target) out <- c(out, segs[[sample.int(length(segs),1)]])
  out[seq_len(T_target)]
}

rs_joint_tests <- function(df_RS, K, x = 5, L = NULL, B = 1000, seed = 123){
  z <- df_RS %>% arrange(.data$date)
  R <- as.integer(z$event); S <- as.integer(z$crisis); Tn <- length(R)
  
  # compute p_obs want (k-specific keep is fine for reporting)
  p_obs <- sapply(K, function(k){
    y <- fw_indicator(R, k, x); keep <- is.finite(y) & (S==1L)
    if (!any(keep)) NA_real_ else mean(y[keep])
  })
  if (any(!is.finite(p_obs))) stop("Not enough data for chosen K.")
  
  # HAC on common sample across K
  Y <- sapply(K, function(k) fw_indicator(R, k, x))
  keep <- apply(is.finite(Y), 1, all)
  Yk <- Y[keep, , drop=FALSE]; Sk <- S[keep]
  
  denom <- mean(Sk); if (!is.finite(denom) || denom <= 0) stop("E[S]=0.")
  Z  <- (Sk * Yk) / denom
  # >>> use common-sample mean for HAC test
  p_obs_hac <- colMeans(Z)
  
  Zc <- scale(Z, center = TRUE, scale = FALSE); T0 <- nrow(Zc)
  if (is.null(L)) L <- default_nw_lag(T0)  # <<< use T0, not Tn
  Gamma0 <- crossprod(Zc)/T0
  w <- 1 - (1:L)/(L+1); Gsum <- matrix(0, ncol(Zc), ncol(Zc))
  for (l in 1:L){
    A <- t(Zc[(l+1):T0, , drop=FALSE]) %*% Zc[1:(T0-l), , drop=FALSE] / T0
    Gsum <- Gsum + w[l]*(A + t(A))
  }
  # >>> divide by T0 to get Var(mean)
  Sigma_hac <- (Gamma0 + Gsum) / T0
  
  pbar <- mean(p_obs_hac); e <- rep(1, length(K))
  ev <- eigen(Sigma_hac, symmetric=TRUE); lam <- pmax(ev$values, .Machine$double.eps)
  Sminushalf <- ev$vectors %*% diag(1/sqrt(lam)) %*% t(ev$vectors)
  eta <- as.numeric(Sminushalf %*% (p_obs_hac - pbar*e))
  stat_hac <- sum(eta^2); pval_hac <- 1 - pchisq(stat_hac, df = length(K)-1L)
  
  # Episode bootstrap on S
  set.seed(seed)
  Pmat <- matrix(NA_real_, nrow=B, ncol=length(K))
  for (b in 1:B){
    Sb <- resample_S_episodes(S, Tn)
    Pmat[b,] <- sapply(K, function(k){
      y <- fw_indicator(R, k, x); keep <- is.finite(y) & (Sb==1L)
      if (!any(keep)) NA_real_ else mean(y[keep])
    })
  }
  Sigma_boot <- stats::cov(Pmat, use="pairwise.complete.obs")
  evb <- eigen(Sigma_boot, symmetric=TRUE); lamb <- pmax(evb$values, .Machine$double.eps)
  Smb <- evb$vectors %*% diag(1/sqrt(lamb)) %*% t(evb$vectors)
  etab <- as.numeric(Smb %*% (p_obs - pbar*e))
  stat_boot <- sum(etab^2); pval_boot <- 1 - pchisq(stat_boot, df = length(K)-1L)
  
  list(p_obs=p_obs, Sigma_hac=Sigma_hac, Sigma_boot=Sigma_boot,
       stat_hac=stat_hac, pval_hac=pval_hac,
       stat_boot=stat_boot, pval_boot=pval_boot)
}

# 4) MASTER RUNNER (build S, curves, bands, peaks, tests, plots)
run_r_given_stress_static <- function(
    df_nber,
    stress_files = NULL,
    qS_all = NULL,
    q = 0.95, min_run = 2L,
    k_seq = 12:90, x_values = c(3,4,5),
    K_case1 = c(38,56), K_case2 = c(36:39,55:58),
    B = 2000, ci_lev = 0.80
){
  # Event series
  event_q <- df_nber %>%
    dplyr::arrange(.data$date) %>%
    dplyr::transmute(date = as.Date(.data$date), event = as.integer(.data$recession))
  
  # Build qS_all if not supplied
  # Build qS_all if not supplied — USE pre-trimmed financial_q
  if (is.null(qS_all)) {
    stopifnot(exists("financial_q"),
              all(c("date","value","source") %in% names(financial_q)))
    q_grid <- as.Date(df_nber$date)
    
    # Apply fixed quantile on the (already quarterly) native series
    qS_all <- financial_q %>%
      arrange(source, date) %>%
      dplyr::group_split(source, .keep = TRUE) %>%
      purrr::map_dfr(function(dfq){
        lab <- unique(dfq$source)
        stress_to_quarterly_indicator(
          df_native = dfq %>% dplyr::select(date, value),
          q        = q,
          min_run  = min_run,
          q_grid   = q_grid
        ) %>% dplyr::mutate(source = lab)
      })
  }
  # normalize to {source,date,cond}
  qS_all <- qS_all %>%
    dplyr::transmute(source = as.character(.data$source),
                     date   = as.Date(.data$date),
                     cond   = as.integer(.data$stress_high)) %>%
    dplyr::arrange(.data$source, .data$date)
  
  
  # CURVES: split-map-bind (avoid group_modify)
  curves <- qS_all %>%
    dplyr::select(.data$source, .data$date, .data$cond) %>%
    split(.$source) %>%
    purrr::imap_dfr(function(dfS, src){
      df_cond <- dplyr::select(dfS, .data$date, .data$cond)
      inter <- intersect(df_cond$date, event_q$date)
      if (length(inter) < 2L)
        return(tibble::tibble(k=integer(), x=integer(), p_hat=double(), se=double(), n=integer(), source=character())[0,])
      out <- purrr::map_dfr(x_values, ~ rs_curve_one(event_q, df_cond, k_seq = k_seq, x = .x))
      if (!nrow(out))
        return(tibble::tibble(k=integer(), x=integer(), p_hat=double(), se=double(), n=integer(), source=character())[0,])
      dplyr::mutate(out, source = src)
    })
  
  # Bands (largest x) + “first peak ≤ 56”
  stopifnot(all(c("k","x","p_hat","se","n","source") %in% names(curves)))
  x_max <- max(x_values); z <- qnorm((1 + ci_lev)/2)
  
  bands_df <- curves %>%
    dplyr::filter(.data$x == x_max, is.finite(.data$p_hat), is.finite(.data$se)) %>%
    dplyr::mutate(lo = pmax(0, .data$p_hat - z*.data$se),
                  hi = pmin(1, .data$p_hat + z*.data$se))
  
  peaks <- first_peak_k(curves, group_vars = c("source","x"))
  
  # Joint tests per source at x = x_max
  jt_tbl <- purrr::map_dfr(unique(qS_all$source), function(src){
    df_S <- qS_all %>% dplyr::filter(.data$source == src) %>%
      dplyr::transmute(date = .data$date, crisis = .data$cond)
    df_RS <- dplyr::inner_join(df_S, event_q, by = "date") %>%
      dplyr::select(.data$date, .data$crisis, .data$event)
    if (sum(df_RS$crisis, na.rm = TRUE) < 8)
      return(tibble::tibble(source = src,
                            case1_boot_pval = NA_real_, case1_hac_pval = NA_real_,
                            case2_boot_pval = NA_real_, case2_hac_pval = NA_real_))
    jt1 <- rs_joint_tests(df_RS, K = K_case1, x = x_max, B = B)
    jt2 <- rs_joint_tests(df_RS, K = K_case2, x = x_max, B = B)
    tibble::tibble(source = src,
                   case1_boot_pval = jt1$pval_boot, case1_hac_pval = jt1$pval_hac,
                   case2_boot_pval = jt2$pval_boot, case2_hac_pval = jt2$pval_hac)
  })
  
  # Plots
  humps <- data.frame(xmin = c(36,55), xmax = c(40,58))
  
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = humps, inherit.aes = FALSE,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                       fill = "red", alpha = 0.10) +
    ggplot2::geom_ribbon(data = bands_df, ggplot2::aes(k, ymin = lo, ymax = hi, group = source),
                         fill = "grey85") +
    ggplot2::geom_line(data = curves, ggplot2::aes(k, p_hat, color = factor(x, levels = x_values)),
                       linewidth = 0.9) +
    ggplot2::facet_wrap(~ source, ncol = 2) +
    ggplot2::scale_color_manual(
      values = setNames(c("#1b9e77","#d95f02","#7570b3")[seq_along(x_values)], x_values),
      name = "Window ±x"
    ) +
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    ggplot2::labs(
      title = "R|S: forward conditional probability by stress source",
      subtitle = paste0("Lines = x∈{", paste(x_values, collapse=", "), "}; shaded = ",
                        100*ci_lev,"% HAC CI at x = ±", x_max),
      x = "k (quarters ahead)", y = "Conditional probability"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  
  
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = humps, inherit.aes = FALSE,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                       fill = "red", alpha = 0.10) +
    ggplot2::geom_ribbon(data = bands_df, ggplot2::aes(k, ymin = lo, ymax = hi, group = source),
                         fill = "grey85") +
    ggplot2::geom_line(data = bands_df, ggplot2::aes(k, p_hat, color = source),
                       linewidth = 0.9) +
    ggplot2::facet_wrap(~ source, ncol = 2) +
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    ggplot2::labs(
      title = paste0("R|S at x = ±", x_max, " (HAC CI shaded)"),
      x = "k (quarters ahead)", y = "Conditional probability", color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  
  
  # Time-series panels (native + quarterly with S and NBER shading)
  # Time-series panels (use pre-trimmed financial_q; no CSV reads) 
  stopifnot(exists("financial_q"),
            all(c("date","value","source") %in% names(financial_q)))
  
  native_all <- financial_q %>%
    dplyr::arrange(source, date) %>%
    dplyr::mutate(freq = "quarterly")
  
  p_native <- if (nrow(native_all)) {
    ggplot2::ggplot(native_all, ggplot2::aes(date, value, color = source)) +
      ggplot2::geom_line(linewidth = 0.7, alpha = 0.9, show.legend = FALSE) +
      ggplot2::facet_wrap(~ source, scales = "free_y", ncol = 1) +
      ggplot2::labs(title = "Financial stress (quarterly)", x = NULL, y = "Index level") +
      ggplot2::theme_minimal(base_size = 12)
  } else NULL
  
  quarterly_levels <- native_all %>%
    dplyr::select(source, date, value) %>%
    dplyr::arrange(source, date)
  
  
  qS_all_plot <- qS_all %>% dplyr::transmute(source=.data$source, date=.data$date, stress_high=.data$cond)
  
  quarterly_merged <- quarterly_levels %>%
    dplyr::left_join(qS_all_plot, by = c("date","source")) %>%
    dplyr::mutate(stress_high = ifelse(is.na(.data$stress_high), 0L, .data$stress_high))
  
  episodes_from_indicator <- function(df_q) {
    df_q <- df_q %>% dplyr::arrange(.data$date)
    z <- df_q$stress_high
    if (!length(z) || all(z == 0L | is.na(z))) return(tibble::tibble(source=df_q$source[1], start=as.Date(NA), end=as.Date(NA))[0,])
    r <- rle(z); ends <- cumsum(r$lengths); starts <- c(1, head(ends, -1) + 1); sel <- which(r$values == 1L)
    tibble::tibble(source = df_q$source[1],
                   start  = df_q$date[starts[sel]],
                   end    = df_q$date[ends[sel]] + months(3) - days(1))
  }
  
  shades <- if (nrow(quarterly_merged)) {
    quarterly_merged %>% dplyr::group_by(.data$source) %>% dplyr::group_modify(~episodes_from_indicator(.x)) %>% dplyr::ungroup()
  } else tibble::tibble(source=character(), start=as.Date(character()), end=as.Date(character()))[0,]
  
  p_quarterly <- if (nrow(quarterly_merged)) {
    ggplot2::ggplot() +
      ggplot2::geom_rect(data = shades, ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = -Inf, ymax = Inf),
                         inherit.aes = FALSE, fill = "grey0", alpha = 0.35) +
      ggplot2::geom_line(data = quarterly_merged, ggplot2::aes(.data$date, .data$value, color = .data$source),
                         linewidth = 0.9, show.legend = FALSE) +
      ggplot2::facet_wrap(~ source, scales = "free_y", ncol = 1) +
      ggplot2::labs(title = "Financial stress (quarterly mean) with high-stress quarters shaded",
                    subtitle = paste0("Shaded = S_t = 1 (q=", q, ", min_run=", min_run, ")"),
                    x = NULL, y = "Quarterly mean level") +
      ggplot2::theme_minimal(base_size = 12)
  } else NULL
  
  nber_shades <- df_nber %>%
    dplyr::mutate(start = .data$date, end = .data$date %m+% months(3) - days(1)) %>%
    dplyr::filter(.data$recession == 1L) %>% dplyr::transmute(start, end)
  
  p_quarterly_nber <- if (nrow(quarterly_merged)) {
    ggplot2::ggplot() +
      ggplot2::geom_rect(data = nber_shades, ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = -Inf, ymax = Inf),
                         inherit.aes = FALSE, fill = "grey40", alpha = 0.20) +
      ggplot2::geom_rect(data = shades, ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = -Inf, ymax = Inf),
                         inherit.aes = FALSE, fill = "grey0", alpha = 0.35) +
      ggplot2::geom_line(data = quarterly_merged, ggplot2::aes(.data$date, .data$value, color = .data$source),
                         linewidth = 0.9, show.legend = FALSE) +
      ggplot2::facet_wrap(~ source, scales = "free_y", ncol = 1) +
      ggplot2::labs(title = "Financial stress (quarterly mean) with high-stress (light) and NBER recessions (dark)",
                    x = NULL, y = "Quarterly mean level") +
      ggplot2::theme_minimal(base_size = 12)
  } else NULL
  
  list(
    curves = curves, bands = bands_df, peaks = peaks, tests = jt_tbl,
    plot_by_x = p1, plot_bigx = p2,
    plot_native = p_native, plot_quarterly = p_quarterly, plot_quarterly_nber = p_quarterly_nber,
    qS_all = qS_all
  )
}


#### 5) CALL
out <- run_r_given_stress_static(df_nber, stress_files = stress_files,
                               qS_all = qS_all,
                               q = 0.9, min_run = 2L,
                               k_seq = 12:90, x_values = c(3,4,5),
                               K_case1 = c(38,56), K_case2 = c(36:39,55:58),
                               B = 10, ci_lev = 0.80)


print(out$plot_bigx)
print(out$plot_by_x)
print(out$plot_native)
print(out$plot_quarterly)
print(out$plot_quarterly_nber)
print(out$peaks)
print(out$tests)








# COMPARISON TABLE: max-prob quarter by x

# 1) NBER peaks
nber_peaks <- first_peak_real %>%
  dplyr::mutate(source = "NBER (R|R)") %>%
  dplyr::select(source, x_band, first_peak_k, first_peak_prob)

# 2) Financial peaks: coerce x -> "±x" to match NBER's x_band
fin_peaks <- out$peaks %>%
  dplyr::mutate(x_band = paste0("±", x)) %>%
  dplyr::select(source, x_band, first_peak_k, first_peak_prob)

# 3) Union of x-band labels to control column order
x_levels <- paste0("±", sort(unique(c(
  as.numeric(gsub("±","", nber_peaks$x_band)),
  suppressWarnings(as.numeric(out$peaks$x))
))))

# 4) Assemble wide comparison table
combined_peaks <- bind_rows(nber_peaks, fin_peaks) %>%
  mutate(
    x_band = factor(x_band, levels = x_levels),
    cell   = sprintf("k = %d | p = %.3f", first_peak_k, first_peak_prob)
  ) %>%
  dplyr::select(source, x_band, cell) %>%
  pivot_wider(names_from = x_band, values_from = cell) %>%
  arrange(source)

# Pretty table
knitr::kable(combined_peaks, align = c("l", rep("c", length(x_levels))),
             caption = "Peak quarter (k) and probability by x, NBER vs. stress indexes")














###############################################################################
## R|S conditional probabilities - rolling stress + joint R-S bootstrap
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(readr)
  library(ggplot2); library(sandwich); library(lubridate); library(slider)
})

`%||%` <- function(a,b) if (length(a)==0 || is.na(a)) b else a
quarter_end <- function(d) lubridate::ceiling_date(d, "quarter") - lubridate::days(1)

# IO: read one native-frequency stress series 
read_stress_local <- function(file, value_col = NULL) {
  x  <- readr::read_csv(file, show_col_types = FALSE)
  nm <- tolower(names(x)); names(x) <- nm
  date_col <- if ("date" %in% nm) "date" else nm[1]
  x <- x %>% mutate(!!date_col := as.Date(.data[[date_col]]))
  
  if (is.null(value_col)) {
    cand <- setdiff(names(x), date_col)
    cand_num <- cand[sapply(x[cand], is.numeric)]
    if (!length(cand_num)) stop("No numeric value column in: ", file)
    value_col <- cand_num[1]
  }
  
  x %>%
    transmute(date = as.Date(.data[[date_col]]),
              value = suppressWarnings(as.numeric(.data[[value_col]]))) %>%
    arrange(.data$date) %>%
    dplyr::filter(is.finite(.data$value))
}


# Stress indexes peaks from curves:
peaks_stress <- first_peak_k(out$curves, group_vars = c("source","x")) %>%
  dplyr::transmute(
    series     = source,
    x          = as.integer(x),
    peak_k     = as.integer(first_peak_k),
    peak_p     = as.numeric(first_peak_prob)
  )

# NBER R|R peaks from crp_real:
peaks_nber <- first_peak_k(crp_real$data, group_vars = "x") %>%
  dplyr::transmute(
    series = "NBER (R|R)",
    x      = as.integer(x),
    peak_k = as.integer(first_peak_k),
    peak_p = as.numeric(first_peak_prob)
  )


# Rolling-quantile stress -> quarterly S_t ∈ {0,1} 
stress_to_quarterly_indicator_rolling <- function(df_native,
                                                  q = 0.9,
                                                  window_years = 12,
                                                  min_run = 2L,
                                                  quarter_rule = c("any","majority"),
                                                  q_grid = NULL) {
  quarter_rule <- match.arg(quarter_rule)
  z <- df_native %>% arrange(.data$date) %>% mutate(date = as.Date(.data$date))
  
  dd <- as.numeric(median(diff(z$date)))
  if (!is.finite(dd) || dd <= 0) stop("Bad/irregular date spacing in native series.")
  n_per_year <- round(365.25 / dd)
  W <- max(8L, floor(window_years * n_per_year))
  
  q_roll <- slider::slide_dbl(
    z$value,
    .f = ~ stats::quantile(.x, probs = q, na.rm = TRUE, names = FALSE),
    .before = W - 1, .complete = FALSE
  )
  S_nat <- as.integer(z$value >= q_roll)
  
  # collapse short spikes
  r <- rle(S_nat); vals <- r$values; lens <- r$lengths
  for (i in seq_along(vals)) if (vals[i] == 1L && lens[i] < min_run) vals[i] <- 0L
  S_nat2 <- inverse.rle(list(values = vals, lengths = lens))
  
  # quarter grid
  if (is.null(q_grid)) {
    q_grid <- seq(lubridate::floor_date(min(z$date), "quarter"),
                  lubridate::floor_date(max(z$date), "quarter"), by = "quarter")
  }
  q_grid <- as.Date(q_grid)
  
  # aggregate to quarter
  S_q <- vapply(q_grid, function(qd){
    q_end <- quarter_end(qd)
    idx <- which(z$date >= qd & z$date <= q_end)
    if (!length(idx)) return(0L)
    if (quarter_rule == "any") as.integer(any(S_nat2[idx] == 1L))
    else                        as.integer(mean(S_nat2[idx]) >= 0.5)
  }, integer(1))
  
  tibble(date = q_grid, stress_high = S_q)
}

# optional diagnostic for rolling threshold 
stress_rolling_diag_plot <- function(df_native, q = 0.95, window_years = 12, min_run = 2L) {
  z <- df_native %>% arrange(.data$date) %>% mutate(date = as.Date(.data$date))
  dd <- as.numeric(median(diff(z$date))); n_per_year <- round(365.25/dd)
  W <- max(8L, floor(window_years * n_per_year))
  
  q_roll <- slider::slide_dbl(
    z$value,
    .f = ~ stats::quantile(.x, probs = q, na.rm = TRUE, names = FALSE),
    .before = W - 1, .complete = FALSE
  )
  S_nat <- as.integer(z$value >= q_roll)
  r <- rle(S_nat); vals <- r$values; lens <- r$lengths
  for (i in seq_along(vals)) if (vals[i] == 1L && lens[i] < min_run) vals[i] <- 0L
  S_nat2 <- inverse.rle(list(values = vals, lengths = lens))
  
  dfp <- tibble(date = z$date, value = z$value, q_roll = q_roll, S = S_nat2)
  rects <- {
    r2 <- rle(dfp$S); ends <- cumsum(r2$lengths); starts <- c(1, head(ends, -1)+1)
    sel <- which(r2$values == 1L)
    if (!length(sel)) tibble(xmin=as.Date(character()), xmax=as.Date(character())) else
      tibble(xmin = dfp$date[starts[sel]], xmax = dfp$date[ends[sel]])
  }
  ggplot(dfp, aes(date, value)) +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey80", alpha = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_line(aes(y = q_roll), linewidth = 0.6, linetype = 2) +
    labs(title = "Rolling-quantile stress diagnostic",
         subtitle = paste0("q=", q, ", window=", window_years, "y, min_run=", min_run,
                           "  (shaded = S=1)"),
         x = NULL, y = "Index / rolling q-quantile") +
    theme_minimal(base_size = 12)
}

# Forward window 
fw_indicator <- function(event01, k, x){
  Tn <- length(event01); y <- integer(Tn)
  for (t in seq_len(Tn)) {
    lo <- t + k - x; hi <- t + k + x
    if (lo < 1L || hi > Tn) { y[t] <- NA_integer_; next }
    y[t] <- as.integer(any(event01[lo:hi] == 1L))
  }
  y
}
default_nw_lag <- function(Tn) floor(4 * (Tn/100)^(2/9))
var_alpha_plus_beta_hac <- function(fit, lag = NULL){
  V <- if (is.null(lag)) sandwich::NeweyWest(fit, prewhite = FALSE, adjust = TRUE)
  else                 sandwich::NeweyWest(fit, lag = lag, prewhite = FALSE, adjust = TRUE)
  e <- rep(1, length(coef(fit)))
  as.numeric(t(e) %*% V %*% e)
}

# single-source conditional curve
rs_curve_one <- function(df_event, df_cond, k_seq = 12:90, x = 5,
                         min_gap = 12L, min_keep = 30L, nw_lag = NULL) {
  a <- df_event %>% arrange(.data$date) %>% mutate(date = as.Date(.data$date))
  b <- df_cond  %>% arrange(.data$date) %>% mutate(date = as.Date(.data$date))
  inter <- intersect(a$date, b$date)
  if (length(inter) < 2L) return(tibble(k=integer(), x=integer(), p_hat=double(), se=double(), n=integer()))
  rng <- range(inter)
  a <- dplyr::filter(a, .data$date >= rng[1] & .data$date <= rng[2])
  b <- dplyr::filter(b, .data$date >= rng[1] & .data$date <= rng[2])
  stopifnot(nrow(a) == nrow(b), all(a$date == b$date))
  
  R <- as.integer(a$event); S <- as.integer(b$cond); Tn <- length(R)
  if (is.null(nw_lag)) nw_lag <- default_nw_lag(Tn)
  
  purrr::map_dfr(k_seq, function(k){
    if (k < min_gap) return(tibble(k=k, x=x, p_hat=NA_real_, se=NA_real_, n=NA_integer_))
    y <- fw_indicator(R, k, x); keep <- is.finite(y)
    n_keep <- sum(keep)
    if (n_keep < min_keep) return(tibble(k=k, x=x, p_hat=NA_real_, se=NA_real_, n=n_keep))
    d <- tibble(y = y[keep], S = S[keep])
    if (min(sum(d$S==1L), sum(d$S==0L)) < 5) return(tibble(k=k, x=x, p_hat=NA_real_, se=NA_real_, n=n_keep))
    fit   <- lm(y ~ S, data = d)
    p_hat <- sum(coef(fit), na.rm = TRUE)
    se    <- sqrt(pmax(var_alpha_plus_beta_hac(fit, nw_lag), 0))
    tibble(k=k, x=x, p_hat=p_hat, se=se, n=n_keep)
  })
}

# Joint uncertainty that respects R–S dependence
compute_p_vec_RS <- function(R, S, K, x){
  sapply(K, function(k){
    y <- fw_indicator(R, k, x)
    keep <- is.finite(y) & (S == 1L)
    if (!any(keep)) NA_real_ else mean(y[keep])
  })
}

# (A) bivariate MBB
mbb_indices <- function(Tn, ell) {
  starts <- sample.int(Tn, size = ceiling(Tn/ell), replace = TRUE)
  idx <- unlist(lapply(starts, function(s) ((s-1) + 0:(ell-1)) %% Tn + 1L))
  idx[seq_len(Tn)]
}
rs_joint_tests_bivar_mbb <- function(df_RS, K, x = 5, B = 1000, ell = NULL, seed = 123) {
  z <- df_RS %>% arrange(.data$date)
  R <- as.integer(z$event); S <- as.integer(z$crisis); Tn <- length(R)
  if (is.null(ell)) ell <- max(8L, round(2 * Tn^(1/3)))
  
  p_obs <- compute_p_vec_RS(R, S, K, x)
  if (any(!is.finite(p_obs))) stop("Not enough data for chosen K.")
  
  set.seed(seed)
  Pmat <- matrix(NA_real_, nrow = B, ncol = length(K))
  for (b in 1:B) {
    idx <- mbb_indices(Tn, min(ell, Tn))
    Rb <- R[idx]; Sb <- S[idx]
    Pmat[b, ] <- compute_p_vec_RS(Rb, Sb, K, x)
  }
  
  Sigma <- stats::cov(Pmat, use = "pairwise.complete.obs")
  pbar  <- mean(p_obs); e <- rep(1, length(K))
  ev <- eigen(Sigma, symmetric = TRUE); lam <- pmax(ev$values, .Machine$double.eps)
  Sinvhalf <- ev$vectors %*% diag(1/sqrt(lam)) %*% t(ev$vectors)
  eta <- as.numeric(Sinvhalf %*% (p_obs - pbar*e))
  stat <- sum(eta^2); pval <- 1 - pchisq(stat, df = length(K) - 1L)
  
  list(stat = stat, df = length(K) - 1L, pval = pval, Sigma_boot = Sigma, p_obs = p_obs)
}

# (B) peak-aligned: cut at union of R and S starts
starts01 <- function(x){
  up <- c(0L, diff(x)); idx <- which(up == 1L)
  if (x[1] == 1L) idx <- c(1L, idx)
  sort(unique(idx))
}
joint_segments_RS <- function(R, S){
  Tn <- length(R)
  cuts <- sort(unique(c(starts01(R), starts01(S))))
  if (!length(cuts)) return(list(list(R = R, S = S)))
  segs <- vector("list", length(cuts))
  for (i in seq_along(cuts)) {
    a <- cuts[i]
    b <- if (i < length(cuts)) cuts[i+1] - 1L else Tn
    segs[[i]] <- list(R = R[a:b], S = S[a:b])
  }
  segs
}
resample_joint_segments <- function(segs, Tn){
  outR <- integer(0); outS <- integer(0)
  while (length(outR) < Tn) {
    s <- segs[[sample.int(length(segs), 1)]]
    outR <- c(outR, s$R); outS <- c(outS, s$S)
  }
  list(R = outR[seq_len(Tn)], S = outS[seq_len(Tn)])
}
rs_joint_tests_bivar_peak <- function(df_RS, K, x = 5, B = 1000, seed = 123){
  z <- df_RS %>% arrange(.data$date)
  R <- as.integer(z$event); S <- as.integer(z$crisis); Tn <- length(R)
  
  p_obs <- compute_p_vec_RS(R, S, K, x)
  if (any(!is.finite(p_obs))) stop("Not enough data for chosen K.")
  
  segs <- joint_segments_RS(R, S)
  set.seed(seed)
  Pmat <- matrix(NA_real_, nrow = B, ncol = length(K))
  for (b in 1:B) {
    rs <- resample_joint_segments(segs, Tn)
    Pmat[b, ] <- compute_p_vec_RS(rs$R, rs$S, K, x)
  }
  
  Sigma <- stats::cov(Pmat, use = "pairwise.complete.obs")
  pbar  <- mean(p_obs); e <- rep(1, length(K))
  ev <- eigen(Sigma, symmetric = TRUE); lam <- pmax(ev$values, .Machine$double.eps)
  Sinvhalf <- ev$vectors %*% diag(1/sqrt(lam)) %*% t(ev$vectors)
  eta <- as.numeric(Sinvhalf %*% (p_obs - pbar*e))
  stat <- sum(eta^2); pval <- 1 - pchisq(stat, df = length(K) - 1L)
  
  list(stat = stat, df = length(K) - 1L, pval = pval, Sigma_boot = Sigma, p_obs = p_obs)
}

# MASTER RUNNER
quarter_end <- function(d) lubridate::ceiling_date(d, "quarter") - lubridate::days(1)

run_r_given_stress_full <- function(
    df_nber,
    stress_files = NULL,     # tibble: file,label (if qS_all not supplied)
    qS_all = NULL,           # {source,date,stress_high} OR {label,date,stress_high} OR {source,date,cond}
    q = 0.9, window_years = 12, min_run = 2L, quarter_rule = c("any","majority"),
    k_seq = 12:90, x_values = c(3,4,5),
    K_case1 = c(38,56), K_case2 = c(36:39,55:58),
    B = 2000, ci_lev = 0.80, block_len = NULL,
    build_diag_plots = FALSE
){
  # freeze args to plain scalars (no NSE lookups later) 
  quarter_rule <- match.arg(quarter_rule)
  q0 <- q; win0 <- window_years; minrun0 <- min_run; rule0 <- quarter_rule
  kseq0 <- k_seq; xvals0 <- x_values; K1 <- K_case1; K2 <- K_case2
  B0 <- B; ci0 <- ci_lev; block0 <- block_len
  
  stopifnot(all(c("date","recession") %in% names(df_nber)))
  event_q <- df_nber %>% dplyr::arrange(.data$date) %>%
    dplyr::transmute(date = as.Date(.data$date), event = as.integer(.data$recession))
  
  if (is.null(qS_all)) {
    stopifnot(exists("financial_q"),
              all(c("date","value","source") %in% names(financial_q)))
    q_grid <- as.Date(df_nber$date)
    
    qS_all <- financial_q %>%
      arrange(source, date) %>%
      dplyr::group_split(source, .keep = TRUE) %>%
      purrr::map_dfr(function(dfq){
        lab <- unique(dfq$source)
        stress_to_quarterly_indicator_rolling(
          df_native    = dfq %>% dplyr::select(date, value),
          q            = q,            # or q0 inside “frozen args” runner
          window_years = window_years, # or win0 if froze args
          min_run      = min_run,      # or minrun0
          quarter_rule = quarter_rule, # or rule0
          q_grid       = q_grid
        ) %>% dplyr::mutate(source = lab)
      })
  }
  
  
  # normalize qS_all => {source,date,cond} 
  stopifnot("date" %in% names(qS_all))
  qS_all <- qS_all %>% dplyr::mutate(date = as.Date(.data$date))
  if (!"source" %in% names(qS_all)) {
    if ("label" %in% names(qS_all)) qS_all <- dplyr::rename(qS_all, source = .data$label)
    else                            qS_all <- dplyr::mutate(qS_all, source = "source_1")
  }
  if (!"cond" %in% names(qS_all)) {
    if ("stress_high" %in% names(qS_all)) qS_all <- dplyr::rename(qS_all, cond = .data$stress_high)
    else stop("qS_all must contain either 'cond' or 'stress_high'.")
  }
  qS_all <- qS_all %>%
    dplyr::transmute(source = as.character(.data$source),
                     date   = .data$date,
                     cond   = as.integer(.data$cond)) %>%
    dplyr::arrange(.data$source, .data$date)
  
  # optional diagnostics
  diag_plots <- list()
  if (build_diag_plots && !is.null(stress_files) && nrow(stress_files)) {
    diag_plots <- purrr::map(seq_len(nrow(stress_files)), function(i){
      f <- stress_files$file[i]; lab <- stress_files$label[i]
      df_native <- read_stress_local(f)
      stress_rolling_diag_plot(df_native, q = q0, window_years = win0, min_run = minrun0) +
        ggplot2::ggtitle(paste0(lab, " — Rolling threshold diagnostic"))
    })
    names(diag_plots) <- stress_files$label
  }
  
  # curves per source and x
  curves <- qS_all %>%
    dplyr::select(source, date, cond) %>%
    split(.$source) %>%
    purrr::imap_dfr(function(dfS, src){
      df_cond <- dplyr::select(dfS, date, cond)
      inter <- intersect(df_cond$date, event_q$date)
      if (length(inter) < 2L)
        return(tibble::tibble(k=integer(), x=integer(), p_hat=double(), se=double(), n=integer(), source=character())[0,])
      out <- purrr::map_dfr(xvals0, ~ rs_curve_one(event_q, df_cond, k_seq = kseq0, x = .x))
      if (!nrow(out))
        return(tibble::tibble(k=integer(), x=integer(), p_hat=double(), se=double(), n=integer(), source=character())[0,])
      dplyr::mutate(out, source = src)
    })
  
  stopifnot(all(c("k","x","p_hat","se","n","source") %in% names(curves)))
  x_max <- max(xvals0); zval <- stats::qnorm((1 + ci0)/2)
  
  bands_df <- curves %>%
    dplyr::filter(.data$x == x_max, is.finite(.data$p_hat), is.finite(.data$se)) %>%
    dplyr::mutate(lo = pmax(0, .data$p_hat - zval*.data$se),
                  hi = pmin(1, .data$p_hat + zval*.data$se))
  
  peaks <- first_peak_k(curves, group_vars = c("source","x"))
  
  # joint tests (bivariate MBB + peak-aligned) 
  jt_tbl <- purrr::map_dfr(unique(qS_all$source), function(src){
    set.seed(1000 + as.integer(factor(src)))
    df_S  <- qS_all %>% dplyr::filter(.data$source == src) %>% dplyr::transmute(date, crisis = cond)
    df_RS <- dplyr::inner_join(df_S, event_q, by = "date")
    df_RS <- df_RS[, c("date","crisis","event")]
    if (sum(df_RS$crisis, na.rm = TRUE) < 8)
      return(tibble::tibble(source = src,
                            case1_mbb_pval = NA_real_, case2_mbb_pval = NA_real_,
                            case1_peak_pval = NA_real_, case2_peak_pval = NA_real_))
    Tn <- nrow(df_RS); ell <- block0 %||% max(8L, round(2 * Tn^(1/3)))
    jt1_mbb  <- rs_joint_tests_bivar_mbb (df_RS, K = K1, x = x_max, B = B0, ell = ell)
    jt2_mbb  <- rs_joint_tests_bivar_mbb (df_RS, K = K2, x = x_max, B = B0, ell = ell)
    jt1_peak <- rs_joint_tests_bivar_peak(df_RS, K = K1, x = x_max, B = B0)
    jt2_peak <- rs_joint_tests_bivar_peak(df_RS, K = K2, x = x_max, B = B0)
    tibble::tibble(source = src,
                   case1_mbb_pval  = jt1_mbb$pval,
                   case2_mbb_pval  = jt2_mbb$pval,
                   case1_peak_pval = jt1_peak$pval,
                   case2_peak_pval = jt2_peak$pval)
  })
  
  # plots 
  humps   <- data.frame(xmin = c(36,55), xmax = c(40,58))
  sub_byx <- paste0("Lines: x ∈ {", paste(xvals0, collapse=", "),
                    "}; shaded: ", 100*ci0, "% HAC CI at x = ±", x_max)
  sub_qts <- paste0("Shaded = S_t = 1 (rolling q=", q0,
                    ", window=", win0, "y, min_run=", minrun0, ", rule=", rule0, ")")
  
  p_by_x <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = humps, inherit.aes = FALSE,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                       fill = "red", alpha = 0.10) +
    ggplot2::geom_ribbon(data = bands_df, ggplot2::aes(k, ymin = lo, ymax = hi, group = source),
                         fill = "grey85") +
    ggplot2::geom_line(data = curves, ggplot2::aes(k, p_hat, color = factor(x, levels = xvals0)),
                       linewidth = 0.9) +
    ggplot2::facet_wrap(~ source, ncol = 2) +
    ggplot2::scale_color_manual(values = setNames(c("#1b9e77","#d95f02","#7570b3")[seq_along(xvals0)], xvals0),
                                name = "Window ±x") +
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    ggplot2::labs(title = "R|S: forward conditional probability by stress source",
                  subtitle = sub_byx,
                  x = "k (quarters ahead)", y = "Conditional probability") +
    ggplot2::theme_minimal(base_size = 12)
  
  p_bigx <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = humps, inherit.aes = FALSE,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                       fill = "red", alpha = 0.10) +
    ggplot2::geom_ribbon(data = bands_df, ggplot2::aes(k, ymin = lo, ymax = hi, group = source),
                         fill = "grey85") +
    ggplot2::geom_line(data = bands_df, ggplot2::aes(k, p_hat, color = source),
                       linewidth = 0.9) +
    ggplot2::facet_wrap(~ source, ncol = 2) +
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    ggplot2::labs(title = paste0("R|S at x = ±", x_max, " (HAC CI shaded)"),
                  x = "k (quarters ahead)", y = "Conditional probability", color = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  
  # Time-series panels
  stopifnot(exists("financial_q"),
            all(c("date","value","source") %in% names(financial_q)))
  
  native_all <- financial_q %>%
    dplyr::arrange(source, date) %>%
    dplyr::mutate(freq = "quarterly")
  
  p_native <- if (nrow(native_all)) {
    ggplot2::ggplot(native_all, ggplot2::aes(date, value, color = source)) +
      ggplot2::geom_line(linewidth = 0.7, alpha = 0.9, show.legend = FALSE) +
      ggplot2::facet_wrap(~ source, scales = "free_y", ncol = 1) +
      ggplot2::labs(title = "Financial stress (quarterly)", x = NULL, y = "Index level") +
      ggplot2::theme_minimal(base_size = 12)
  } else NULL
  
  quarterly_levels <- native_all %>%
    dplyr::select(source, date, value) %>%
    dplyr::arrange(source, date)
  
  qS_all_plot <- qS_all %>% dplyr::transmute(source=.data$source, date=.data$date, stress_high=.data$cond)
  
  quarterly_merged <- quarterly_levels %>%
    dplyr::left_join(qS_all_plot, by = c("date","source")) %>%
    dplyr::mutate(stress_high = ifelse(is.na(.data$stress_high), 0L, .data$stress_high))
  
  episodes_from_indicator <- function(df_q) {
    df_q <- df_q %>% dplyr::arrange(.data$date)
    z <- df_q$stress_high
    if (!length(z) || all(z == 0L | is.na(z))) {
      return(tibble::tibble(source=df_q$source[1], start=as.Date(character()), end=as.Date(character()))[0,])
    }
    r <- rle(z); ends <- cumsum(r$lengths); starts <- c(1, head(ends, -1) + 1); sel <- which(r$values == 1L)
    tibble::tibble(source = df_q$source[1],
                   start  = df_q$date[starts[sel]],
                   end    = quarter_end(df_q$date[ends[sel]]))
  }
  
  shades <- if (nrow(quarterly_merged)) {
    quarterly_merged %>% dplyr::group_by(.data$source) %>%
      dplyr::group_modify(~episodes_from_indicator(.x)) %>%
      dplyr::ungroup()
  } else tibble::tibble(source=character(), start=as.Date(character()), end=as.Date(character()))[0,]
  
  sub_qts <- paste0("Shaded = S_t = 1 (rolling q=", q0,
                    ", window=", win0, "y, min_run=", minrun0, ", rule=", rule0, ")")
  
  p_quarterly <- if (nrow(quarterly_merged)) {
    ggplot2::ggplot() +
      ggplot2::geom_rect(data = shades, ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = -Inf, ymax = Inf),
                         inherit.aes = FALSE, fill = "grey0", alpha = 0.35) +
      ggplot2::geom_line(data = quarterly_merged, ggplot2::aes(.data$date, .data$value, color = .data$source),
                         linewidth = 0.9, show.legend = FALSE) +
      ggplot2::facet_wrap(~ source, scales = "free_y", ncol = 1) +
      ggplot2::labs(title = "Financial stress (quarterly mean) with high-stress quarters shaded",
                    subtitle = sub_qts, x = NULL, y = "Quarterly mean level") +
      ggplot2::theme_minimal(base_size = 12)
  } else NULL
  
  nber_shades <- df_nber %>%
    dplyr::mutate(start = .data$date, end = quarter_end(.data$date)) %>%
    dplyr::filter(.data$recession == 1L) %>%
    dplyr::transmute(start, end)
  
  p_quarterly_nber <- if (nrow(quarterly_merged)) {
    ggplot2::ggplot() +
      ggplot2::geom_rect(data = nber_shades, ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = -Inf, ymax = Inf),
                         inherit.aes = FALSE, fill = "grey40", alpha = 0.20) +
      ggplot2::geom_rect(data = shades, ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = -Inf, ymax = Inf),
                         inherit.aes = FALSE, fill = "grey0", alpha = 0.35) +
      ggplot2::geom_line(data = quarterly_merged, ggplot2::aes(.data$date, .data$value, color = .data$source),
                         linewidth = 0.9, show.legend = FALSE) +
      ggplot2::facet_wrap(~ source, scales = "free_y", ncol = 1) +
      ggplot2::labs(title = "Financial stress (quarterly mean) with high-stress (light) and NBER recessions (dark)",
                    x = NULL, y = "Quarterly mean level") +
      ggplot2::theme_minimal(base_size = 12)
  } else NULL
  
  list(
    curves = curves, bands = bands_df, peaks = peaks, tests = jt_tbl,
    plot_by_x = p_by_x, plot_bigx = p_bigx,
    plot_native = p_native, plot_quarterly = p_quarterly, plot_quarterly_nber = p_quarterly_nber,
    diag_plots = diag_plots, qS_all = qS_all
  )
}
# out$curves: source, k, x, p_hat, se, n
peaks_stress <- peak_rows(out$curves, group_vars = c("source","x"), k_cap = 56L) %>%
  dplyr::transmute(
    series      = source,
    x           = as.integer(x),
    peak_k      = as.integer(k),
    peak_p      = as.numeric(p_hat),
    n_at_peak   = as.integer(n),
    se_at_peak  = as.numeric(se)
  )
# crp_real$data: k, x, p_hat
peaks_nber <- peak_rows(crp_real$data, group_vars = "x", k_cap = 56L) %>%
  dplyr::transmute(
    series = "NBER (R|R)",
    x       = as.integer(x),
    peak_k  = as.integer(k),
    peak_p  = as.numeric(p_hat)
  )




#### CALL
out <- run_r_given_stress_full(
  df_nber,  stress_files = NULL,   # <- prevents re-reading/re-aggregating
  qS_all       = qS_all, # <- the {source,date,stress_high} 
  q = 0.9, window_years = 12, min_run = 2L, quarter_rule = "majority",
  k_seq = 12:90, x_values = c(3,4,5),
  K_case1 = c(38,56), K_case2 = c(36:39,55:58),
  B = 1000, block_len = 24, ci_lev = 0.80,
  build_diag_plots = TRUE
)

print(out$plot_bigx); print(out$plot_by_x); print(out$tests)





make_wide_table <- function(df_stress, df_nber,
                            k_min = NULL, k_max = NULL,
                            p_min = NULL, p_max = NULL,
                            keep_out_of_bounds = FALSE,
                            flag = FALSE) {
  
  out_long <- dplyr::bind_rows(
    df_stress %>% dplyr::select(series, x, peak_k, peak_p),
    df_nber   %>% dplyr::select(series, x, peak_k, peak_p)
  ) %>%
    dplyr::arrange(series, x)
  
  # If empty, return a minimal empty table with expected columns
  if (nrow(out_long) == 0L) {
    x_levels <- paste0("\u00B1", sort(unique(c(df_stress$x, df_nber$x))))
    out <- tibble::tibble(series = character())
    if (length(x_levels)) {
      for (nm in x_levels) out[[nm]] <- character()
    }
    return(out)
  }
  
  # Build mask safely; use which(ok) for subsetting
  ok <- rep(TRUE, nrow(out_long))
  if (!is.null(k_min)) ok <- ok & (out_long$peak_k >= k_min)
  if (!is.null(k_max)) ok <- ok & (out_long$peak_k <= k_max)
  if (!is.null(p_min)) ok <- ok & (out_long$peak_p >= p_min)
  if (!is.null(p_max)) ok <- ok & (out_long$peak_p <= p_max)
  
  # Compose cell text now (used in both branches)
  base_cell <- sprintf("k = %d | p = %.3f", out_long$peak_k, out_long$peak_p)
  
  if (!keep_out_of_bounds) {
    out_long <- out_long[which(ok), , drop = FALSE]
    out_long$cell <- base_cell[which(ok)]
  } else {
    out_long$cell <- base_cell
    out_long$cell[!ok] <- if (flag) paste0(out_long$cell[!ok], " \u2020") else NA_character_
  }
  
  # If everything got filtered out, still return an empty wide with columns
  if (nrow(out_long) == 0L) {
    x_levels <- paste0("\u00B1", sort(unique(c(df_stress$x, df_nber$x))))
    out <- tibble::tibble(series = character())
    if (length(x_levels)) for (nm in x_levels) out[[nm]] <- character()
    return(out)
  }
  
  out_long <- out_long %>% dplyr::mutate(x_band = paste0("\u00B1", x))
  
  # Column order, even if some x are missing post-filter
  x_levels <- paste0("\u00B1", sort(unique(c(df_stress$x, df_nber$x))))
  
  out_long %>%
    dplyr::select(series, x_band, cell) %>%
    tidyr::pivot_wider(names_from = x_band, values_from = cell) %>%
    dplyr::mutate(series = as.character(series)) %>%
    dplyr::relocate(series) %>%
    # keep only expected x columns (those that existed pre-filter)
    { dplyr::select(., series, dplyr::all_of(intersect(x_levels, names(.)))) }
}
peaks_wide <- make_wide_table(
  df_stress = peaks_stress,
  df_nber   = peaks_nber,
  k_min = 7, k_max = 56,
  keep_out_of_bounds = TRUE,
  flag = TRUE
)

# Pretty print
knitr::kable(peaks_wide, align = c("l", rep("c", ncol(peaks_wide)-1)),
             caption = "Peak quarter (k) and probability by x, NBER vs. stress indexes")











###############################################################################
## STANDARD S|S  -  P( S future window | S_t = 1 )
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(ggplot2); library(sandwich)
})

# Helpers
make_window_indicator <- function(r, k, x){
  Tn <- length(r); y <- integer(Tn)
  for (t in seq_len(Tn)) {
    lo <- t + k - x; hi <- t + k + x
    if (lo < 1L || hi > Tn) { y[t] <- NA_integer_; next }
    y[t] <- as.integer(any(r[lo:hi] == 1L))
  }
  y
}

default_nw_lag <- function(Tn) floor(4 * (Tn/100)^(2/9))

# HAC var(α+β) but never prewhite; return NA cleanly if anything goes wrong
nw_var_alpha_plus_beta_safe <- function(fit, lag = NULL){
  p <- length(coef(fit))
  V <- tryCatch(
    {
      if (is.null(lag)) {
        suppressWarnings(sandwich::NeweyWest(fit, prewhite = FALSE, adjust = TRUE))
      } else {
        suppressWarnings(sandwich::NeweyWest(fit, lag = lag, prewhite = FALSE, adjust = TRUE))
      }
    },
    error = function(e) matrix(NA_real_, p, p)
  )
  if (!is.matrix(V) || any(!is.finite(V))) return(NA_real_)
  e <- rep(1, p)
  as.numeric(t(e) %*% V %*% e)
}

# Build static S_t per series
# S_t = 1{ value_t >= q-quantile over the full (quarterly) sample }, then
# collapse spikes shorter than min_run to 0.
build_static_S <- function(df_q, q = 0.95, min_run = 2L){
  z <- df_q %>% arrange(date) %>% transmute(date = as.Date(date), value = as.numeric(value))
  thr <- stats::quantile(z$value, probs = q, na.rm = TRUE, names = FALSE)
  S0  <- as.integer(z$value >= thr)
  r   <- rle(S0); vals <- r$values; lens <- r$lengths
  for (i in seq_along(vals)) if (vals[i] == 1L && lens[i] < min_run) vals[i] <- 0L
  S   <- inverse.rle(list(values = vals, lengths = lens))
  tibble(date = z$date, S = as.integer(S))
}

# Optionally build S for all sources in `financial_q`
build_static_S_all <- function(financial_q, q = 0.95, min_run = 2L){
  financial_q %>%
    dplyr::arrange(source, date) %>%
    group_by(source) %>%
    group_modify(~ build_static_S(.x, q = q, min_run = min_run)) %>%
    ungroup()
}

# One-series S|S curve 
# Event = S (same series), Condition = S (same series)
# p_hat(k,x) = α+β from OLS(y^{k,x} ~ 1 + S_t), HAC se for α+β
ss_curve_one <- function(df_S, k_seq = 12:90, x = 5,
                         min_gap = 12L, min_keep = 30L, nw_lag = NULL){
  z <- df_S %>% arrange(date) %>% transmute(date = as.Date(date), S = as.integer(S))
  R <- z$S                       # "event" series = S
  C <- z$S                       # "condition"   = S (same)
  Tn <- length(R)
  if (is.null(nw_lag)) nw_lag <- default_nw_lag(Tn)
  
  map_dfr(k_seq, function(k){
    if (k < min_gap) return(tibble(k=k, x=x, p_hat=NA_real_, se=NA_real_, n=NA_integer_))
    y <- make_window_indicator(R, k, x)
    keep <- is.finite(y)
    n_keep <- sum(keep)
    if (n_keep < min_keep) return(tibble(k=k, x=x, p_hat=NA_real_, se=NA_real_, n=n_keep))
    d <- tibble(y = y[keep], S = C[keep])
    # need both classes present for OLS to be identified
    if (min(sum(d$S==1L), sum(d$S==0L)) < 5) return(tibble(k=k, x=x, p_hat=NA_real_, se=NA_real_, n=n_keep))
    fit   <- lm(y ~ S, data = d)
    p_hat <- sum(coef(fit), na.rm = TRUE)                    # α + β
    se    <- sqrt(pmax(nw_var_alpha_plus_beta_safe(fit, nw_lag), 0))
    tibble(k=k, x=x, p_hat=p_hat, se=se, n=n_keep)
  })
}

# Joint tests across K
# (A) HAC (multivariate NW) on Z_t(k) = S_t * y^{k,x} / E[S_t]
joint_test_hac_S <- function(df_S, K, x = 5, L = NULL){
  z <- df_S %>% arrange(date) %>% transmute(S = as.integer(S))
  S <- z$S
  Y <- sapply(K, function(k) make_window_indicator(S, k, x))   # T x m
  keep <- apply(is.finite(Y), 1, all)
  Yk <- Y[keep, , drop=FALSE]; Sk <- S[keep]
  denom <- mean(Sk); if (!is.finite(denom) || denom <= 0) stop("E[S]=0; cannot condition on S=1.")
  Z  <- (Sk * Yk) / denom
  
  p_obs <- colMeans(Z)                     # use common-sample means
  Zc <- scale(Z, center = TRUE, scale = FALSE); T0 <- nrow(Zc)
  if (is.null(L)) L <- default_nw_lag(T0)  # <<< T0, not raw Tn
  
  Gamma0 <- crossprod(Zc)/T0
  w <- 1 - (1:L)/(L+1); Gsum <- matrix(0, ncol(Zc), ncol(Zc))
  for (l in 1:L){
    A <- t(Zc[(l+1):T0, , drop=FALSE]) %*% Zc[1:(T0-l), , drop=FALSE] / T0
    Gsum <- Gsum + w[l] * (A + t(A))
  }
  Sigma_hat <- (Gamma0 + Gsum) / T0       # <<< divide by T0
  
  pbar <- mean(p_obs); e <- rep(1, length(K))
  ev <- eigen(Sigma_hat, symmetric=TRUE); lam <- pmax(ev$values, .Machine$double.eps)
  Sminushalf <- ev$vectors %*% diag(1/sqrt(lam)) %*% t(ev$vectors)
  eta <- as.numeric(Sminushalf %*% (p_obs - pbar*e))
  stat <- sum(eta^2); pval <- 1 - pchisq(stat, df = length(K)-1L)
  list(stat = stat, df = length(K)-1L, pval = pval, Sigma_hat = Sigma_hat, p_obs = p_obs)
}

# (B) Episode bootstrap on S (peak-to-peak analogue for S=1 starts)
find_S_starts <- function(S){
  up <- c(0L, diff(S)); idx <- which(up == 1L)
  if (length(idx) == 0L && S[1] == 1L) idx <- 1L
  idx
}
resample_S_episodes <- function(S, T_target){
  idx <- find_S_starts(S); if (!length(idx)) stop("No S=1 episodes.")
  segs <- lapply(seq_along(idx), function(i){
    a <- idx[i]; b <- if (i < length(idx)) idx[i+1]-1L else length(S)
    S[a:b]
  })
  out <- integer(0)
  while (length(out) < T_target) out <- c(out, segs[[sample.int(length(segs),1)]])
  out[seq_len(T_target)]
}
joint_test_bootstrap_S <- function(df_S, K, x = 5, B = 1000, seed = 123){
  z <- df_S %>% arrange(date) %>% transmute(S = as.integer(S))
  S <- z$S; Tn <- length(S)
  p_obs <- sapply(K, function(k){
    y <- make_window_indicator(S, k, x); keep <- is.finite(y) & (S == 1L)
    if (!any(keep)) NA_real_ else mean(y[keep])
  })
  if (any(!is.finite(p_obs))) stop("Not enough data for chosen K.")
  set.seed(seed)
  Pmat <- matrix(NA_real_, nrow = B, ncol = length(K))
  for (b in 1:B){
    Sb <- resample_S_episodes(S, Tn)
    Pmat[b, ] <- sapply(K, function(k){
      y <- make_window_indicator(S, k, x); keep <- is.finite(y) & (Sb == 1L)
      if (!any(keep)) NA_real_ else mean(y[keep])
    })
  }
  Sigma <- stats::cov(Pmat, use = "pairwise.complete.obs")
  pbar <- mean(p_obs); e <- rep(1, length(K))
  ev <- eigen(Sigma, symmetric=TRUE); lam <- pmax(ev$values, .Machine$double.eps)
  Sinvhalf <- ev$vectors %*% diag(1/sqrt(lam)) %*% t(ev$vectors)
  eta <- as.numeric(Sinvhalf %*% (p_obs - pbar*e))
  stat <- sum(eta^2); pval <- 1 - pchisq(stat, df = length(K)-1L)
  list(stat = stat, df = length(K)-1L, pval = pval, Sigma_boot = Sigma, p_obs = p_obs)
}

# Main runner 
run_ss_standard <- function(
    financial_q,
    sources     = NULL,          # default: all in financial_q
    q           = 0.95,
    min_run     = 2L,
    k_seq       = 12:90,
    x_values    = c(3,4,5),
    ci_lev      = 0.80,
    K_case1     = c(38,56),
    K_case2     = c(36:39,55:58),
    B           = 2000
){
  stopifnot(all(c("date","value","source") %in% names(financial_q)))
  fq <- financial_q %>% arrange(source, date) %>% mutate(date = as.Date(date))
  if (is.null(sources)) sources <- unique(fq$source)
  
  # Build S_t per source (static quantile)
  qS_all <- fq %>%
    dplyr::filter(source %in% sources) %>%
    group_by(source) %>%
    group_modify(~ build_static_S(.x, q = q, min_run = min_run)) %>%
    ungroup()
  
  # Curves for each source and x
  curves <- qS_all %>%
    split(.$source) %>%
    imap_dfr(function(dfS, src){
      out <- map_dfr(x_values, ~ ss_curve_one(dfS %>% dplyr::select(date, S),
                                              k_seq = k_seq, x = .x))
      if (!nrow(out)) return(tibble(k=integer(), x=integer(), p_hat=double(), se=double(), n=integer(), source=character())[0,])
      mutate(out, source = src)
    })
  
  # Bands at largest x
  x_max <- max(x_values); z <- qnorm((1 + ci_lev)/2)
  bands <- curves %>%
    dplyr::filter(x == x_max, is.finite(p_hat), is.finite(se)) %>%
    dplyr::mutate(lo = pmax(0, p_hat - z*se), hi = pmin(1, p_hat + z*se))
  
  peaks <- first_peak_k(curves, group_vars = c("source","x"))
  
  # Joint tests per source, at x = x_max
  tests <- map_dfr(unique(qS_all$source), function(src){
    df_S <- qS_all %>% dplyr::filter(source == src) %>% dplyr::select(date, S)
    if (sum(df_S$S, na.rm = TRUE) < 8)
      return(tibble(source = src,
                    case1_boot_pval = NA_real_, case1_hac_pval = NA_real_,
                    case2_boot_pval = NA_real_, case2_hac_pval = NA_real_))
    jt1b <- joint_test_bootstrap_S(df_S, K = K_case1, x = x_max, B = B)
    jt2b <- joint_test_bootstrap_S(df_S, K = K_case2, x = x_max, B = B)
    jt1h <- joint_test_hac_S     (df_S, K = K_case1, x = x_max)
    jt2h <- joint_test_hac_S     (df_S, K = K_case2, x = x_max)
    tibble(source = src,
           case1_boot_pval = jt1b$pval, case1_hac_pval = jt1h$pval,
           case2_boot_pval = jt2b$pval, case2_hac_pval = jt2h$pval)
  })
  
  # Plots
  humps <- data.frame(xmin = c(36,55), xmax = c(40,58))
  p_by_x <- ggplot() +
    geom_rect(data = humps, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "red", alpha = 0.10) +
    geom_ribbon(data = bands, aes(k, ymin = lo, ymax = hi, group = source),
                fill = "grey85") +
    geom_line(data = curves, aes(k, p_hat, color = factor(x, levels = sort(unique(x)))),
              linewidth = 0.9) +
    facet_wrap(~ source, ncol = 2) +
    scale_color_manual(
      values = setNames(c("#1b9e77","#d95f02","#7570b3")[seq_along(sort(unique(curves$x)))],
                        sort(unique(curves$x))),
      name = "Window ±x"
    ) +
    coord_cartesian(ylim = c(0,1)) +
    labs(title = "S|S: forward conditional probability by stress source",
         subtitle = paste0("Lines: x ∈ {", paste(sort(unique(curves$x)), collapse=", "),
                           "}; shaded: ", 100*ci_lev, "% HAC CI at largest x"),
         x = "k (quarters ahead)", y = "Conditional probability") +
    theme_minimal(base_size = 12)
  
  p_bigx <- ggplot() +
    geom_rect(data = humps, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "red", alpha = 0.10) +
    geom_ribbon(data = bands, aes(k, ymin = lo, ymax = hi, group = source),
                fill = "grey85") +
    geom_line(data = bands, aes(k, (lo + hi)/2, color = source),
              linewidth = 0.9) +
    facet_wrap(~ source, ncol = 2) +
    coord_cartesian(ylim = c(0,1)) +
    labs(title = "S|S at largest x (HAC CI shaded)",
         x = "k (quarters ahead)", y = "Conditional probability", color = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  
  list(curves = curves, bands = bands, peaks = peaks, tests = tests,
       plot_by_x = p_by_x, plot_bigx = p_bigx,
       meta = list(q = q, min_run = min_run, k_seq = k_seq, x_values = x_values,
                   ci_lev = ci_lev, K_case1 = K_case1, K_case2 = K_case2))
}



#### CALL
sources_use <- unique(financial_q$source)

out_ss <- run_ss_standard(
  financial_q,
  sources  = sources_use,
  q        = 0.9,          #  0.8/0/9 fixed sample quantile (standard, NOT rolling)
  min_run  = 2L,
  k_seq    = 12:90,
  x_values = c(3,4,5),
  ci_lev   = 0.80,
  K_case1  = c(38,56),
  K_case2  = c(36:39,55:58),
  B        = 1000
)

print(out_ss$plot_by_x)
print(out_ss$plot_bigx)
out_ss$peaks
out_ss$tests






































################################################################################
## LOCAL SURFACE (heatmap)
################################################################################

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(ggplot2) })

# Helper: assert the grouping column exists and get its symbol
.ensure_group <- function(df, group_var){
  stopifnot(is.character(group_var), length(group_var)==1L, group_var %in% names(df))
  rlang::sym(group_var)
}

# 1) Surface builder (robust keys + robust grid column creation)
make_local_surface <- function(curves,
                               group_var = "window_label",
                               x_sel = NULL,
                               k_min = 12L, k_max = 90L,
                               drop_all_na_groups = TRUE) {
  gsym <- .ensure_group(curves, group_var)
  
  z <- curves %>%
    dplyr::filter(is.finite(.data$p_hat)) %>%
    dplyr::mutate(!!gsym := as.character(.data[[group_var]]))
  
  if (is.null(x_sel)) x_sel <- max(z$x, na.rm = TRUE)
  
  # preserve factor order if present, otherwise data order
  lvl <- if (is.factor(curves[[group_var]])) {
    levels(curves[[group_var]])
  } else {
    z %>% dplyr::arrange(!!gsym) %>% dplyr::distinct(!!gsym) %>% dplyr::pull(!!gsym)
  }
  
  # grid with correct column NAME (no quasiquotation pitfalls)
  grid <- tidyr::crossing(
    !!!setNames(list(factor(lvl, levels = lvl)), group_var),
    k = as.integer(seq(k_min, k_max))
  )
  
  z_sub <- z %>%
    dplyr::filter(.data$x == x_sel) %>%
    dplyr::select(!!gsym, k, p_hat, se, n)
  
  # >>> KEY FIX: pass an UNNAMED by with the *value* of group_var
  surf <- dplyr::left_join(grid, z_sub, by = c(group_var, "k")) %>%
    dplyr::mutate(x_sel = as.integer(x_sel))
  
  if (drop_all_na_groups) {
    keep <- surf %>%
      dplyr::group_by(!!gsym) %>%
      dplyr::summarise(all_na = all(!is.finite(p_hat)), .groups = "drop") %>%
      dplyr::filter(!all_na) %>% dplyr::pull(!!gsym)
    surf <- surf %>%
      dplyr::filter(.data[[group_var]] %in% keep) %>%
      dplyr::mutate(!!gsym := droplevels(factor(.data[[group_var]], levels = lvl)))
  }
  surf
}

# 2) Ridge extractor (first peak in [7,56]) — unchanged logic, safer factoring
make_ridge_from_curves <- function(curves,
                                   group_var = "window_label",
                                   x_sel = NULL,
                                   k_min_peak = 7L, k_max_peak = 56L){
  gsym <- .ensure_group(curves, group_var)
  z <- curves %>%
    dplyr::filter(is.finite(p_hat)) %>%
    dplyr::mutate(!!gsym := as.character(.data[[group_var]]))
  
  if (is.null(x_sel)) x_sel <- max(z$x, na.rm = TRUE)
  
  rid <- z %>%
    dplyr::filter(.data$x == x_sel) %>%
    first_peak_k(group_vars = group_var,
                 k_col  = "k", p_col = "p_hat",
                 k_min  = k_min_peak, k_max = k_max_peak) %>%
    dplyr::rename(peak_k = first_peak_k, peak_p = first_peak_prob)
  
  rid[[group_var]] <- factor(rid[[group_var]],
                             levels = z %>% dplyr::arrange(!!gsym) %>% dplyr::pull(!!gsym) %>% unique())
  rid$x_sel <- as.integer(x_sel)
  rid
}

# 3) Plotters (unchanged API). Minor internal robustness on aes
plot_local_surface <- function(surf, ridge = NULL,
                               group_var = "window_label",
                               title = "Local surface of conditional probability",
                               subtitle = NULL,
                               p_limits = c(0,1),
                               se_cap = NULL,
                               n_min  = NULL) {
  gsym <- .ensure_group(surf, group_var)
  df <- surf
  
  df$alpha <- 1
  if (!is.null(se_cap) && "se" %in% names(df))
    df$alpha <- pmin(df$alpha, ifelse(!is.na(df$se), pmax(0.2, 1 - (df$se/se_cap)), 0.2))
  if (!is.null(n_min) && "n" %in% names(df))
    df$alpha <- pmin(df$alpha, ifelse(!is.na(df$n), pmax(0.2, pmin(1, df$n / n_min)), 0.2))
  
  p <- ggplot(df, aes(x = k, y = !!gsym, fill = p_hat, alpha = alpha)) +
    geom_tile() +
    scale_fill_viridis_c(limits = p_limits, option = "C", na.value = "grey90") +
    scale_alpha_identity() +
    labs(title = title,
         subtitle = subtitle %||% paste0("x = ±", unique(df$x_sel)),
         x = "k (quarters since recession start)", y = group_var,
         fill = "p̂(k)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  if (!is.null(ridge) && nrow(ridge)) {
    p <- p + geom_point(data = ridge,
                        aes(x = peak_k, y = .data[[group_var]]),
                        inherit.aes = FALSE, shape = 21, size = 2.2,
                        fill = "white", color = "black", stroke = 0.4)
  }
  p
}

sizer_panel <- function(surf,
                        group_var = "window_label",
                        k_min = min(surf$k, na.rm=TRUE)+1L,
                        k_max = max(surf$k, na.rm=TRUE)-1L,
                        se_cap = NULL, n_min = NULL) {
  gsym <- .ensure_group(surf, group_var)
  df <- surf %>%
    dplyr::filter(k >= k_min, k <= k_max) %>%
    dplyr::arrange(!!gsym, k) %>%
    dplyr::group_by(!!gsym) %>%
    dplyr::mutate(dp = dplyr::lead(p_hat) - dplyr::lag(p_hat)) %>%
    dplyr::ungroup()
  
  df$sign <- ifelse(!is.finite(df$dp), NA_character_,
                    ifelse(df$dp > 0, "up", ifelse(df$dp < 0, "down", "flat")))
  df$alpha <- 1
  if (!is.null(se_cap) && "se" %in% names(df))
    df$alpha <- pmin(df$alpha, ifelse(!is.na(df$se), pmax(0.2, 1 - (df$se/se_cap)), 0.2))
  if (!is.null(n_min) && "n" %in% names(df))
    df$alpha <- pmin(df$alpha, ifelse(!is.na(df$n), pmax(0.2, pmin(1, df$n / n_min)), 0.2))
  
  ggplot(df, aes(x = k, y = !!gsym, alpha = alpha)) +
    geom_tile(aes(fill = sign)) +
    scale_fill_manual(values = c("up"="#1b9e77","down"="#d95f02","flat"="grey70"),
                      na.value = "grey90") +
    scale_alpha_identity() +
    labs(title = "SiZer-style slope map of p̂(k)",
         subtitle = paste0("green = rising; orange = falling; grey = flat; x = ±", unique(df$x_sel)),
         x = "k", y = group_var, fill = "sign(Δp̂)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}


################################################################################
##   (A) Rolling R|R windows: out_rr_roll$curves with 'window_label'
##   (B) Per-source R|S
################################################################################

## (A) Rolling surface over windows (x = largest by default)

surf_roll  <- make_local_surface(out_rr_roll$curves, group_var = "window_label",
                                 x_sel = NULL, k_min = 12, k_max = 90)
ridge_roll <- make_ridge_from_curves(out_rr_roll$curves, group_var = "window_label",
                                     x_sel = unique(surf_roll$x_sel)[1], k_min_peak = 7, k_max_peak = 56)
print(plot_local_surface(surf_roll, ridge_roll, group_var = "window_label",
                         title = "Rolling local surface: p̂(k | R_t=1) by window",
                         subtitle = "Heatmap across rolling windows; dot = first peak in [7,56]",
                         se_cap = 0.15, n_min = 40))
print(sizer_panel(surf_roll, group_var = "window_label", se_cap = 0.15, n_min = 40))



## (B) Per-source surface (R|S), x = ±5 
x_choice <- 5L
surf_src <- make_local_surface(out$curves, group_var = "source",
                               x_sel = x_choice, k_min = 12, k_max = 90)
ridge_src <- make_ridge_from_curves(out$curves, group_var = "source",
                                    x_sel = x_choice, k_min_peak = 7, k_max_peak = 56)

p_src_surface <- plot_local_surface(
  surf_src, ridge_src, group_var = "source",
  title = paste0("Local surface by stress source (R|S), x = ±", x_choice),
  subtitle = "Heatmap over k; dot = first peak in [7,56]",
  se_cap = 0.15, n_min = 40
)
print(p_src_surface)



 
 
 
 
 
 
 
 
 

