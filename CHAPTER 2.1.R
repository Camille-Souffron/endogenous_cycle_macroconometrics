################################################################################
# CHAPTER 2 - UNIVARIATE ANALYSIS (Section 1)
################################################################################

# Camille Souffron - MASTER THESIS APE (PSE)

#NB: all packages in the DATA_LOADING file



# ONLY IF NOT LOG HOURS. IF LOG HOURS: KEEP DATA_MERGED_REAL THE FROM DATA_LOADING FILE
data_merged_real <- data_merged_fin
data_merged_real <- data_merged_fin %>% rename(log_hours = finance)


###################################
# Choose frequency:    "q" = quarterly (Cutoff:1/80) ; delta= 0.05 ; N = 40
#                      "m" = monthly (Cutoff:1/240) ; delta = 0.0171 ; N = 120
#                      "w" = weekly (Cutoff:1/1040) ; delta = 0.00393 ; N = 520
# Cutoff > 20y
# Truncation N: last 10y
# delta: decay of 5% per quarter (1 - (1 - 0.05)^(4 / obs_per_year_real) => local half‑life ≈ 13.5q
# Horizon: 1000q ≈ 250yrs => ~28 cycles (9yrs each): 1000, 3000, 13000
# Horizon in power/CI plots :	1000,	3000,	13000
################################### 
# If LOG TS (standard deviation): = 100, otherwise, = 1
LOG = 1
# LOG = 1
###
freq_real <- "q"                    # "q", "m", or "w"

obs_per_year_real <- switch(freq_real,
                            q = 4,
                            m = 12,
                            w = 52)        
obs_unit_real     <- switch(freq_real,
                            q = "quarter",
                            m = "month",
                            w = "week")

##### Set the discount factor (memory half-life decay)
#ROBUSTNESS CHECK
# delta <- 1 - (1 - 0.001)^(4 / obs_per_year_real)
# delta <- 1 - (1 - 0.01)^(4 / obs_per_year_real)
delta <- 1 - (1 - 0.05)^(4 / obs_per_year_real)
# delta <- 1 - (1 - 0.1)^(4 / obs_per_year_real)
# delta <- 1 - (1 - 0.2)^(4 / obs_per_year_real)

N <- 10 * obs_per_year_real # to BUILD the first H_{t-1}


##########################################################################################
# PICK A DETRENDING METHOD 
#   "ideal20y" => high‑pass removing > 20 years
# ROBUSTNESS CHECKS:
#   "bandpass"=> Baxter-King band‑pass (6-32 Q)
#   "hp1600"  => Hodrick–Prescott (λ = 1600) (HIGH λ: full cycle, LOW λ: full trend)
#   "poly2"–"poly5" => 2nd–5th‑order polynomial trends
#   "linear"  => simple linear time trend
#   "log"     => log(t) trend
#   "none"    => no detrending
#   "ssa_full"   => SSA on levels (ex post), returns the oscillatory group as cycle
##########################################################################################

detrend_method <- "ideal20y"



to_quarterly_mean <- function(df, value_col) {
  stopifnot("date" %in% names(df), value_col %in% names(df))
  d <- df[, c("date", value_col)]
  if (!inherits(d$date, "Date")) d$date <- as.Date(d$date)
  
  # Quarter tag from the date (works for weekly or monthly dates)
  d$qtr <- as.yearqtr(d$date)
  
  # Mean within each quarter
  out <- aggregate(d[[value_col]], by = list(qtr = d$qtr),
                   FUN = function(x) mean(x, na.rm = TRUE))
  names(out)[2] <- value_col
  
  # Use end-of-quarter as the quarterly date
  out$date <- as.Date(out$qtr, frac = 1)
  out$qtr <- NULL
  
  # Order and return (date, value)
  out[order(out$date), c("date", value_col)]
}

# Example:
data_merged_real <- to_quarterly_mean(data_merged_real, "log_hours")
head(data_merged_real)


freq_real <- "q"                    

obs_per_year_real <- switch(freq_real,
                            q = 4,
                            m = 12,
                            w = 52)        
obs_unit_real     <- switch(freq_real,
                            q = "quarter",
                            m = "month",
                            w = "week")
forecast_starter <- 3*obs_per_year_real #3y = 1961:Q1
# how many prior obs needed before estimating
N_accumulate  <- 10 * obs_per_year_real           # 10‑year accumulation
N_ar          <- 2                           # AR(2) needs 2 lags
warmup        <- max(N_accumulate, N_ar)



#######-
# Plot raw series
#######
plot(
  data_merged_real$date,
  data_merged_real$log_hours,
  type = "l",
  xlab = "Date",
  ylab = "",
  main = ""
)






############################
# UNIT ROOT TESTS 
############################

run_unit_root_tests <- function(x, model = c("level","trend")) {
  # ensure we use the new signature
  model <- match.arg(model)
  x_num <- na.omit(as.numeric(x))
  if(length(x_num) < 20) stop("Series too short for unit‐root tests.")
  
  # 1) Augmented Dickey–Fuller via urca::ur.df
  #    "drift" = intercept only;  "trend" = intercept + time trend
  adf_type <- if(model=="level") "drift" else "trend"
  adf_res  <- ur.df(x_num, type=adf_type, selectlags="AIC")
  
  # 2) Phillips–Perron via urca::ur.pp
  pp_model <- if(model=="level") "constant" else "trend"
  pp_res   <- ur.pp(x_num, type="Z-alpha", model=pp_model, lags="short")
  
  # 3) KPSS via tseries::kpss.test
  #    null="Level" → I(0) around a constant;  "Trend" → I(0) around a trend
  kpss_null <- if(model=="level") "Level" else "Trend"
  kpss_res  <- kpss.test(x_num, null=kpss_null)
  
  # 4) DF–GLS via urca::ur.ers
  dfgls_res <- ur.ers(x_num,
                      type   = "DF-GLS",
                      model  = pp_model,    # match the PP spec
                      lag.max= trunc(4 * (length(x_num)/100)^(1/4))
  )
  
  # Return the summaries of each
  list(
    ADF   = summary(adf_res),
    PP    = summary(pp_res),
    KPSS  = kpss_res,
    DFGLS = summary(dfgls_res)
  )
}

# raw spread (unfiltered)
spread <- data_merged_real$log_hours
# test around a constant mean
res_level <- run_unit_root_tests(spread, model="level")
# test around a deterministic trend
res_trend <- run_unit_root_tests(spread, model="trend")

print(res_level)
print(res_trend)




########################################################
# High-Pass Filter Function
########################################################
# This function zeros out all Fourier components with absolute frequency lower than cutoff_freq.
ideal_highpass_filter <- function(x, cutoff_freq = 1 / (20 * obs_per_year_real)) {
  N <- length(x)
  # 1) Compute the DFT of x.
  X <- fft(x)
  
  # 2) Determine the frequency associated with each DFT component.
  # For k = 0,...,N-1, the frequency is k/N for k <= N/2 and (k-N)/N for k > N/2.
  k <- 0:(N-1)
  freq_k <- ifelse(k <= N/2, k/N, (k-N)/N)
  
  # 3) Create a mask that retains frequencies with absolute value >= cutoff_freq.
  keep <- abs(freq_k) >= cutoff_freq
  
  # 4) Zero out low-frequency coefficients.
  X_filtered <- X
  X_filtered[!keep] <- 0
  
  # 5) Inverse FFT to return the high-frequency time series.
  x_high <- Re(fft(X_filtered, inverse = TRUE) / N)
  return(x_high)
}


############################ SSA helpers (required by get_cycle: 'ssa_full' / 'ssa_full_rt') 
.require_Rssa <- function() {
  if (!requireNamespace("Rssa", quietly = TRUE))
    stop("Package 'Rssa' is required for SSA. Do: install.packages('Rssa')")
}

# Safe accessor for cfg with defaults 
.get_cfg <- function(name, default) {
  if (exists("cfg", inherits = TRUE) && !is.null(cfg[[name]])) cfg[[name]] else default
}

# Self-contained SSA (levels) decomposer 
.ssa_full_decompose_levels <- function(x,
                                       L        = 32L,
                                       max_K    = 20L,
                                       trend_min_period_obs = 60L,
                                       trend_lowshare_min   = 0.80,
                                       pair_wcor_min        = 0.60,
                                       pair_centroid_tol    = 0.10) {
  if (!requireNamespace("Rssa", quietly = TRUE))
    stop("Package 'Rssa' is required: install.packages('Rssa')")
  
  x <- as.numeric(x); n <- length(x)
  if (n < max(L + 8L, 50L)) stop("SSA needs longer series (n ≥ max(L+8, 50)).")
  
  S <- Rssa::ssa(x, L = min(L, n - 1L), kind = "1d-ssa")
  K <- min(max_K, length(S$sigma))
  
  # reconstruct first K singletons
  RC <- vector("list", K)
  for (k in seq_len(K)) RC[[k]] <- Rssa::reconstruct(S, groups = list(k = k))$k
  
  # weighted correlation (for pairing)
  W <- suppressWarnings(Rssa::wcor(S, groups = as.list(seq_len(K))))
  
  # local helpers
  .spec1 <- function(y) {
    sp <- spec.pgram(y, spans = NULL, pad = 0, taper = 0, detrend = FALSE, plot = FALSE)
    list(freq = as.numeric(sp$freq), spec = pmax(as.numeric(sp$spec), 0))
  }
  .stats <- function(y) {
    sp <- .spec1(y); f <- sp$freq; S <- sp$spec
    wS <- S / sum(S, na.rm = TRUE)
    centroid <- sum(f * wS, na.rm = TRUE)
    list(
      centroid = centroid,
      lowshare = function(f_cut) sum(wS[f <= f_cut], na.rm = TRUE)
    )
  }
  .is_pair <- function(i, j, pair_wcor_min, pair_centroid_tol, f_nyq = 0.5) {
    wcor_ij <- suppressWarnings(W[i, j])
    if (!is.finite(wcor_ij) || wcor_ij < pair_wcor_min) return(FALSE)
    d <- abs(stats_list[[i]]$centroid - stats_list[[j]]$centroid)
    d <= pair_centroid_tol * f_nyq
  }
  
  stats_list <- lapply(RC, .stats)
  f_cut <- 1 / trend_min_period_obs
  
  # trend-like singletons (heavy very-low-freq share)
  trend_idxs <- which(vapply(seq_len(K), function(k) {
    lf <- try(stats_list[[k]]$lowshare(f_cut), silent = TRUE)
    is.finite(lf) && lf >= trend_lowshare_min
  }, logical(1)))
  
  # oscillatory pairs among the rest (adjacent quasi-sin/cos with high wcor)
  used <- rep(FALSE, K); used[trend_idxs] <- TRUE
  pairs <- list(); k <- 1L
  while (k < K) {
    if (!used[k] && !used[k + 1L] &&
        .is_pair(k, k + 1L, pair_wcor_min, pair_centroid_tol)) {
      pairs[[length(pairs) + 1L]] <- c(k, k + 1L)
      used[k] <- used[k + 1L] <- TRUE
      k <- k + 2L; next
    }
    k <- k + 1L
  }
  
  grp_trend <- sort(unique(trend_idxs))
  grp_cycle <- sort(unique(unlist(pairs)))
  grp_noise <- setdiff(seq_len(K), union(grp_trend, grp_cycle))
  
  rec <- function(idxs) if (length(idxs)) Rssa::reconstruct(S, groups = list(g = idxs))$g else rep(0, n)
  
  list(
    trend  = as.numeric(rec(grp_trend)),
    cycle  = as.numeric(rec(grp_cycle)),
    noise  = as.numeric(rec(grp_noise)),
    groups = list(trend = grp_trend, cycle = grp_cycle, noise = grp_noise),
    S      = S
  )
}

# Expanding-sample (real-time) SSA wrapper (causal) 
.ssa_full_realtime_levels <- function(x,
                                      L        = 32L,
                                      max_K    = 20L,
                                      trend_min_period_obs = 60L,
                                      trend_lowshare_min   = 0.80,
                                      pair_wcor_min        = 0.60,
                                      pair_centroid_tol    = 0.10) {
  x <- as.numeric(x); n <- length(x)
  cyc_rt   <- rep(NA_real_, n)
  trend_rt <- rep(NA_real_, n)
  start_t <- max(L + 8L, 50L)
  if (n < start_t) return(list(trend = trend_rt, cycle = cyc_rt))
  for (t in seq.int(start_t, n)) {
    decomp <- .ssa_full_decompose_levels(
      x[1:t], L, max_K, trend_min_period_obs, trend_lowshare_min,
      pair_wcor_min, pair_centroid_tol
    )
    trend_rt[t] <- tail(decomp$trend, 1)
    cyc_rt[t]   <- tail(decomp$cycle, 1)
  }
  list(trend = trend_rt, cycle = cyc_rt)
}





################################################################################################################
# EXTRACTION OF THE CYCLICAL COMPONENT
################################################################################################################
# always return BOTH cycle & trend
get_cycle_trend <- function(x, method) {
  stopifnot(is.numeric(x), length(x) > 5L)
  
  # helper for methods that yield only a cycle estimate
  as_cycle_with_residual_trend <- function(cyc, x) {
    list(cycle = as.numeric(cyc), trend = as.numeric(x - cyc),
         meta  = list(source = "residual_trend"))
  }
  
  switch(method,
         
         ## Brick-Wall High-Pass based on years N 
         ideal20y  = { cutoff <- 1 / (20  * obs_per_year_real)
         cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         
         ideal100y = { cutoff <- 1 / (100 * obs_per_year_real)
         cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         
         ideal90y  = { cutoff <- 1 / ( 90 * obs_per_year_real)
         cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         
         ideal70y  = { cutoff <- 1 / ( 70 * obs_per_year_real)
         cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         
         ideal60y  = { cutoff <- 1 / ( 60 * obs_per_year_real)
         cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         
         ideal50y  = { cutoff <- 1 / ( 50 * obs_per_year_real)
         cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         
         ## Baxter-King Band-Pass (6-32 Q) 
         bandpass = { 
           if (!requireNamespace("mFilter", quietly = TRUE))
             stop("Package 'mFilter' is required: install.packages('mFilter')")
           
           # convert year bounds to periods given the sampling freq
           pl <- as.integer(round(yrs_low  * obs_per_year_real))   # lower period (e.g., 6 qtrs)
           pu <- as.integer(round(yrs_high * obs_per_year_real))   # upper period (e.g., 32 qtrs)
           
           # "Regular" BK for quarterly business cycles typically uses K = 12 (loses 12 obs each end)
           K  <- if (obs_per_year_real == 4) 12L else max(12L, pu %/% 3L)
           
           xts <- ts(as.numeric(x), frequency = obs_per_year_real)
           bk  <- mFilter::bkfilter(xts, pl = pl, pu = pu, nfix = K)
           
           # BK already returns both cycle and trend; cycle has NA at the ends (length K on each side)
           list(
             cycle = as.numeric(bk$cycle),
             trend = as.numeric(bk$trend),
             meta  = list(source = "bkfilter", pl = pl, pu = pu, K = K)
           )
         },
         ## Hodrick–Prescott (λ = 1600) 
         hp1600 = {
           out <- hpfilter(ts(x, freq = obs_per_year_real), freq = 1600, type = "lambda")
           # hpfilter already returns both
           list(cycle = as.numeric(out$cycle),
                trend = as.numeric(out$trend),
                meta  = list(source = "hpfilter"))
         },
         
         ##  Polynomial / linear / log trends (cycle = residuals) 
         poly2 = { t <- seq_along(x); cyc <- resid(lm(x ~ t + I(t^2)))
         as_cycle_with_residual_trend(cyc, x) },
         poly3 = { t <- seq_along(x); cyc <- resid(lm(x ~ t + I(t^2) + I(t^3)))
         as_cycle_with_residual_trend(cyc, x) },
         poly4 = { t <- seq_along(x); cyc <- resid(lm(x ~ t + I(t^2) + I(t^3) + I(t^4)))
         as_cycle_with_residual_trend(cyc, x) },
         poly5 = { t <- seq_along(x); cyc <- resid(lm(x ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5)))
         as_cycle_with_residual_trend(cyc, x) },
         linear = { t <- seq_along(x); cyc <- resid(lm(x ~ t))
         as_cycle_with_residual_trend(cyc, x) },
         log    = { t <- seq_along(x); cyc <- resid(lm(x ~ log(t)))
         as_cycle_with_residual_trend(cyc, x) },
         
         ## SSA on levels (returns BOTH, no pre-filtering!) 
         ssa_full = {
           .require_Rssa()
           dec <- .ssa_full_decompose_levels(
             x,
             L        = .get_cfg("ssa_full_L",                  32L),
             max_K    = .get_cfg("ssa_full_max_components",     20L),
             trend_min_period_obs = .get_cfg("ssa_trend_min_period_obs",    60L),
             trend_lowshare_min   = .get_cfg("ssa_trend_lowfreq_share_min", 0.80),
             pair_wcor_min        = .get_cfg("ssa_pair_wcor_min",           0.60),
             pair_centroid_tol    = .get_cfg("ssa_pair_centroid_tol",       0.10)
           )
           list(cycle = as.numeric(dec$cycle),
                trend = as.numeric(dec$trend),
                meta  = list(source = "ssa_full"))
         },
         
         ## No detrending 
         none = {
           list(cycle = as.numeric(x),
                trend = rep(0, length(x)),
                meta  = list(source = "none"))
         },
         
         stop("Unknown detrending method: ", method)
  )
}



# detrend_method <- "ideal20y"

# Baxter-King bandpass settings (lower=6 quarters, upper=32 quarters)
yrs_low  <- 1.5  # 6Q
yrs_high <- 8  # 32Q
bp_bounds <- c(1/(yrs_high*obs_per_year_real),
               1/(yrs_low *obs_per_year_real)) * 2


# Apply and score
res_ct <- get_cycle_trend(data_merged_real$log_hours, detrend_method)

data_merged_real$h     <- LOG * res_ct$cycle  # scale cycle *after* extraction (if logged TS: LOG=100, otherwise, =1)
data_merged_real$trend <- res_ct$trend        # always the proper trend for the method

cat("Method =", detrend_method,
    "→ corr(h, trend) =", round(cor(data_merged_real$h, data_merged_real$trend), 4), "\n")


# Re‑draw panels (a)/(b) exactly as before
p1 <- ggplot(data_merged_real, aes(date)) +
  geom_line(aes(y=log_hours), color="black", size=0.8) +
  geom_line(aes(y=trend),   color="red",   size=0.9) +
  labs(title="(a) Level and Trend", x="Date", y="") +
  theme_minimal()

p2 <- ggplot(data_merged_real, aes(date, h)) +
  geom_line(color="black", size=0.8) +
  labs(title="(b) Cyclical Component", x="Date", y="") +
  theme_minimal()

grid.arrange(p1, p2, ncol=2)


sum(is.na(data_merged_real$log_hours))   
which(is.na(data_merged_real$log_hours))  












########################
# UNIT ROOT TESTS on the cyclical component "h"
########################
head(data_merged_real$h)
summary(data_merged_real$h)

cycle_unit_root_results <- run_unit_root_tests(data_merged_real$h)

cat("=== ADF on cyclical component h ===\n")
print(cycle_unit_root_results$ADF)

cat("\n=== PP on cyclical component h ===\n")
print(cycle_unit_root_results$PP)

cat("\n=== KPSS on cyclical component h ===\n")
print(cycle_unit_root_results$KPSS)

cat("\n=== DF–GLS on cyclical component h ===\n")
print(cycle_unit_root_results$DFGLS)





##################
# Construct the Accumulation Series and Lagged Variables
##################

# Truncated accumulator: X_t^(N) = δ Σ_{j=0}^{N-1} (1-δ)^j x_{t-j}
calc_accumulation_truncated <- function(x, delta, N) {
  w <- delta * (1 - delta)^(0:(N - 1))
  as.numeric(stats::filter(x, w, method = "convolution", sides = 1))
}
data_merged_real <- data_merged_real[order(data_merged_real$date), ]
data_merged_real$X      <- calc_accumulation_truncated(data_merged_real$h, delta, N)
data_merged_real$X_lag1 <- dplyr::lag(data_merged_real$X, 1)

# AR lags for h
data_merged_real <- data_merged_real %>%
  mutate(
    x_lag1 = dplyr::lag(h, 1),
    x_lag2 = dplyr::lag(h, 2)
  )

# WINDOWING HELPER: subset AFTER features are built on the full timeline 
prepare_window <- function(df_full,
                           start_date = NULL,
                           end_date   = NULL,
                           N_required = 10 * obs_per_year_real) {
  # df_full already has h, X, X_lag1, x_lag1, x_lag2 on the full sample
  df_full <- df_full[order(df_full$date), ]
  
  # Default automated start: first date with enough pre-sample history
  # (N_required for the accumulator and 2 lags for AR)
  warmup <- max(N_required, 2)
  if (is.null(start_date)) {
    start_date <- df_full$date[warmup + 1L]
  }
  if (is.null(end_date)) {
    end_date <- max(df_full$date, na.rm = TRUE)
  }
  
  # 1) Subset to the requested window
  df_win <- subset(df_full, date >= start_date & date <= end_date)
  
  # 2) Drop only rows inside the window that still have missing lags
  df_win <- tidyr::drop_na(df_win, x_lag1, x_lag2, X_lag1)
  
  list(
    df_samp      = df_win,
    start_date_f = min(df_win$date),
    end_date_f   = max(df_win$date),
    horizon_data = nrow(df_win)
  )
}

## Build windowed sample once
win <- prepare_window(data_merged_real)   # auto-picks first feasible start
df_samp       <- win$df_samp
start_date_f  <- win$start_date_f
horizon_data  <- win$horizon_data
forecast_dates<- df_samp$date


forecast_origin <- df_samp$date[1+ forecast_starter]










###################################
# Create Nonlinear Transformation Variables
###################################
# Construct additional variables for various model specifications.
# Build a single regression row from the state (x_{t-1}, x_{t-2}, X_{t-1}).
# (AR(2), Linear, Minimal, Intermediate, Full), with no duplicate aliases.
make_row <- function(x_lag1, x_lag2, X_lag1, include_intercept = TRUE) {
  # input validation
  if (length(x_lag1) != 1L || length(x_lag2) != 1L || length(X_lag1) != 1L)
    stop("make_row expects scalars for x_lag1, x_lag2, and X_lag1.")
  if (any(!is.finite(c(x_lag1, x_lag2, X_lag1))))
    stop("Inputs must be finite numerics.")
  
  # primitives
  xl1 <- as.numeric(x_lag1)   # x_{t-1}
  xl2 <- as.numeric(x_lag2)   # x_{t-2}
  XL  <- as.numeric(X_lag1)   # X_{t-1}
  
  # powers 
  xl1_sq  <- xl1^2
  xl2_sq  <- xl2^2
  XL_sq   <- XL^2
  xl1_cub <- xl1^3
  xl2_cub <- xl2^3
  XL_cub  <- XL^3
  
  # assemble 
  out <- data.frame(
    # linear terms
    x_lag1  = xl1,
    x_lag2  = xl2,
    X_lag1  = XL,
    
    # quadratic terms
    x_lag1_sq = xl1_sq,
    x_lag2_sq = xl2_sq,
    X_lag1_sq = XL_sq,
    
    # cubic terms
    x_lag1_cub = xl1_cub,
    x_lag2_cub = xl2_cub,
    X_lag1_cub = XL_cub,
    
    # interactions (2nd & 3rd order)
    x_lag1_x_lag2     = xl1 * xl2,
    x_lag1_X_lag1     = xl1 * XL,
    x_lag2_X_lag1     = xl2 * XL,
    x_lag1_sq_x_lag2  = xl1_sq * xl2,
    x_lag1_x_lag2_sq  = xl1 * xl2_sq,
    x_lag1_sq_X_lag1  = xl1_sq * XL,   # keep *_lag1 naming
    x_lag2_sq_X_lag1  = xl2_sq * XL,
    x_lag1_X_lag1_sq  = xl1 * XL_sq,   # keep *_lag1 naming
    x_lag2_X_lag1_sq  = xl2 * XL_sq,
    x_lag1_x_lag2_X_lag1 = xl1 * xl2 * XL,
    check.names = FALSE
  )
  
  if (isTRUE(include_intercept)) {
    out <- cbind("(Intercept)" = 1, out, check.names = FALSE)
  }
  out
}



# Inspect the Final Processed Data
str(data_merged_real)
summary(data_merged_real)

# Plot the Processed Series
# Overlay the cyclical component h (in red) and the accumulation series X (in blue)
ggplot(data_merged_real, aes(x = date)) +
  geom_line(aes(y = h), color = "red", size = 1) +
  geom_line(aes(y = X), color = "blue", size = 1) +
  labs(title = "Processed Series: Cyclical Component h (Red) and Accumulated Series X (Blue)",
       x = "Date", y = "Value") +
  theme_minimal()









############################
# STAR models
############################

# helpers 
# Smooth transitions
G_logistic <- function(z, eta, c) {               # eta = log(gamma) ⇒ gamma = exp(eta) > 0
  gamma <- exp(eta)
  1 / (1 + exp(-gamma * (z - c)))
}
G_estar <- function(z, eta, c) {                  # ESTAR: symmetric around c
  gamma <- exp(eta)
  1 - exp(-gamma * (z - c)^2)
}

# Transition variable z_{t-1}. Start with z = x_{t-1} (can later pass a linear combo)
make_z <- function(df, z_var = "x_lag1", v = c(1,0,0)) {
  if (identical(z_var, "x_lag1")) return(df$x_lag1)
  if (identical(z_var, "combo"))  return(v[1]*df$x_lag1 + v[2]*df$x_lag2 + v[3]*df$X_lag1)
  stop("Unknown z_var: use 'x_lag1' or 'combo'")
}

# STAR specification evaluator (vectorized)
#   spec = "standard_lstar", "minimal_lstar", "minimal_estar"
#   theta layout:
#     standard_lstar:  [a0, a1, b10,b11, b20,b21, b30,b31, eta, c]
#     minimal_*     :  [a0,      b10,     b20,     b30,    d1, d2, d3, eta, c]
star_eval <- function(theta, x1, x2, X1, z, spec = c("standard_lstar","minimal_lstar","minimal_estar")) {
  spec <- match.arg(spec)
  if (spec == "standard_lstar") {
    a0 <- theta[1]; a1 <- theta[2]
    b10 <- theta[3];  b11 <- theta[4]
    b20 <- theta[5];  b21 <- theta[6]
    b30 <- theta[7];  b31 <- theta[8]
    eta <- theta[9];  c  <- theta[10]
    G <- G_logistic(z, eta, c)
    mu <- a0 + a1*G
    b1 <- b10 + b11*G
    b2 <- b20 + b21*G
    b3 <- b30 + b31*G
    return(mu + b1*x1 + b2*x2 + b3*X1)
  } else {
    a0 <- theta[1]
    b10 <- theta[2]; b20 <- theta[3]; b30 <- theta[4]
    d1  <- theta[5]; d2  <- theta[6]; d3  <- theta[7]
    eta <- theta[8]; c   <- theta[9]
    G   <- if (spec == "minimal_lstar") G_logistic(z, eta, c) else G_estar(z, eta, c)
    mu <- a0
    b1 <- b10 + d1*G
    b2 <- b20 + d2*G
    b3 <- b30 + d3*G
    return(mu + b1*x1 + b2*x2 + b3*X1)
  }
}


## STAR estimation
# Returns a list with: spec, coef (named), se (if Hessian ok), fitted, resid, logLik, AIC, BIC,
#                      z_var, v, gfun ("LSTAR"/"ESTAR")
estimate_star <- function(df,
                          spec = c("standard_lstar","minimal_lstar","minimal_estar"),
                          z_var = "x_lag1", v = c(1,0,0),
                          starts_gamma = c(0.5, 1, 2, 5, 10),
                          starts_c = c("median", "mean", "0"),
                          trace = TRUE, maxit = 2000, compute_se = FALSE) {
  spec <- match.arg(spec)
  
  # data used
  y  <- df$h
  x1 <- df$x_lag1
  x2 <- df$x_lag2
  X1 <- df$X_lag1
  z  <- make_z(df, z_var = z_var, v = v)
  
  ok <- is.finite(y) & is.finite(x1) & is.finite(x2) & is.finite(X1) & is.finite(z)
  y  <- y[ok]; x1 <- x1[ok]; x2 <- x2[ok]; X1 <- X1[ok]; z <- z[ok]
  n  <- length(y)
  stopifnot(n >= 50)
  
  # linear warm start 
  lin  <- stats::lm(y ~ x1 + x2 + X1)
  a0_  <- unname(coef(lin)[1])
  b10_ <- unname(coef(lin)["x1"]); if (is.na(b10_)) b10_ <- 0
  b20_ <- unname(coef(lin)["x2"]); if (is.na(b20_)) b20_ <- 0
  b30_ <- unname(coef(lin)["X1"]); if (is.na(b30_)) b30_ <- 0
  if (is.na(a0_)) a0_ <- 0
  
  # starting grid
  c0s <- switch(starts_c[1],
                median = median(z, na.rm = TRUE),
                mean   = mean(z,   na.rm = TRUE),
                `0`    = 0,
                median(z, na.rm = TRUE))
  inits <- list()
  if (spec == "standard_lstar") {
    for (g0 in starts_gamma) for (c0 in c(c0s, stats::quantile(z, c(.25,.75), na.rm = TRUE))) {
      inits[[length(inits) + 1L]] <- c(
        a0 = a0_, a1 = 0,
        b10 = b10_, b11 = 0,
        b20 = b20_, b21 = 0,
        b30 = b30_, b31 = 0,
        eta = log(g0), c = c0
      )
    }
  } else {
    for (g0 in starts_gamma) for (c0 in c(c0s, stats::quantile(z, c(.25,.75), na.rm = TRUE))) {
      inits[[length(inits) + 1L]] <- c(
        a0 = a0_,
        b10 = b10_, b20 = b20_, b30 = b30_,
        d1 = 0, d2 = 0, d3 = 0,
        eta = log(g0), c = c0
      )
    }
  }
  
  k   <- length(inits[[1L]])          # number of parameters
  se  <- rep(NA_real_, k)
  
  # objective
  obj <- function(theta) {
    yhat <- star_eval(theta, x1, x2, X1, z, spec)
    sum((y - yhat)^2)
  }
  
  # optimize over starts
  best <- list(val = Inf, par = NULL, conv = NA_integer_)
  for (st in inits) {
    opt <- try(stats::optim(st, obj, method = "BFGS",
                            control = list(maxit = maxit, reltol = 1e-10)), silent = TRUE)
    if (!inherits(opt, "try-error") && is.finite(opt$value) && opt$value < best$val) {
      best <- list(val = opt$value, par = opt$par, conv = opt$convergence)
      if (isTRUE(trace)) message(sprintf("STAR(%s): SSE = %.6f", spec, opt$value))
    }
  }
  
  theta_hat <- best$par
  if (is.null(theta_hat)) stop("STAR optimization failed for spec = ", spec)
  
  # fitted, resid, criteria
  yhat <- star_eval(theta_hat, x1, x2, X1, z, spec)
  res  <- y - yhat
  s2   <- mean(res^2)
  k    <- length(theta_hat)
  ll   <- -0.5 * n * (log(2*pi*s2) + 1)         # Gaussian QMLE
  AIC  <- n * log(s2) + 2 * k
  BIC  <- n * log(s2) + k * log(n)
  
  #  Jacobian-based SEs (sandwich under homoskedastic QMLE)
  if (compute_se && requireNamespace("numDeriv", quietly = TRUE)) {
    G <- try(numDeriv::jacobian(function(th) star_eval(th, x1, x2, X1, z, spec), theta_hat),
             silent = TRUE)
    if (!inherits(G, "try-error")) {
      JtJ <- crossprod(G)
      V   <- try(chol2inv(chol(JtJ)) * s2, silent = TRUE)
      if (!inherits(V, "try-error")) se <- sqrt(pmax(diag(V), 0))
    }
  }
  
  nm <- if (spec == "standard_lstar") {
    c("a0","a1","b10","b11","b20","b21","b30","b31","eta","c")
  } else {
    c("a0","b10","b20","b30","d1","d2","d3","eta","c")
  }
  names(theta_hat) <- nm
  names(se)        <- nm
  
  structure(list(
    spec   = spec,
    z_var  = z_var,
    v      = v,
    gfun   = if (spec == "minimal_estar") "ESTAR" else "LSTAR",
    coef   = theta_hat,
    se     = se,
    fitted = yhat,
    resid  = res,
    logLik = ll, AIC = AIC, BIC = BIC,
    nobs   = n
  ), class = "star_fit")
}

## Shared utilities for STAR testing/diagnostics 

# Small helpers
vec <- function(x) as.numeric(x)
rss <- function(e) sum(vec(e)^2)
safe_pkg <- function(p) { if (!requireNamespace(p, quietly=TRUE)) stop("Please install ", p) }

# Compute vcov(theta) for a fitted STAR 
vcov_star <- function(star_obj, df, tol = 1e-8, ridge = 0) {
  safe_pkg("numDeriv")
  # rebuild sample
  y  <- df$h; x1 <- df$x_lag1; x2 <- df$x_lag2; X1 <- df$X_lag1
  z  <- if (star_obj$z_var == "x_lag1") x1 else with(star_obj, v[1]*x_lag1 + v[2]*x_lag2 + v[3]*X_lag1)
  ok <- is.finite(y) & is.finite(x1) & is.finite(x2) & is.finite(X1) & is.finite(z)
  y  <- y[ok]; x1 <- x1[ok]; x2 <- x2[ok]; X1 <- X1[ok]; z <- z[ok]
  
  spec <- star_obj$spec
  f <- function(theta) star_eval(theta, x1, x2, X1, z, spec)
  
  # Jacobian wrt parameters
  G <- numDeriv::jacobian(f, star_obj$coef)   # n x k
  s2 <- mean((y - f(star_obj$coef))^2)
  
  JtJ <- crossprod(G)                          # k x k
  if (ridge > 0) {
    JtJ <- JtJ + diag(ridge * max(diag(JtJ), 1), ncol(JtJ))
  }
  
  V <- try(chol2inv(chol(JtJ)) * s2, silent = TRUE)
  if (inherits(V, "try-error")) {
    # Eigen-based Moore–Penrose pseudoinverse of J'J
    ev <- eigen(JtJ, symmetric = TRUE)
    lam <- ev$values
    P   <- ev$vectors
    if (!any(is.finite(lam))) stop("Jacobian has no finite eigenvalues (severe non-ID).")
    keep <- lam > (max(lam, na.rm = TRUE) * tol)
    rnk  <- sum(keep)
    if (rnk == 0) stop("Jacobian rank is zero (flat transition ⇒ non-identification).")
    JtJ_pinv <- P[, keep, drop=FALSE] %*% diag(1/lam[keep], nrow = rnk) %*% t(P[, keep, drop=FALSE])
    V <- JtJ_pinv * s2
    attr(V, "used_pinv") <- TRUE
    attr(V, "rank") <- rnk
  } else {
    attr(V, "used_pinv") <- FALSE
    attr(V, "rank") <- ncol(JtJ)
  }
  
  dimnames(V) <- list(names(star_obj$coef), names(star_obj$coef))
  V
}


# STAR restriction matrices for nested linearity
R_minimal_lstar <- function(theta_names) {
  # H0: d1 = d2 = d3 = 0   (nonlinearity switches off)
  idx <- match(c("d1","d2","d3"), theta_names)
  R <- matrix(0, nrow=3, ncol=length(theta_names))
  R[cbind(1:3, idx)] <- 1
  R
}
R_standard_lstar <- function(theta_names) {
  # H0: a1 = b11 = b21 = b31 = 0
  idx <- match(c("a1","b11","b21","b31"), theta_names)
  R <- matrix(0, nrow=4, ncol=length(theta_names))
  R[cbind(1:4, idx)] <- 1
  R
}

# Quasi-Wald (engineering) test; WARNING: unidentified (eta,c) under H0
wald_star <- function(star_obj, df, which = c("minimal_lstar","standard_lstar")) {
  which <- match.arg(which)
  V <- vcov_star(star_obj, df)
  th <- star_obj$coef
  if (which=="minimal_lstar") {
    R <- R_minimal_lstar(names(th)); q <- 3L
  } else {
    R <- R_standard_lstar(names(th)); q <- 4L
  }
  Rth  <- drop(R %*% th)
  RVRT <- R %*% V %*% t(R)
  stat <- drop(crossprod(Rth, solve(RVRT, Rth)))
  pval <- pchisq(stat, df=q, lower.tail=FALSE)
  list(stat=stat, df=q, p.value=pval, note="Wald under weak ID is only indicative; prefer bootstrap LR below.")
}

# Bootstrap LR for STAR nonlinearity (robust to (eta,c) nuisance) 
# Tests Linear vs STAR by parametric bootstrap: simulate from Linear (null), re-estimate STAR each draw.
star_bootstrap_LR <- function(df, star_spec = c("minimal_lstar","standard_lstar"),
                              B = 499, seed = 123) {
  star_spec <- match.arg(star_spec)
  set.seed(seed)
  # Fit null (Linear) 
  m0 <- lm(h ~ x_lag1 + x_lag2 + X_lag1, data=df)
  y  <- df$h
  e0 <- resid(m0)
  X  <- model.matrix(m0)
  yhat0 <- drop(X %*% coef(m0))
  sse0 <- sum(resid(m0)^2)
  
  # Fit STAR on original data to get observed LR
  m1 <- estimate_star(df, spec = if (star_spec=="minimal_lstar") "minimal_lstar" else "standard_lstar",
                      z_var="x_lag1", trace=FALSE)
  sse1 <- rss(m1$resid)
  n    <- length(y)
  LR_obs <- n * log(sse0/sse1)
  
  # Helper: simulate under null preserving AR structure
  sim_once <- function() {
    e <- sample(e0, replace=TRUE)
    y_sim <- yhat0 + e
    df_sim <- df
    # rebuild lagged features and accumulator consistently
    df_sim$h <- y_sim
    # rebuild X (accumulator) causally with delta 
    # only need to recompute lags
    df_sim$x_lag1 <- dplyr::lag(df_sim$h, 1)
    df_sim$x_lag2 <- dplyr::lag(df_sim$h, 2)
    # Use existing X_lag1 in df (constructed from original h) as an approximation.
    # For full rigor, recompute X from h using calc_accumulation_truncated() and delta/N.
    # (Left as-is to keep bootstrap lightweight.)
    # Refit STAR on simulated series
    fit_star <- try(estimate_star(df_sim, spec = if (star_spec=="minimal_lstar") "minimal_lstar" else "standard_lstar",
                                  z_var="x_lag1", trace=FALSE, compute_se = TRUE), silent=TRUE)
    if (inherits(fit_star, "try-error")) return(NA_real_)
    sse1b <- rss(fit_star$resid)
    nb    <- sum(is.finite(y_sim))
    nb * log(sum((y_sim - fitted(m0))^2)/sse1b) # approximate LR on sim
  }
  
  LRs <- replicate(B, sim_once())
  LRs <- LRs[is.finite(LRs)]
  pval <- mean(LRs >= LR_obs)
  list(LR_obs=LR_obs, LRs=LRs, p.value=pval, note="Parametric bootstrap LR: Linear (null) vs STAR (alt).")
}

#  Local stability for STAR
# Solve for steady state(s) h* and compute Jacobian eigenvalues at each.
star_steady_states <- function(star_obj, delta, grid = NULL) {
  th   <- star_obj$coef
  spec <- star_obj$spec
  zfun <- function(h) if (star_obj$z_var=="x_lag1") h else star_obj$v[1]*h + star_obj$v[2]*h + star_obj$v[3]*h
  
  Gfun <- switch(spec,
                 "standard_lstar" = function(h) 1/(1 + exp(-exp(th["eta"]) * (zfun(h) - th["c"]))),
                 "minimal_lstar"  = function(h) 1/(1 + exp(-exp(th["eta"]) * (zfun(h) - th["c"]))),
                 "minimal_estar"  = function(h) 1 - exp(-exp(th["eta"]) * (zfun(h) - th["c"])^2)
  )
  
  mu   <- function(h) if (spec=="standard_lstar") th["a0"] + th["a1"]*Gfun(h) else th["a0"]
  b1   <- function(h) if (spec=="standard_lstar") th["b10"] + th["b11"]*Gfun(h) else th["b10"] + th["d1"]*Gfun(h)
  b2   <- function(h) if (spec=="standard_lstar") th["b20"] + th["b21"]*Gfun(h) else th["b20"] + th["d2"]*Gfun(h)
  b3   <- function(h) if (spec=="standard_lstar") th["b30"] + th["b31"]*Gfun(h) else th["b30"] + th["d3"]*Gfun(h)
  
  Ffix <- function(h) mu(h) + (b1(h) + b2(h) + b3(h)) * h - h  # fixed-point eq: h = mu + (b1+b2+b3) h
  
  # search interval from sample scale (fallback if no df_samp here)
  if (is.null(grid)) grid <- seq(-3, 3, length.out = 601)
  vals <- sapply(grid, Ffix)
  sgn  <- sign(vals)
  idx  <- which(diff(sgn)!=0 & is.finite(vals[-length(vals)]) & is.finite(vals[-1]))
  roots <- numeric(0)
  for (j in idx) {
    lo <- grid[j]; hi <- grid[j+1]
    r  <- try(uniroot(Ffix, lower=lo, upper=hi)$root, silent=TRUE)
    if (!inherits(r,"try-error")) roots <- c(roots, unname(r))
  }
  roots <- unique(round(roots, 10))
  if (!length(roots)) roots <- 0  # at least inspect zero
  
  # Jacobian and eigenvalues at (h*, h*, H*), with H*=h* (steady state)
  J_of <- function(hs) {
    b1h <- b1(hs); b2h <- b2(hs); b3h <- b3(hs)
    matrix(c(
      b1h,         b2h,         b3h,
      1,           0,           0,
      delta*b1h,   delta*b2h,   (1-delta)+delta*b3h
    ), 3,3, byrow=TRUE)
  }
  out <- lapply(roots, function(hs) {
    J <- J_of(hs)
    ev <- eigen(J, only.values = TRUE)$values
    list(h_star=hs, J=J, eigenvalues=ev,
         rho = max(Mod(ev)),
         imag_max = max(abs(Im(ev))))
  })
  out
}





########################################################
##### MODEL ESTIMATION
########################################################
add_features <- function(df) {
  df %>%
    mutate(
      # squares
      x_lag1_sq = x_lag1^2,
      x_lag2_sq = x_lag2^2,
      X_lag1_sq = X_lag1^2,
      # cubes
      x_lag1_cub = x_lag1^3,
      x_lag2_cub = x_lag2^3,
      X_lag1_cub = X_lag1^3,
      # interactions (canonical names)
      x_lag1_x_lag2        = x_lag1 * x_lag2,
      x_lag1_X_lag1        = x_lag1 * X_lag1,
      x_lag2_X_lag1        = x_lag2 * X_lag1,
      x_lag1_sq_x_lag2     = x_lag1_sq * x_lag2,
      x_lag1_x_lag2_sq     = x_lag1 * x_lag2_sq,
      x_lag1_sq_X_lag1     = x_lag1_sq * X_lag1,
      x_lag2_sq_X_lag1     = x_lag2_sq * X_lag1,
      x_lag1_X_lag1_sq     = x_lag1 * X_lag1_sq,
      x_lag2_X_lag1_sq     = x_lag2 * X_lag1_sq,
      x_lag1_x_lag2_X_lag1 = x_lag1 * x_lag2 * X_lag1
    )
}

## Now add features on the windowed data
df_samp <- add_features(df_samp)

# quick guard 
need <- c("x_lag1","x_lag2","X_lag1",
          "x_lag1_sq","x_lag2_sq","X_lag1_sq",
          "x_lag1_cub","x_lag2_cub","X_lag1_cub",
          "x_lag1_x_lag2","x_lag1_X_lag1","x_lag2_X_lag1",
          "x_lag1_sq_x_lag2","x_lag1_x_lag2_sq",
          "x_lag1_sq_X_lag1","x_lag2_sq_X_lag1",
          "x_lag1_X_lag1_sq","x_lag2_X_lag1_sq",
          "x_lag1_x_lag2_X_lag1")
stopifnot(all(need %in% names(df_samp)))


# Model 1: AR(2) Model
model_ar2 <- lm(h ~ x_lag1 + x_lag2, data = df_samp)
summary(model_ar2)

# Model 2: Linear Model (AR(2) + Accumulation Term)
model_linear <- lm(h ~ x_lag1 + x_lag2 + X_lag1, data = df_samp)
summary(model_linear)

# Model 3: Minimal Nonlinear Model (adds cubic term in x_lag1)
model_minimal <- lm(h ~ x_lag1 + x_lag2 + X_lag1 + x_lag1_cub, data = df_samp)
summary(model_minimal)

# Model 4: Intermediate Nonlinear Model
# Includes cubic terms in x_lag1 and x_lag2, cubic term in accumulation,
# and two cross-products.
model_intermediate <- lm(
  h ~ x_lag1 + x_lag2 + X_lag1 +
    x_lag1_cub + x_lag2_cub + X_lag1_cub +
    x_lag1_sq_X_lag1 + x_lag1_X_lag1_sq,
  data = df_samp)

summary(model_intermediate)

# Model 5: Full Nonlinear Model
# Includes all linear, quadratic, and cubic terms and interactions.
model_full <- lm(h ~ x_lag1 + x_lag2 + X_lag1 +
                   x_lag1_sq + x_lag2_sq + X_lag1_sq +
                   x_lag1_x_lag2 + x_lag1_X_lag1 + x_lag2_X_lag1 +
                   x_lag1_cub + x_lag2_cub + X_lag1_cub +
                   x_lag1_sq_x_lag2 + x_lag1_x_lag2_sq +
                   x_lag1_sq_X_lag1 + x_lag2_sq_X_lag1 +
                   x_lag1_X_lag1_sq + x_lag2_X_lag1_sq +
                   x_lag1_x_lag2_X_lag1,
                 data = df_samp)
summary(model_full)




#  STAR models (share same regressors as Linear/Minimal etc.) 
# z = x_{t-1} to start;
model_star_std  <- estimate_star(df_samp, spec = "standard_lstar", z_var = "x_lag1")
model_star_min  <- estimate_star(df_samp, spec = "minimal_lstar",  z_var = "x_lag1")
model_estar_min <- estimate_star(df_samp, spec = "minimal_estar",  z_var = "x_lag1")

print(model_star_std)
print(model_star_min)
print(model_estar_min)







########################################################
############## DIAGNOSTICS for all models (LM + STAR variants)
########################################################
if (!requireNamespace("lmtest", quietly = TRUE)) stop("Please install.packages('lmtest')")
if (!requireNamespace("tseries", quietly = TRUE)) stop("Please install.packages('tseries')")

# Collect all models by name
models <- list(
  AR2          = model_ar2,
  Linear       = model_linear,
  Minimal      = model_minimal,
  Intermediate = model_intermediate,
  Full         = model_full,
  STAR_std     = model_star_std,
  STAR_min     = model_star_min,
  ESTAR_min    = model_estar_min
)

# Helper: compute DW on arbitrary residual vector
dw_stat_from_resid <- function(e) sum(diff(e)^2, na.rm=TRUE) / sum(e^2, na.rm=TRUE)

# Helper: Ljung–Box p-value at specified lag
lb_p <- function(e, lag = 20) {
  tryCatch({
    Box.test(e, lag = lag, type = "Ljung-Box")$p.value
  }, error = function(...) NA_real_)
}

# Build summary row for a single model
summarize_model <- function(name, mod) {
  is_lm   <- inherits(mod, "lm")
  is_star <- inherits(mod, "star_fit")
  if (!is_lm && !is_star) stop("Unknown model class for: ", name)
  
  # Residuals & fitted
  resid_vec  <- if (is_lm) residuals(mod) else mod$resid
  fitted_vec <- if (is_lm) fitted(mod)    else mod$fitted
  resid_vec  <- as.numeric(resid_vec)
  fitted_vec <- as.numeric(fitted_vec)
  n <- sum(is.finite(resid_vec))
  
  # R2 / Adj-R2
  if (is_lm) {
    s <- summary(mod)
    R2    <- unname(s$r.squared)
    adjR2 <- unname(s$adj.r.squared)
    k     <- length(coef(mod)) - 1L
  } else {
    fs <- fit_stats_from_resid(fitted_vec, resid_vec, k = length(mod$coef))
    R2    <- fs$R2
    adjR2 <- fs$Adj_R2
    k     <- length(mod$coef)
  }
  
  # AIC / BIC
  AICv <- if (is_lm) AIC(mod) else mod$AIC
  BICv <- if (is_lm) BIC(mod) else mod$BIC
  
  # Durbin–Watson
  if (is_lm) {
    dw <- tryCatch({ lmtest::dwtest(mod) }, error = function(...) NULL)
    DW_stat <- if (is.null(dw)) dw_stat_from_resid(resid_vec) else unname(dw$statistic)
    DW_p    <- if (is.null(dw)) NA_real_                       else unname(dw$p.value)
  } else {
    DW_stat <- dw_stat_from_resid(resid_vec)
    DW_p    <- NA_real_   # no valid p-value under nonlinear STAR
  }
  
  # Breusch–Godfrey up to lag 4 (LM only)
  if (is_lm) {
    bg <- tryCatch({ lmtest::bgtest(mod, order = 4) }, error = function(...) NULL)
    BG_p <- if (is.null(bg)) NA_real_ else unname(bg$p.value)
  } else {
    BG_p <- NA_real_
  }
  
  # Ljung–Box at lags 12 and 20 (on residuals)
  LB12_p <- lb_p(resid_vec, lag = 12)
  LB20_p <- lb_p(resid_vec, lag = 20)
  
  # Normality: Shapiro–Wilk (only reliable up to 5000 obs)
  SW_p <- if (n >= 3 && n <= 5000) {
    tryCatch(stats::shapiro.test(resid_vec)$p.value, error = function(...) NA_real_)
  } else NA_real_
  
  # Jarque–Bera 
  JB_p <- tryCatch(tseries::jarque.bera.test(resid_vec)$p.value, error = function(...) NA_real_)
  
  data.frame(
    model = name,
    n = n,
    k = k,
    R2 = R2, Adj_R2 = adjR2,
    AIC = AICv, BIC = BICv,
    DW = as.numeric(DW_stat), DW_p = as.numeric(DW_p),
    BG_p = as.numeric(BG_p),
    LB12_p = as.numeric(LB12_p),
    LB20_p = as.numeric(LB20_p),
    Shapiro_p = as.numeric(SW_p),
    JB_p = as.numeric(JB_p),
    class = if (is_lm) "lm" else "STAR",
    stringsAsFactors = FALSE
  )
}

#  Helper: R2, Adj-R2, and DW from fitted+resid 
fit_stats_from_resid <- function(fitted_vec, resid_vec, k) {
  y <- as.numeric(fitted_vec) + as.numeric(resid_vec)
  e <- as.numeric(resid_vec)
  n <- sum(is.finite(y))
  
  sse <- sum(e^2, na.rm = TRUE)
  sst <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  
  R2 <- if (sst > 0) 1 - sse / sst else NA_real_
  Adj_R2 <- if (is.finite(R2) && n > k + 1) 1 - (sse / (n - k)) / (sst / (n - 1)) else NA_real_
  
  # Durbin–Watson from residuals
  DW <- if (sse > 0) sum(diff(e)^2, na.rm = TRUE) / sse else NA_real_
  
  list(R2 = R2, Adj_R2 = Adj_R2, DW = DW)
}

lb_pval <- function(e, lag = 20) {
  tryCatch(Box.test(e, lag = lag, type = "Ljung-Box")$p.value,
           error = function(...) NA_real_)
}

diag_table <- do.call(rbind, Map(summarize_model, names(models), models))

#  printing
diag_table_print <- within(diag_table, {
  R2        <- round(R2, 4)
  Adj_R2    <- round(Adj_R2, 4)
  AIC       <- round(AIC, 2)
  BIC       <- round(BIC, 2)
  DW        <- round(DW, 3)
  DW_p      <- signif(DW_p, 3)
  BG_p      <- signif(BG_p, 3)
  LB12_p    <- signif(LB12_p, 3)
  LB20_p    <- signif(LB20_p, 3)
  Shapiro_p <- signif(Shapiro_p, 3)
  JB_p      <- signif(JB_p, 3)
})

print(diag_table_print, row.names = FALSE)








################## STAR Restriction test
# STAR restriction tests 
# (i) Indicative Wald (fast, but eta/c not ID under H0)
wald_std <- wald_star(model_star_std,  df_samp, which="standard_lstar")
wald_min <- wald_star(model_star_min,  df_samp, which="minimal_lstar")
wald_est <- wald_star(model_estar_min, df_samp, which="minimal_lstar")  # same restrictions d1=d2=d3=0

cat("\n[Indicative] Wald Standard LSTAR (a1=b11=b21=b31=0):  chi2(",
    wald_std$df, ")=", round(wald_std$stat,3), "  p=", signif(wald_std$p.value,3), "\n", sep="")
cat("[Indicative] Wald Minimal LSTAR (d1=d2=d3=0):         chi2(",
    wald_min$df, ")=", round(wald_min$stat,3), "  p=", signif(wald_min$p.value,3), "\n", sep="")
cat("[Indicative] Wald ESTAR-min (d1=d2=d3=0):             chi2(",
    wald_est$df, ")=", round(wald_est$stat,3), "  p=", signif(wald_est$p.value,3), "\n", sep="")




##############################################################
# WALD TEST
##############################################################

###############################
# 1) AR(2) vs. Linear
###############################
# AR(2) Model lacks X_lag1.
hyp_ar2_linear <- c("X_lag1 = 0")
wald_ar2_linear <- linearHypothesis(lm(h ~ x_lag1 + x_lag2 + X_lag1, data = df_samp), 
                                    hyp_ar2_linear)
cat("Wald Test: AR(2) vs. Linear\n")
print(wald_ar2_linear)

###############################
# 2) AR(2) vs. Minimal 
###############################
# AR(2) lacks X_lag1 and x_lag1_cub.
hyp_ar2_minimal <- c("X_lag1 = 0", "x_lag1_cub = 0")
wald_ar2_minimal <- linearHypothesis(lm(h ~ x_lag1 + x_lag2 + X_lag1 + x_lag1_cub, data = df_samp),
                                     hyp_ar2_minimal)
cat("Wald Test: AR(2) vs. Minimal\n")
print(wald_ar2_minimal)

###############################
# 3) AR(2) vs. Intermediate 
###############################
# AR(2) lacks X_lag1, x_lag1_cub, and also x_lag2_cub, X_lag1_cub, x_lag1_sq_X, x_lag1_X_sq.
hyp_ar2_intermed <- c(
  "X_lag1 = 0", "x_lag1_cub = 0", "x_lag2_cub = 0", "X_lag1_cub = 0",
  "x_lag1_sq_X_lag1 = 0", "x_lag1_X_lag1_sq = 0")
wald_ar2_intermed <- linearHypothesis(
  lm(h ~ x_lag1 + x_lag2 + X_lag1 +
       x_lag1_cub + x_lag2_cub + X_lag1_cub +
       x_lag1_sq_X_lag1 + x_lag1_X_lag1_sq,
     data = df_samp),
  hyp_ar2_intermed)

cat("Wald Test: AR(2) vs. Intermediate\n")
print(wald_ar2_intermed)

###############################
# 4) AR(2) vs. Full
###############################
# AR(2) lacks X_lag1 and all nonlinear terms in the Full model.
hyp_ar2_full <- c( 
  "X_lag1 = 0",                    # Accumulation term
  # Quadratic terms:
  "x_lag1_sq = 0", "x_lag2_sq = 0", "X_lag1_sq = 0", 
  "x_lag1_x_lag2 = 0", "x_lag1_X_lag1 = 0", "x_lag2_X_lag1 = 0", 
  # Cubic terms:
  "x_lag1_cub = 0", "x_lag2_cub = 0", "X_lag1_cub = 0", 
  # Mixed interaction terms:
  "x_lag1_sq_x_lag2 = 0", "x_lag1_x_lag2_sq = 0", 
  "x_lag1_sq_X_lag1 = 0", "x_lag2_sq_X_lag1 = 0", 
  "x_lag1_X_lag1_sq = 0", "x_lag2_X_lag1_sq = 0", 
  "x_lag1_x_lag2_X_lag1 = 0"
)
wald_ar2_full <- linearHypothesis(lm(h ~ x_lag1 + x_lag2 + X_lag1 +
                                       x_lag1_sq + x_lag2_sq + X_lag1_sq +
                                       x_lag1_x_lag2 + x_lag1_X_lag1 + x_lag2_X_lag1 +
                                       x_lag1_cub + x_lag2_cub + X_lag1_cub +
                                       x_lag1_sq_x_lag2 + x_lag1_x_lag2_sq +
                                       x_lag1_sq_X_lag1 + x_lag2_sq_X_lag1 +
                                       x_lag1_X_lag1_sq + x_lag2_X_lag1_sq +
                                       x_lag1_x_lag2_X_lag1, data = df_samp),
                                  hyp_ar2_full)
cat("Wald Test: AR(2) vs. Full\n")
print(wald_ar2_full)

###############################
# 5) Linear vs. Minimal
###############################
# Minimal adds x_lag1_cub relative to Linear.
hyp_linear_minimal <- c("x_lag1_cub = 0")
wald_linear_minimal <- linearHypothesis(lm(h ~ x_lag1 + x_lag2 + X_lag1 + x_lag1_cub, data = df_samp),
                                        hyp_linear_minimal)
cat("Wald Test: Linear vs. Minimal\n")
print(wald_linear_minimal)

###############################
# 6) Linear vs. Intermediate 
###############################
# Intermediate adds (x_lag1_cub, x_lag2_cub, X_lag1_cub, x_lag1_sq_X, x_lag1_X_sq) beyond Linear.
hyp_linear_intermed <- c(
  "x_lag1_cub = 0", "x_lag2_cub = 0", "X_lag1_cub = 0",
  "x_lag1_sq_X_lag1 = 0", "x_lag1_X_lag1_sq = 0")
wald_linear_intermed <- linearHypothesis(
  lm(h ~ x_lag1 + x_lag2 + X_lag1 +
       x_lag1_cub + x_lag2_cub + X_lag1_cub +
       x_lag1_sq_X_lag1 + x_lag1_X_lag1_sq,
     data = df_samp),
  hyp_linear_intermed)
cat("Wald Test: Linear vs. Intermediate\n")
print(wald_linear_intermed)

###############################
# 7) Linear vs. Full 
###############################
# Linear has x_lag1, x_lag2, X_lag1. Test all other terms in Full.
hyp_linear_full <- c("x_lag1_sq = 0", "x_lag2_sq = 0", "X_lag1_sq = 0",
                     "x_lag1_x_lag2 = 0", "x_lag1_X_lag1 = 0", "x_lag2_X_lag1 = 0",
                     "x_lag1_cub = 0", "x_lag2_cub = 0", "X_lag1_cub = 0",
                     "x_lag1_sq_x_lag2 = 0", "x_lag1_x_lag2_sq = 0",
                     "x_lag1_sq_X_lag1 = 0", "x_lag2_sq_X_lag1 = 0",
                     "x_lag1_X_lag1_sq = 0", "x_lag2_X_lag1_sq = 0",
                     "x_lag1_x_lag2_X_lag1 = 0")
wald_linear_full <- linearHypothesis(lm(h ~ x_lag1 + x_lag2 + X_lag1 +
                                          x_lag1_sq + x_lag2_sq + X_lag1_sq +
                                          x_lag1_x_lag2 + x_lag1_X_lag1 + x_lag2_X_lag1 +
                                          x_lag1_cub + x_lag2_cub + X_lag1_cub +
                                          x_lag1_sq_x_lag2 + x_lag1_x_lag2_sq +
                                          x_lag1_sq_X_lag1 + x_lag2_sq_X_lag1 +
                                          x_lag1_X_lag1_sq + x_lag2_X_lag1_sq +
                                          x_lag1_x_lag2_X_lag1, data = df_samp),
                                     hyp_linear_full)
cat("Wald Test: Linear vs. Full\n")
print(wald_linear_full)

###############################
# 8) Minimal vs. Intermediate 
###############################
# Minimal has x_lag1_cub; Intermediate adds x_lag2_cub, X_lag1_cub, x_lag1_sq_X, x_lag1_X_sq.
hyp_minimal_intermed <- c(
  "x_lag2_cub = 0", "X_lag1_cub = 0",
  "x_lag1_sq_X_lag1 = 0", "x_lag1_X_lag1_sq = 0")
wald_minimal_intermed <- linearHypothesis(
  lm(h ~ x_lag1 + x_lag2 + X_lag1 +
       x_lag1_cub + x_lag2_cub + X_lag1_cub +
       x_lag1_sq_X_lag1 + x_lag1_X_lag1_sq,
     data = df_samp),
  hyp_minimal_intermed)

cat("Wald Test: Minimal vs. Intermediate\n")
print(wald_minimal_intermed)

###############################
# 9) Minimal vs. Full 
###############################
# Minimal vs. Full: Minimal has x_lag1_cub; Full includes additionally:
hyp_minimal_full <- c("x_lag1_sq = 0", "x_lag2_sq = 0", "X_lag1_sq = 0",
                      "x_lag1_x_lag2 = 0", "x_lag1_X_lag1 = 0", "x_lag2_X_lag1 = 0",
                      "x_lag2_cub = 0", "X_lag1_cub = 0",
                      "x_lag1_sq_x_lag2 = 0", "x_lag1_x_lag2_sq = 0",
                      "x_lag1_sq_X_lag1 = 0", "x_lag2_sq_X_lag1 = 0",
                      "x_lag1_X_lag1_sq = 0", "x_lag2_X_lag1_sq = 0",
                      "x_lag1_x_lag2_X_lag1 = 0")
wald_minimal_full <- linearHypothesis(lm(h ~ x_lag1 + x_lag2 + X_lag1 +
                                           x_lag1_sq + x_lag2_sq + X_lag1_sq +
                                           x_lag1_x_lag2 + x_lag1_X_lag1 + x_lag2_X_lag1 +
                                           x_lag1_cub + x_lag2_cub + X_lag1_cub +
                                           x_lag1_sq_x_lag2 + x_lag1_x_lag2_sq +
                                           x_lag1_sq_X_lag1 + x_lag2_sq_X_lag1 +
                                           x_lag1_X_lag1_sq + x_lag2_X_lag1_sq +
                                           x_lag1_x_lag2_X_lag1, data = df_samp),
                                      hyp_minimal_full)
cat("Wald Test: Minimal vs. Full\n")
print(wald_minimal_full)

###############################
# 10) Intermediate vs. Full 
###############################
# Intermediate lacks everything in Full except for x_lag1_cub, x_lag2_cub, X_lag1_cub,
# x_lag1_sq_X, and x_lag1_X_sq. Test all remaining terms.
hyp_intermed_full <- c("x_lag1_sq = 0", "x_lag2_sq = 0", "X_lag1_sq = 0",
                       "x_lag1_x_lag2 = 0", "x_lag1_X_lag1 = 0", "x_lag2_X_lag1 = 0",
                       "x_lag1_sq_x_lag2 = 0", "x_lag1_x_lag2_sq = 0",
                       "x_lag1_sq_X_lag1 = 0", "x_lag2_sq_X_lag1 = 0",
                       "x_lag1_X_lag1_sq = 0", "x_lag2_X_lag1_sq = 0",
                       "x_lag1_x_lag2_X_lag1 = 0")
wald_intermed_full <- linearHypothesis(lm(h ~ x_lag1 + x_lag2 + X_lag1 +
                                            x_lag1_sq + x_lag2_sq + X_lag1_sq +
                                            x_lag1_x_lag2 + x_lag1_X_lag1 + x_lag2_X_lag1 +
                                            x_lag1_cub + x_lag2_cub + X_lag1_cub +
                                            x_lag1_sq_x_lag2 + x_lag1_x_lag2_sq +
                                            x_lag1_sq_X_lag1 + x_lag2_sq_X_lag1 +
                                            x_lag1_X_lag1_sq + x_lag2_X_lag1_sq +
                                            x_lag1_x_lag2_X_lag1, data = df_samp),
                                       hyp_intermed_full)
cat("Wald Test: Intermediate vs. Full\n")
print(wald_intermed_full)






#############################################################
## LASSO with AIC / BIC / AICc / CV selector
#############################################################
#  0) SETTINGS
response_col   <- "h"       
selection_mode <- "CV"      # choices: "AIC", "BIC", "AICc", "CV"

# Build the Full design (no intercept column!)
X <- as.matrix(df_samp[, c(
  "x_lag1","x_lag2","X_lag1",
  "x_lag1_sq","x_lag2_sq","X_lag1_sq",
  "x_lag1_x_lag2","x_lag1_X_lag1","x_lag2_X_lag1",
  "x_lag1_cub","x_lag2_cub","X_lag1_cub",
  "x_lag1_sq_x_lag2","x_lag1_x_lag2_sq",
  "x_lag1_sq_X_lag1","x_lag2_sq_X_lag1",
  "x_lag1_X_lag1_sq","x_lag2_X_lag1_sq",
  "x_lag1_x_lag2_X_lag1"
)])

y <- df_samp$h
cc <- complete.cases(X, y); X <- X[cc, , drop=FALSE]; y <- y[cc]

# 1) Design matrix (no intercept), manual standardization 
X_raw <- model.matrix(model_full, data = df_samp)[, -1, drop = FALSE]
y_all <- df_samp[[response_col]]

# keep only complete rows
cc    <- stats::complete.cases(X_raw, y_all)
X_raw <- X_raw[cc, , drop = FALSE]
y     <- y_all[cc]

# scale columns (avoid divide-by-zero)
scaler <- apply(X_raw, 2, sd)
scaler[!is.finite(scaler) | scaler == 0] <- 1
X <- sweep(X_raw, 2, scaler, "/")
n <- length(y)

# 2) Fit a dense LASSO path 
fit <- glmnet(
  x = X, y = y,
  alpha = 1,
  standardize = FALSE,   # we already scaled X
  intercept   = TRUE,    # unpenalized intercept
  nlambda     = 3000,
  lambda.min.ratio = max(1e-6, 0.01/length(y))
)

# 3) Information-criteria along the path 
preds <- predict(fit, newx = X)             # n x nlambda
rss   <- colSums((y - preds)^2)             # length = nlambda
df    <- fit$df + 1L                        # +1 for intercept
aic   <- n * log(rss / n) + 2 * df
bic   <- n * log(rss / n) + log(n) * df

# AICc is undefined when n - df - 1 <= 0; set those to +Inf
den   <- pmax(n - df - 1, 1e-9)
aicc  <- aic + 2 * df * (df + 1) / den
aicc[ n - df - 1 <= 0 ] <- Inf

lambda_aic  <- fit$lambda[ which.min(aic)  ]
lambda_bic  <- fit$lambda[ which.min(bic)  ]
lambda_aicc <- fit$lambda[ which.min(aicc) ]

cat(sprintf("λ (AIC)=%.6g, λ (BIC)=%.6g, λ (AICc)=%.6g\n",
            lambda_aic, lambda_bic, lambda_aicc))

# 4) Blocked K-fold CV (time-series friendly) 
K <- 10
foldid <- cut(seq_len(nrow(X)), breaks = K, labels = FALSE)

cvfit <- cv.glmnet(
  x = X, y = y,
  alpha = 1,
  standardize = FALSE,
  intercept   = TRUE,
  nlambda     = 5000,
  lambda.min.ratio = max(1e-6, 0.01/length(y)),
  foldid = foldid
)
lambda_cv <- cvfit$lambda.min
cat(sprintf("λ (CV.min)=%.6g  (λ.1se=%.6g)\n", cvfit$lambda.min, cvfit$lambda.1se))

# 5) Helper: Post-LASSO OLS at a chosen λ
get_post_lasso <- function(lambda_use) {
  bh <- as.matrix(coef(fit, s = lambda_use))
  surv <- setdiff(rownames(bh)[bh[, 1] != 0], "(Intercept)")
  if (length(surv) == 0) {
    # intercept-only fallback
    mod <- lm(reformulate(NULL, response = response_col), data = df_samp[cc, , drop = FALSE])
    return(list(model = mod, survivors = character(0), lambda = lambda_use))
  }
  form <- reformulate(surv, response = response_col)
  mod  <- lm(form, data = df_samp[cc, , drop = FALSE])
  list(model = mod, survivors = surv, lambda = lambda_use)
}

# 6) Pick which λ to use and run post-LASSO OLS
lambda_pick <- switch(toupper(selection_mode),
                      "AIC"  = lambda_aic,
                      "BIC"  = lambda_bic,
                      "AICc" = lambda_aicc,
                      "CV"   = lambda_cv,
                      stop("selection_mode must be one of: 'AIC','BIC','AICc','CV'")
)

post <- get_post_lasso(lambda_pick)

cat(sprintf("\nSelection = %s  →  λ=%.6g\n", toupper(selection_mode), post$lambda))
cat("Survivors:\n"); print(post$survivors)

# quick summary
sm <- summary(post$model)
print(sm$coefficients)
cat(sprintf("\nR^2 = %.3f   Adj.R^2 = %.3f\n", sm$r.squared, sm$adj.r.squared))

# (optional) also keep a CV-based and an AIC-based model side-by-side
post_AIC <- get_post_lasso(lambda_aic)
post_CV  <- get_post_lasso(lambda_cv)

model_lasso <- post$model        # or: post_AIC$model, post_CV$model

coef_lasso_vec <- coef(model_lasso)
lasso_coef_mat <- matrix(coef_lasso_vec,
                         ncol = 1,
                         dimnames = list(names(coef_lasso_vec), NULL))









######### Summaries: R², Adjusted R² and Durbin-Watson Test ######### 

# First, for the models estimated via lm (AR(2), Linear, Minimal, Intermediate, Full):
model_list <- list(
  "AR(2)" = model_ar2,
  "Linear" = model_linear,
  "Minimal" = model_minimal,
  "Intermediate" = model_intermediate,
  "Full" = model_full
)

results_df <- data.frame(
  Model = character(),
  R2 = numeric(),
  Adj_R2 = numeric(),
  DW = numeric(),
  stringsAsFactors = FALSE
)

for (nm in names(model_list)) {
  m <- model_list[[nm]]
  ss <- summary(m)
  r2_val <- ss$r.squared
  adjr2_val <- ss$adj.r.squared
  # Use dwtest from lmtest package:
  dw_val <- dwtest(m, alternative = "two.sided")$statistic
  results_df <- rbind(results_df, data.frame(
    Model = nm,
    R2 = round(r2_val, 3),
    Adj_R2 = round(adjr2_val, 3),
    DW = round(dw_val, 2),
    stringsAsFactors = FALSE
  ))
}


# 1) Residual diagnostics & fit stats 

star_models <- list("STAR-std"=model_star_std, "STAR-min"=model_star_min, "ESTAR-min"=model_estar_min)
for (nm in names(star_models)) {
  mm <- star_models[[nm]]
  k  <- length(mm$coef)
  fs <- fit_stats_from_resid(mm$fitted, mm$resid, k)
  results_df <- rbind(
    results_df,
    data.frame(Model = nm,
               R2    = round(fs$R2, 3),
               Adj_R2= round(fs$Adj_R2, 3),
               DW    = round(fs$DW, 2),
               stringsAsFactors = FALSE)
  )
  cat("\n", nm, ": Ljung-Box p(20) =", round(lb_pval(mm$resid, lag=20), 3), "\n")
}
print(results_df)

# FOR LASSO: 
resid_lasso <- resid(model_lasso)
sum_lasso   <- summary(model_lasso)
R2_lasso <- sum_lasso$r.squared
n_lasso <- length(resid_lasso)                     # number of observations
p_lasso <- length(coef(model_lasso)) - 1           # number of predictors (drop intercept)
# Now compute Adjusted R² and Durbin–Watson
Adj_R2_lasso <- 1 - (1 - R2_lasso) * ((n_lasso - 1) / (n_lasso - p_lasso - 1))
DW_lasso     <- sum(diff(resid_lasso)^2) / sum(resid_lasso^2)
# Append to existing results_df
results_df <- rbind(results_df, data.frame(
  Model   = "LASSO",
  R2      = round(R2_lasso,  3),
  Adj_R2  = round(Adj_R2_lasso, 3),
  DW      = round(DW_lasso,   2),
  stringsAsFactors = FALSE
))


print(results_df)








################################################################################
# LOCAL STABILITY & SS ANALYSIS 
################################################################################

# Maximum Eigenvalue of the Local Linear Approximation
`%||%` <- function(a,b) if (length(a)==0||is.na(a)) b else a
ρ       <- function(M)   max(Mod(eigen(M,only.values=TRUE)$values))
eig <- function(M) eigen(M, only.values = TRUE)$values


# 1. ONE-ROW structure every regressor 
make_row <- function(x) {
  xl1 <- x; xl2 <- x; XL <- x
  xl1_sq <- xl1^2; xl2_sq <- xl2^2; XL_sq <- XL^2
  xl1_cub <- xl1^3; xl2_cub <- xl2^3; XL_cub <- XL^3
  
  data.frame(
    x_lag1  = xl1,
    x_lag2  = xl2,
    X_lag1  = XL,
    # quadratic
    x_lag1_sq = xl1_sq,
    x_lag2_sq = xl2_sq,
    X_lag1_sq = XL_sq,
    # cubic
    x_lag1_cub = xl1_cub,
    x_lag2_cub = xl2_cub,
    X_lag1_cub = XL_cub,
    # interactions (canonical)
    x_lag1_x_lag2       = xl1 * xl2,
    x_lag1_X_lag1       = xl1 * XL,
    x_lag2_X_lag1       = xl2 * XL,
    x_lag1_sq_x_lag2    = xl1_sq * xl2,
    x_lag1_x_lag2_sq    = xl1 * xl2_sq,
    x_lag1_sq_X_lag1    = xl1_sq * XL,
    x_lag2_sq_X_lag1    = xl2_sq * XL,
    x_lag1_X_lag1_sq    = xl1 * XL_sq,
    x_lag2_X_lag1_sq    = xl2 * XL_sq,
    x_lag1_x_lag2_X_lag1= xl1 * xl2 * XL
  )
}



# 2  Put the six specifications in a list
specs <- list(
  AR2          = model_ar2,
  Linear       = model_linear,
  Minimal      = model_minimal,
  Intermediate = model_intermediate,
  Full         = model_full,
  LASSO        = coef(model_lasso)          # numeric vector
)

# 3  Search interval for roots  (μ ± 5σ of the sample)
xv <- df_samp$h[is.finite(df_samp$h)]
μ   <- if (length(xv)) mean(xv) else 0
σ   <- if (length(xv)&&sd(xv)>0) sd(xv) else 1
lo  <- μ - 5*σ
hi  <- μ + 5*σ
grid<- seq(lo, hi, length.out = 601)


# 4   Main loop

for (nm in names(specs)) {
  
  cat("\n================  MODEL:", nm, "  ================\n")
  
  obj      <- specs[[nm]]
  is_vec   <- is.numeric(obj)              # TRUE for post-LASSO
  beta     <- if (is_vec) obj else coef(obj)
  beta     <- beta[!is.na(beta)]
  terms    <- names(beta)
  
  # 4A) λ_max at (0,0,0)
  b1 <- beta["x_lag1"] %||% 0
  b2 <- beta["x_lag2"] %||% 0
  b3 <- beta["X_lag1"] %||% 0
  J0 <- matrix(c(
    b1,           b2,           b3,
    1,            0,            0,
    delta*b1, delta*b2, (1-delta)+delta*b3
  ), 3, 3, byrow = TRUE)
  
  # extract full spectrum
  ev0 <- eig(J0)
  cat("Eigenvalues at x = 0 :", round(ev0, 4), "\n")
  cat("  → max |λ|         :", round(max(Mod(ev0)), 4), "\n")
  cat("  → Im(λ)           :", round(Im(ev0), 4), "\n\n")
  
  
  # purely linear?
  if (!any(grepl("(_sq|_cub|I\\()", terms))){
    cat("Linear specification ⇒ unique SS at zero.\n")
    next
  }
  
  # 4B) Define scalar-safe F(x)
  # 4B vector-safe F(x) and G(x)=F(x)-x
  col_or_zero <- function(df, nm, n) {
    if (nm=="(Intercept)") {
      rep(1, n)
    } else if (nm %in% names(df)) {
      df[[nm]]
    } else {
      rep(0, n)
    }
  }
  
  # 4B) vector-safe F(x)  (always returns length(x))
  F <- function(x_vec) {
    sapply(x_vec, function(xx) {
      row <- make_row(xx)
      
      # start with intercept, if present
      h <- if ("(Intercept)" %in% names(beta)) beta["(Intercept)"] else 0
      
      # add every other regressor actually in 'beta'
      for (nm in setdiff(names(beta), "(Intercept)")) {
        h <- h + beta[nm] * row[[nm]]
      }
      h
    })
  }
  G <- function(x) F(x) - x                      # same length as x_vec
  
  # same vectorisation as F
  
  # 4C) ROBUST root search on [lo, hi] 
  # 4C)  analytic roots of G(x) 
  #
  #  G(x) = c0 + c1 x + c2 x^2 + c3 x^3      (degree ≤ 3 by construction)
  #  Evaluate G at 4 distinct points and solve for c0…c3, then polyroot().
  #
  xsamp <- c(0, 1, 2, 3)
  Gvals <- G(xsamp)                       # length-4 vector
  V     <- cbind(1, xsamp, xsamp^2, xsamp^3)   # Vandermonde 4×4
  coeff <- solve(V, Gvals)                # c0, c1, c2, c3
  
  roots_all <- polyroot(coeff)            # complex roots
  roots     <- Re(roots_all[abs(Im(roots_all)) < 1e-8])   # real ≈ roots
  
  # keep real roots inside [lo, hi] and not numerically 0
  roots <- unique(round(roots[roots >= lo & roots <= hi & abs(roots) > 1e-12], 10))
  
  if (!length(roots)) {
    cat("No non-zero SS in [", round(lo,2), ", ", round(hi,2), "].\n", sep="")
    next
  }
  cat("Non-zero SS :", paste(format(roots, digits=6), collapse = ", "), "\n")
  
  
  # 4D)  Jacobian & λ_max at each root (unchanged)
  # 4D)  Jacobian & λ_max at each non-zero SS  
  #
  #  Need a version of F that accepts three separate state variables
  #  (x_{t-1}, x_{t-2}, X_{t-1}).  We build the regression row explicitly.
  #
  F3 <- function(x1, x2, X1) {
    xl1 <- x1; xl2 <- x2; XL <- X1
    xl1_sq <- xl1^2; xl2_sq <- xl2^2; XL_sq <- XL^2
    xl1_cu <- xl1^3; xl2_cu <- xl2^3; XL_cu <- XL^3
    
    row_vals <- c(
      `(Intercept)`         = 1,
      x_lag1                = xl1,
      x_lag2                = xl2,
      X_lag1                = XL,
      x_lag1_sq             = xl1_sq,
      x_lag2_sq             = xl2_sq,
      X_lag1_sq             = XL_sq,
      x_lag1_cub            = xl1_cu,
      x_lag2_cub            = xl2_cu,
      X_lag1_cub            = XL_cu,
      x_lag1_x_lag2         = xl1 * xl2,
      x_lag1_X_lag1         = xl1 * XL,
      x_lag2_X_lag1         = xl2 * XL,
      x_lag1_sq_x_lag2      = xl1_sq * xl2,
      x_lag1_x_lag2_sq      = xl1 * xl2_sq,
      x_lag1_sq_X_lag1      = xl1_sq * XL,   # canonical
      x_lag2_sq_X_lag1      = xl2_sq * XL,
      x_lag1_X_lag1_sq      = xl1 * XL_sq,   # canonical
      x_lag2_X_lag1_sq      = xl2 * XL_sq,
      x_lag1_x_lag2_X_lag1  = xl1 * xl2 * XL
    )
    sum(beta[names(beta)] * row_vals[names(beta)])
  }
  
  H <- function(z) {
    x1 <- z[1]; x2 <- z[2]; X1 <- z[3]
    x_next <- F3(x1, x2, X1)
    X_next <- (1 - delta) * X1 + delta * x_next
    c(x_next, x1, X_next)
  }
  
  for (xs in roots) {
    JJ  <- jacobian(H, c(xs, xs, xs))
    for (xs in roots) {
      JJ  <- jacobian(H, c(xs, xs, xs))
      
      # full eigen‑spectrum
      evs <- eig(JJ)
      cat("x* =", format(xs, digits = 6), "\n")
      cat("  Eigenvalues:       ", round(evs, 6), "\n")
      cat("  → max |λ|:         ", round(max(Mod(evs)), 6), "\n")
      cat("  → Im(λ):           ", round(Im(evs), 6), "\n\n")
    }
  }
}

# Gather the zero‐root λ_max for each spec into a named vector
obs_lambda_zero <- sapply(names(specs), function(nm) {
  obj    <- specs[[nm]]
  is_vec <- is.numeric(obj)
  beta   <- if (is_vec) obj else coef(obj)
  b1     <- beta["x_lag1"] %||% 0
  b2     <- beta["x_lag2"] %||% 0
  b3     <- beta["X_lag1"] %||% 0
  J0     <- matrix(c(
    b1,           b2,           b3,
    1,            0,            0,
    delta*b1, delta*b2, (1-delta)+delta*b3
  ), 3, 3, byrow=TRUE)
  ρ(J0)
})
# Now obs_lambda_zero["Minimal"], ["Intermediate"], ["Full"] are observed values.


# Extended Table: Panel A (Fit) + Panel B (Local Stability)
panelA <- results_df

# Compute imaginary parts of the eigenvalues at the zero SS for each spec
# we already have ρ() and eig() defined above
obs_lambda_zero_Im <- sapply(names(specs), function(nm) {
  obj    <- specs[[nm]]
  is_vec <- is.numeric(obj)
  beta   <- if (is_vec) obj else coef(obj)
  b1     <- beta["x_lag1"] %||% 0
  b2     <- beta["x_lag2"] %||% 0
  b3     <- beta["X_lag1"] %||% 0
  J0     <- matrix(c(
    b1,           b2,           b3,
    1,            0,            0,
    delta*b1, delta*b2, (1-delta)+delta*b3
  ), 3, 3, byrow=TRUE)
  ev0    <- eig(J0)
 max(abs(Im(ev0)))
})

# Panel B: max eigenvalue at x*=0 from local‐stability block & its imaginary‐part
panelB <- data.frame(
  Model               = names(obs_lambda_zero),
  Lambda_max_at_0     = as.numeric(obs_lambda_zero),
  Imag_part_max_at_0  = as.numeric(obs_lambda_zero_Im),
  row.names           = NULL,
  stringsAsFactors    = FALSE
)

# (Optional) merge into one “extended” table 
extended3 <- merge(panelA, panelB, by="Model")

# Print them side by side 

cat("\n**Panel A – Fit Statistics**\n")
kable(panelA,     digits=3, align="lccc")

cat("\n**Panel B – Local Stability (λmax & Im at h* = 0)**\n")
kable(panelB, digits = 3, align = "lccc")

cat("\n**Extended Table – Combined**\n")
kable(extended3, digits = 3, align = "lccccc")




# STAR local stability: steady states and spectra 
delta_used <- delta 

star_ls <- lapply(star_models, function(m) star_steady_states(m, delta=delta_used))
# Summarize into a tidy data.frame (one row per (model, steady state))
panelB_star <- do.call(rbind, lapply(names(star_ls), function(nm) {
  items <- star_ls[[nm]]
  do.call(rbind, lapply(seq_along(items), function(i) {
    it <- items[[i]]
    data.frame(
      Model        = nm,
      SteadyState  = i,
      h_star       = round(it$h_star, 6),
      Lambda_max   = round(it$rho, 6),
      Imag_partMax = round(it$imag_max, 6),
      stringsAsFactors = FALSE
    )
  }))
}))
cat("\n**Panel B (STAR) – Local Stability (per steady state)**\n")
print(panelB_star)

# Pick the central SS (closest to zero) for each STAR model
pick_central <- function(df) df[which.min(abs(df$h_star)), , drop = FALSE]
panelB_star_central <- do.call(rbind, lapply(split(panelB_star, panelB_star$Model), pick_central))
# Rename/select to match panelB’s schema
panelB_star_central_to_merge <- panelB_star_central %>%
  dplyr::transmute(
    Model,
    Lambda_max_at_0    = Lambda_max,
    Imag_part_max_at_0 = Imag_partMax
  )
# Now stack
panelB_ext <- dplyr::bind_rows(panelB, panelB_star_central_to_merge)
# including STAR rows
# panelA 
extended3_ext <- merge(panelA, panelB_ext, by = "Model", all = TRUE)
# Print
cat("\n**Extended Table — with STAR (central SS)**\n")
kable(extended3_ext[order(extended3_ext$Model), ], digits = 3, align = "lccccc")

























################################################################################
# DETERMINISTIC FORECAST
################################################################################


# 0)  GLOBAL CONSTANTS
horizon <- 1000 * obs_per_year_real / 4      # simulate 1 000 quarters ≈ 250 y

# horizon_data = 1000

# Choose what to show
# Models that have (h, H) → eligible for the common phase portrait
phase_candidates <- c("Linear","Minimal","Intermediate","Full","LASSO", "STAR_std", "STAR_min", "ESTAR_min") # LASSO

# Models on the phase portrait:
include_phase <- setdiff(phase_candidates, c("Full"))  # <— example: drop LASSO

# bi-plot time-series: pick any labels from AR(2), Linear, Minimal, Intermediate, Full, LASSO
include_ts_a <- c("AR(2)","Minimal","Linear", "STAR_std")      # left panel
include_ts_b <- c("Intermediate","LASSO", "STAR_std", "STAR_min", "ESTAR_min")   # right panel (# Full, LASSO)

## whenn need a single‑step date increment
increment <- if (freq_real=="q") months(3) else
  if (freq_real=="m") months(1) else weeks(1)


# Safety Check AFTER
cat(obs_per_year_real, obs_unit_real, delta, N, horizon, "\n")



# 1)  Build initial h‑history, lags and starting H_{t-1}

# rename the arg for clarity
build_initial_objects <- function(df, origin_date) {
  stopifnot(origin_date %in% df$date)
  i0 <- which(df$date == origin_date)
  list(
    initial_state = c(
      h     = df$h[i0],
      h_lag = df$x_lag1[i0]
    ),
    H_start = df$X_lag1[i0]   # <- use precomputed X_{t-1}
  )
}

init_objs <- build_initial_objects(data_merged_real, origin_date = forecast_origin)
initial_h     <- init_objs$initial_state["h"]
initial_h_lag <- init_objs$initial_state["h_lag"]
H_tm1_start   <- init_objs$H_start



# 2)  Generic core that keeps the TRUE recursion for H_t

# PURE RECURSIVE core for simulations 
simulate_core_recursive <- function(h_tm1, h_tm2, H_tm1, horizon, delta, law_h) {
  h_out <- H_out <- numeric(horizon)
  for (t in 1:horizon) {
    h_t <- law_h(h_tm1, h_tm2, H_tm1)            # model's law of motion
    H_t <- (1 - delta) * H_tm1 + delta * h_t     # EWMA recursion (∞ lag)
    h_out[t] <- h_t; H_out[t] <- H_t
    h_tm2 <- h_tm1; h_tm1 <- h_t; H_tm1 <- H_t   # roll lags
  }
  data.frame(Time = 1:horizon, h = h_out, H = H_out)
}



# 3)  Model-specific wrappers

## 3.1  AR(2)  ––  no S_t term
simulate_forecast_ar2 <- function(model, h_tm1, h_tm2, horizon)
{
  co <- coef(model)
  b0 <- co["(Intercept)"]; b1 <- co["x_lag1"]; b2 <- co["x_lag2"]
  
  h_out <- numeric(horizon)
  for(t in 1:horizon) {
    h_t <- b0 + b1 * h_tm1 + b2 * h_tm2
    h_out[t] <- h_t
    h_tm2 <- h_tm1; h_tm1 <- h_t
  }
  data.frame(Time = 1:horizon, h = h_out)        # no H column
}

## 3.2  Linear
simulate_forecast_linear <- function(model, h0, h_1, H_1, horizon, delta)
{
  co <- coef(model)
  b0 <- co["(Intercept)"]; b1 <- co["x_lag1"]
  b2 <- co["x_lag2"];     b3 <- co["X_lag1"]
  
  law_h <- function(h1, h2, H1) b0 + b1 * h1 + b2 * h2 + b3 * H1
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, law_h)
}

## 3.3  Minimal
simulate_forecast_minimal <- function(model, h0, h_1, H_1, horizon, delta)
{
  co <- coef(model)
  b0 <- co["(Intercept)"]; b1 <- co["x_lag1"]; b2 <- co["x_lag2"]
  b3 <- co["X_lag1"];      b4 <- co["x_lag1_cub"]
  
  law_h <- function(h1, h2, H1) b0 + b1 * h1 + b2 * h2 + b3 * H1 + b4 * h1^3
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, law_h)
}

## 3.4  Intermediate
# small helper so missing coeffs default to zero
`%||%` <- function(a, b) if (is.na(a) || length(a) == 0) b else a

simulate_forecast_intermediate <- function(model, h0, h_1, H_1, horizon, delta) {
  co <- coef(model)
  b0 <- co["(Intercept)"] %||% 0
  b1 <- co["x_lag1"]      %||% 0
  b2 <- co["x_lag2"]      %||% 0
  b3 <- co["X_lag1"]      %||% 0
  b4 <- co["x_lag1_cub"]  %||% 0
  b5 <- co["x_lag2_cub"]  %||% 0
  b6 <- co["X_lag1_cub"]  %||% 0
  b7 <- co["x_lag1_sq_X_lag1"] %||% 0
  b8 <- co["x_lag1_X_lag1_sq"] %||% 0
  
  law_h <- function(h1, h2, H1)
    b0 + b1*h1 + b2*h2 + b3*H1 +
    b4*h1^3 + b5*h2^3 + b6*H1^3 +
    b7*(h1^2)*H1 + b8*h1*(H1^2)
  
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, law_h)
}

## 3.5  Full
simulate_forecast_full <- function(model, h0, h_1, H_1, horizon, delta) {
  co <- coef(model); z <- function(nm) co[nm] %||% 0
  b0  <- z("(Intercept)"); b1 <- z("x_lag1"); b2 <- z("x_lag2"); b3 <- z("X_lag1")
  b4  <- z("x_lag1_sq");   b5 <- z("x_lag2_sq");   b6  <- z("X_lag1_sq")
  b7  <- z("x_lag1_x_lag2"); b8 <- z("x_lag1_X_lag1"); b9 <- z("x_lag2_X_lag1")
  b10 <- z("x_lag1_cub");  b11 <- z("x_lag2_cub");  b12 <- z("X_lag1_cub")
  b13 <- z("x_lag1_sq_x_lag2"); b14 <- z("x_lag1_x_lag2_sq")
  b15 <- z("x_lag1_sq_X_lag1"); b16 <- z("x_lag2_sq_X_lag1")
  b17 <- z("x_lag1_X_lag1_sq"); b18 <- z("x_lag2_X_lag1_sq")
  b19 <- z("x_lag1_x_lag2_X_lag1")
  
  law_h <- function(h1, h2, H1)
    b0 + b1*h1 + b2*h2 + b3*H1 +
    b4*h1^2 + b5*h2^2 + b6*H1^2 +
    b7*h1*h2 + b8*h1*H1 + b9*h2*H1 +
    b10*h1^3 + b11*h2^3 + b12*H1^3 +
    b13*(h1^2)*h2 + b14*h1*(h2^2) +
    b15*(h1^2)*H1 + b16*(h2^2)*H1 +
    b17*h1*(H1^2) + b18*h2*(H1^2) +
    b19*h1*h2*H1
  
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, law_h)
}

## 3.6  LASSO (same structure as Minimal if cubic chosen)
simulate_forecast_lasso <- function(coef_vec, h0, h_1, H_1, horizon, delta) {
  # optional: drop NAs & ensure names are present
  coef_vec <- coef_vec[!is.na(coef_vec)]
  
  make_features <- function(h1, h2, H1) {
    df <- data.frame(
      x_lag1 = h1, x_lag2 = h2, X_lag1 = H1
    )
    df$x_lag1_sq            <- df$x_lag1^2
    df$x_lag2_sq            <- df$x_lag2^2
    df$X_lag1_sq            <- df$X_lag1^2
    df$x_lag1_x_lag2        <- df$x_lag1 * df$x_lag2
    df$x_lag1_X_lag1        <- df$x_lag1 * df$X_lag1
    df$x_lag2_X_lag1        <- df$x_lag2 * df$X_lag1
    df$x_lag1_cub           <- df$x_lag1^3
    df$x_lag2_cub           <- df$x_lag2^3
    df$X_lag1_cub           <- df$X_lag1^3
    df$x_lag1_sq_x_lag2     <- df$x_lag1^2 * df$x_lag2
    df$x_lag1_x_lag2_sq     <- df$x_lag1 * df$x_lag2^2
    df$x_lag1_sq_X_lag1     <- df$x_lag1^2 * df$X_lag1
    df$x_lag2_sq_X_lag1     <- df$x_lag2^2 * df$X_lag1
    df$x_lag1_X_lag1_sq     <- df$x_lag1 * df$X_lag1^2
    df$x_lag2_X_lag1_sq     <- df$x_lag2 * df$X_lag1^2
    df$x_lag1_x_lag2_X_lag1 <- df$x_lag1 * df$x_lag2 * df$X_lag1
    df
  }
  
  step_fun <- function(h1, h2, H1) {
    full_df <- make_features(h1, h2, H1)
    keep    <- intersect(names(coef_vec), c("(Intercept)", names(full_df)))
    mm      <- cbind("(Intercept)" = 1, full_df)[, keep, drop = FALSE]
    as.numeric(as.matrix(mm) %*% coef_vec[keep])
  }
  
  # >>> this was `law_h` before
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, step_fun)
}


## 3.x  STAR deterministic forecast 
simulate_forecast_star <- function(star_obj, h0, h_1, H_1, horizon, delta) {
  stopifnot(inherits(star_obj, "star_fit"))
  th   <- star_obj$coef
  spec <- star_obj$spec
  zvar <- star_obj$z_var
  v    <- star_obj$v
  
  # Per-step law of motion h_t = f(h_{t-1}, h_{t-2}, H_{t-1})
  law_h <- function(h1, h2, H1) {
    z <- if (zvar == "x_lag1") h1 else v[1]*h1 + v[2]*h2 + v[3]*H1
    if (spec == "standard_lstar") {
      a0 <- th["a0"]; a1 <- th["a1"]
      b10 <- th["b10"]; b11 <- th["b11"]
      b20 <- th["b20"]; b21 <- th["b21"]
      b30 <- th["b30"]; b31 <- th["b31"]
      eta <- th["eta"];  c  <- th["c"]
      G <- 1/(1+exp(-exp(eta)*(z-c)))
      (a0 + a1*G) + (b10 + b11*G)*h1 + (b20 + b21*G)*h2 + (b30 + b31*G)*H1
    } else {
      a0 <- th["a0"]
      b10 <- th["b10"]; b20 <- th["b20"]; b30 <- th["b30"]
      d1 <- th["d1"];   d2  <- th["d2"];   d3  <- th["d3"]
      eta <- th["eta"]; c   <- th["c"]
      G <- if (spec == "minimal_lstar") 1/(1+exp(-exp(eta)*(z-c))) else (1 - exp(-exp(eta)*(z-c)^2))
      a0 + (b10 + d1*G)*h1 + (b20 + d2*G)*h2 + (b30 + d3*G)*H1
    }
  }
  
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, law_h)
}




# 4)  Run every simulation
forecast_ar2 <- simulate_forecast_ar2(model_ar2,
                                      initial_h, initial_h_lag, horizon)

forecast_linear <- simulate_forecast_linear(model_linear,
                                            initial_h, initial_h_lag, H_tm1_start,
                                            horizon, delta)

forecast_minimal <- simulate_forecast_minimal(model_minimal,
                                              initial_h, initial_h_lag, H_tm1_start,
                                              horizon, delta)

forecast_inter <- simulate_forecast_intermediate(model_intermediate,
                                                 initial_h, initial_h_lag, H_tm1_start,
                                                 horizon, delta)

forecast_full <- simulate_forecast_full(model_full,
                                        initial_h, initial_h_lag, H_tm1_start,
                                        horizon, delta)

forecast_lasso <- simulate_forecast_lasso(coef_lasso_vec,
                                          initial_h, initial_h_lag, H_tm1_start,
                                          horizon, delta)

forecast_star_std  <- simulate_forecast_star(model_star_std,
                                             initial_h, initial_h_lag, H_tm1_start,
                                             horizon, delta)
forecast_star_min  <- simulate_forecast_star(model_star_min,
                                             initial_h, initial_h_lag, H_tm1_start,
                                             horizon, delta)
forecast_estar_min <- simulate_forecast_star(model_estar_min,
                                             initial_h, initial_h_lag, H_tm1_start,
                                             horizon, delta)

# sanity: first three rows of each
lapply(list(AR2=forecast_ar2, Linear=forecast_linear, Minimal=forecast_minimal,
            Intermediate=forecast_inter, Full=forecast_full, LASSO=forecast_lasso,
            Star_std = forecast_star_std, Star_min = forecast_star_min,
            Estar_min = forecast_estar_min),
       head, 3)


# Build sample from the same origin used to initialize the simulation
df_plot <- data_merged_real %>% arrange(date) %>% dplyr::filter(date >= forecast_origin)
if (nrow(df_plot) > horizon_data) df_plot <- df_plot[1:horizon_data, ]
real_df <- df_plot %>% dplyr::rename(x_data = h) %>% 
  dplyr::select(date, x_data) %>% dplyr::rename(value = x_data) %>%
  mutate(model = "Data")


# Now we can build dates safely
# requires: library(lubridate)
origin_date <- as.Date(forecast_origin)
# Use %m+% for month/quarter steps (rolls to last valid day), plain + for weeks
from <- if (freq_real %in% c("q", "m")) {
  origin_date %m+% increment
} else {
  origin_date + increment
}
forecast_dates <- seq(
  from       = from,               # first period after origin
  by         = obs_unit_real,      # "quarter"/"month"/"week"
  length.out = horizon_data
)


stopifnot(forecast_origin %in% data_merged_real$date)
row0 <- which(data_merged_real$date == forecast_origin)
stopifnot(is.finite(data_merged_real$x_lag1[row0]))

trim_and_stamp <- function(df, horizon_data, dates) {
  out <- df[1:horizon_data, , drop = FALSE]
  out$date <- dates
  out
}

forecast_ar2      <- trim_and_stamp(forecast_ar2,      horizon_data, forecast_dates)
forecast_linear   <- trim_and_stamp(forecast_linear,   horizon_data, forecast_dates)
forecast_minimal  <- trim_and_stamp(forecast_minimal,  horizon_data, forecast_dates)
forecast_inter    <- trim_and_stamp(forecast_inter,    horizon_data, forecast_dates)
forecast_full     <- trim_and_stamp(forecast_full,     horizon_data, forecast_dates)
forecast_lasso    <- trim_and_stamp(forecast_lasso,    horizon_data, forecast_dates)

forecast_star_std   <- trim_and_stamp(forecast_star_std,   horizon_data, forecast_dates)
forecast_star_min   <- trim_and_stamp(forecast_star_min,   horizon_data, forecast_dates)
forecast_estar_min  <- trim_and_stamp(forecast_estar_min,  horizon_data, forecast_dates)



##############################
# SERIES TIME PLOT AND COMPARISON : 3D (STATE SPACE)
##############################

# --- 4. Prepare Forecast Data for Plotting ---
add_state_space <- function(df_forecast) {
  df_forecast$h_lag <- c(NA, df_forecast$h[-nrow(df_forecast)])
  df_forecast <- df_forecast[-1, ]  # remove the first row with NA
  return(df_forecast)
}

# For AR(2) (which does not have accumulation), add h_lag:
forecast_ar2_state <- add_state_space(forecast_ar2)
# For models with accumulation (Linear, Minimal, Intermediate, Full, LASSO):
forecast_linear_state  <- add_state_space(forecast_linear)
forecast_minimal_state <- add_state_space(forecast_minimal)
forecast_inter_state   <- add_state_space(forecast_inter)
forecast_full_state    <- add_state_space(forecast_full)
forecast_lasso_state   <- add_state_space(forecast_lasso)

# print forecast
print(forecast_ar2_state)
print(forecast_linear_state)
print(forecast_minimal_state)
print(forecast_inter_state)
print(forecast_full_state)
print(forecast_lasso_state)
# --- 5. Create a Combined Data Frame for Overlay Plots ---
data_real <- df_samp %>% dplyr::rename(x_data = h)

# Create real_df: select date and x_data, rename x_data to value, and add a model identifier.
real_df <- data_real %>% 
  dplyr::select(date, x_data) %>% 
  dplyr::rename(value = x_data) %>% 
  mutate(model = "Data")

# Inspect the resulting data frame
head(real_df)



############################################################
# 3D + 2D State‐Space for All Models
############################################################
forecast_list <- list(
  Linear       = forecast_linear ,
  Minimal      = forecast_minimal ,
  Intermediate = forecast_inter ,
  Full         = forecast_full ,
  LASSO        = forecast_lasso ,
  STAR   = forecast_star_std ,
  STAR_min   = forecast_star_min ,
  ESTAR  = forecast_estar_min
)

for(model_name in names(forecast_list)) {
  df <- forecast_list[[model_name]]
  
  # extract vectors
  h <- df$h
  H <- df$H
  n <- length(h)
  
  # 1) Build the 3‐dimensional state‐space data frame
  state3d <- data.frame(
    h_tm2 = h[   1:(n-2) ],   # h at t-2
    H_tm1 = H[   2:(n-1) ],   # H at t-1
    h_tm1 = h[   2:(n-1) ]    # h at t-1
  )
  
  # 2) 3D Scatterplot (stars) 
  scatterplot3d(
    x     = state3d$h_tm2,
    y     = state3d$H_tm1,
    z     = state3d$h_tm1,
    pch   = 8,           # star symbol
    color = "black",
    angle = 20,          # front-on view
    main  = paste0("State Space (", model_name, " Model)"),
    xlab  = expression(h[t-2]),
    ylab  = expression(H[t-1]),
    zlab  = expression(h[t-1])
  )
  
  # 3) Build the 2‐dimensional projection data frame 
  proj2d <- data.frame(
    h_t = h,
    H_t = H
  )
  
  # 4) 2D path + stars 
  p <- ggplot(proj2d, aes(x = h_t, y = H_t)) +
    geom_path(size = 0.8) +
    geom_point(shape = 8, size = 1) +
    labs(
      title = paste0("Projection onto (h_t, H_t) - ", model_name),
      x     = expression(h[t]),
      y     = expression(H[t])
    ) +
    theme_minimal()
  
  print(p)
  
  # pause for a moment for one plot at a time
  readline(prompt = "Press <Enter> for next model…")
}







###### SUPERPOSITION – split into (A) baseline models vs (B) STAR models ######
# Which models to show in each panel
base_models <- c("Linear","Minimal","Intermediate","Full","LASSO")
star_models <- c("STAR","STAR_min","ESTAR")

# Only keep the ones that actually exist in forecast_list
base_present <- intersect(base_models, names(forecast_list))
star_present <- intersect(star_models, names(forecast_list))

# Helper to bind and rename
make_combined <- function(lst, keep) {
  bind_rows(lst[keep], .id = "model") %>%
    transmute(model, h_t = h, H_t = H) %>%
    drop_na(h_t, H_t)
}

combined_proj_base <- make_combined(forecast_list, base_present)
combined_proj_star <- make_combined(forecast_list, star_present)

model_cols <- c(
  Data         = "black",
  "AR(2)"      = "#1f78b4",
  Minimal      = "#33a02c",
  Linear       = "#fb9a99",
  Intermediate = "#e31a1c",
  Full         = "#ff7f00",
  LASSO        = "#6a3d9a",
  STAR         = "#008b8b",
  STAR_min     = "#b15928",
  ESTAR        = "#a6cee3"
)

# drop the old pal_base / pal_star and use ONE palette everywhere 
pal_models <- scale_color_manual(name = "Model", values = model_cols)

# (optional but robust): lock factor levels so ggplot won’t reshuffle legends
combined_proj_base$model <- factor(combined_proj_base$model, levels = names(model_cols))
combined_proj_star$model <- factor(combined_proj_star$model, levels = names(model_cols))

# (A) Baseline models — phase plot with same colors
p_base <- ggplot(combined_proj_base, aes(x = h_t, y = H_t, color = model)) +
  geom_path(linewidth = 1) +
  geom_point(shape = 8, size = 1, alpha = 0.6) +
  pal_models +                                    # <- same colors as time-series
  labs(title = "Superposed 2D Trajectories – Baseline Models",
       x = expression(h[t]), y = expression(H[t])) +
  theme_minimal() +
  theme(legend.position = "right")

# (B) STAR models — phase plot with same colors
p_star <- ggplot(combined_proj_star, aes(x = h_t, y = H_t, color = model)) +
  geom_path(linewidth = 1) +
  geom_point(shape = 8, size = 1, alpha = 0.6) +
  pal_models +                                    # <- same colors as time-series
  labs(title = "Superposed 2D Trajectories – STAR Models",
       x = expression(h[t]), y = expression(H[t])) +
  theme_minimal() +
  theme(legend.position = "right")

# Show side-by-side
gridExtra::grid.arrange(p_base, p_star, ncol = 2)




##############################
# SERIES TIME PLOT AND COMPARISON : 2D 
##############################

# For forecasts, we extract the forecasted cyclical variable from column 's' # (Note: AR(2) only has s; the others have h and H. Here we plot the cyclical variable x_t = s.) 
forecast_ar2_long <- forecast_ar2 %>% dplyr::select(date, h) %>% dplyr::rename(value = h) %>% mutate(model = "AR(2)") 
forecast_linear_long <- forecast_linear %>% dplyr::select(date, h) %>% dplyr::rename(value = h) %>% mutate(model = "Linear") 
forecast_minimal_long <- forecast_minimal %>% dplyr::select(date, h) %>% dplyr::rename(value = h) %>% mutate(model = "Minimal") 
forecast_inter_long <- forecast_inter %>% dplyr::select(date, h) %>% dplyr::rename(value = h) %>% mutate(model = "Intermediate") 
forecast_full_long <- forecast_full %>% dplyr::select(date, h) %>% dplyr::rename(value = h) %>% mutate(model = "Full") 
forecast_lasso_long <- forecast_lasso %>% dplyr::select(date, h) %>% dplyr::rename(value = h) %>% mutate(model = "LASSO")

forecast_star_std_long <- forecast_star_std %>% dplyr::select(date, h) %>% dplyr::rename(value = h) %>% mutate(model = "STAR")
forecast_star_min_long <- forecast_star_min %>% dplyr::select(date, h) %>% dplyr::rename(value = h) %>% mutate(model = "STAR_min")
forecast_estar_min_long <- forecast_estar_min %>% dplyr::select(date, h) %>% dplyr::rename(value = h) %>% mutate(model = "ESTAR")


# 1) Colours (make STAR labels consistent with hyphens)
model_cols <- c(
  Data         = "black",
  "AR(2)"      = "#1f78b4",
  Minimal      = "#33a02c",
  Linear       = "#fb9a99",
  Intermediate = "#e31a1c",
  Full         = "#ff7f00",
  LASSO        = "#6a3d9a",
  STAR   = "#008b8b",
  STAR_min   = "#b15928",
  ESTAR  = "#a6cee3"
)

# 2) Regular two plots (baseline models only) 
panel_a <- ggplot(
  bind_rows(real_df, forecast_ar2_long, forecast_minimal_long, forecast_linear_long),
  aes(x = date, y = value, color = model)
) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = model_cols) +
  labs(title = "Panel (a): Data, AR(2), Minimal & Linear",
       x = "Date", y = "Cyclical Variable", color = NULL) +
  theme_minimal()

panel_b <- ggplot(
  bind_rows(real_df, forecast_inter_long, forecast_full_long, forecast_lasso_long),
  aes(x = date, y = value, color = model)
) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = model_cols) +
  labs(title = "Panel (b): Data, Intermediate, Full & LASSO",
       x = "Date", y = "Cyclical Variable", color = NULL) +
  theme_minimal()

gridExtra::grid.arrange(panel_a, panel_b, ncol = 2)

# 3) One plot with STAR models only (plus Data) 
star_only_long <- bind_rows(
  real_df,
  forecast_star_std_long,
  forecast_star_min_long,
  forecast_estar_min_long
)

p_star <- ggplot(star_only_long, aes(x = date, y = value, color = model)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = model_cols[c("Data","STAR","STAR_min","ESTAR")]) +
  labs(title = "Panel (c): STAR Family",
       x = "Date", y = "Cyclical Variable", color = NULL) +
  theme_minimal()

print(p_star)





 
 
 
 
 




 
 
 
################################################################################
###### TRAJECTORY CLASSIFIER ######
# LC (persistent), Focus (damped), Direct (no persistent),
################################################################################

 # CONFIG 
 BURN_IN_MULT        <- 20    # burn-in = 20 * obs_per_year_real
 MIN_PEAKS_TAIL      <- 6     # min tail maxima to consider oscillation
 PERIOD_CV_MAX       <- 0.12  # tail period CV gate
 SPECTRAL_SHARE_MIN  <- 0.12  # share around spectral peak (±3 bins)
 SPECTRAL_FLAT_MAX   <- 0.90  # spectral flatness bound
 AMP_MIN_TAIL        <- 0.02  # min tail amplitude (scale of h)
 LC_DAMP_ABS_MAX     <- 0.010 # |median log amplitude ratio| for LC
 FOCUS_DAMP_MIN      <- 0.020 # median decay per cycle for Focus
 FOCUS_SIGN_PVAL     <- 0.05  # Wilcoxon one-sided p-value (negative)
 EXP_BOUND           <- 100  # 100 or 1e6 # explosive cutoff for |h_t|
 
 `%||%` <- function(a, b) if (is.null(a) || length(a)==0L || isTRUE(is.na(a))) b else a
 
 # HELPERS
 finite_ok <- function(x) all(is.finite(x)) && max(abs(x)) <= EXP_BOUND
 
 get_extrema <- function(x) {
   dx <- diff(x); s <- sign(dx); s[s==0] <- NA
   i_max <- which(head(s,-1) > 0 & tail(s,-1) < 0) + 1L
   i_min <- which(head(s,-1) < 0 & tail(s,-1) > 0) + 1L
   list(max = i_max, min = i_min)
 }
 
 spectral_peak_frequency <- function(x, spans=c(5,5), min_freq=1/400) {
   sp <- spec.pgram(x, taper=0.1, fast=TRUE, detrend=TRUE, plot=FALSE, spans=spans)
   f  <- as.numeric(sp$freq); ps <- as.numeric(sp$spec)
   ok <- is.finite(f) & is.finite(ps) & f >= min_freq & f <= 0.5
   if (!any(ok)) return(NA_real_)
   f[ok][ which.max(ps[ok]) ]
 }
 
 spectral_gate <- function(x_tail, band_bins=3,
                           share_min=SPECTRAL_SHARE_MIN, sfm_max=SPECTRAL_FLAT_MAX) {
   sp <- spec.pgram(x_tail, taper=0.1, fast=TRUE, detrend=TRUE, plot=FALSE, spans=c(5,5))
   ps <- as.numeric(sp$spec); ps <- ps[is.finite(ps) & ps > 0]
   if (length(ps) < 10L) return(FALSE)
   k  <- which.max(ps)
   lo <- max(1, k - band_bins); hi <- min(length(ps), k + band_bins)
   share_band <- sum(ps[lo:hi]) / sum(ps)
   sfm <- exp(mean(log(ps))) / mean(ps)
   (share_band >= share_min) && (sfm <= sfm_max)
 }
 
 bandpass_filter <- function(x, f_lo, f_hi, order=2) {
   if (!requireNamespace("signal", quietly=TRUE)) return(x)
   f_lo <- max(1e-6, min(f_lo, 0.49)); f_hi <- max(f_lo*1.05, min(f_hi, 0.5))
   bf <- signal::butter(order, c(f_lo,f_hi)/0.5, type="pass")
   as.numeric(signal::filtfilt(bf, x))
 }
 
 hilbert_envelope <- function(x) {
   if (requireNamespace("pracma", quietly=TRUE) &&
       "hilbert" %in% getNamespaceExports("pracma")) {
     Mod(pracma::hilbert(x))
   } else {
     N <- length(x); Xf <- fft(x); h <- numeric(N)
     if (N %% 2 == 0) { h[c(1,N/2+1)] <- 1; h[2:(N/2)] <- 2
     } else { h[1] <- 1; h[2:((N+1)/2)] <- 2 }
     Mod(fft(Xf*h, inverse=TRUE) / N)
   }
 }
 
 # MODEL ABSTRACTION
 is_star_fit <- function(x) inherits(x, "star_fit") ||
   (is.list(x) && all(c("spec","z_var","coef") %in% names(x)))
 
 .G_logistic <- function(z, eta, c) 1/(1 + exp(-exp(eta) * (z - c)))
 .G_estar    <- function(z, eta, c) 1 - exp(-exp(eta) * (z - c)^2)
 
 # Build F(h_{t-1}, h_{t-2}, H_{t-1}) for lm or STAR
 .build_F3 <- function(fit) {
   if (!is_star_fit(fit)) {
     # polynomial/linear lm
     co <- coef(fit); z <- function(nm) ifelse(is.na(co[nm]), 0, co[nm])
     function(x1,x2,X1) {
       b0  <- z("(Intercept)")
       b1  <- z("x_lag1");     b2  <- z("x_lag2");     b3  <- z("X_lag1")
       b4  <- z("x_lag1_sq");  b5  <- z("x_lag2_sq");  b6  <- z("X_lag1_sq")
       b7  <- z("x_lag1_x_lag2"); b8 <- z("x_lag1_X_lag1"); b9 <- z("x_lag2_X_lag1")
       b10 <- z("x_lag1_cub"); b11 <- z("x_lag2_cub"); b12 <- z("X_lag1_cub")
       b13 <- z("x_lag1_sq_x_lag2"); b14 <- z("x_lag1_x_lag2_sq")
       b15 <- z("x_lag1_sq_X_lag1"); b16 <- z("x_lag2_sq_X_lag1")
       b17 <- z("x_lag1_X_lag1_sq"); b18 <- z("x_lag2_X_lag1_sq")
       b19 <- z("x_lag1_x_lag2_X_lag1")
       b0 + b1*x1 + b2*x2 + b3*X1 +
         b4*x1^2 + b5*x2^2 + b6*X1^2 +
         b7*x1*x2 + b8*x1*X1 + b9*x2*X1 +
         b10*x1^3 + b11*x2^3 + b12*X1^3 +
         b13*x1^2*x2 + b14*x1*x2^2 +
         b15*x1^2*X1 + b16*x2^2*X1 +
         b17*x1*X1^2 + b18*x2*X1^2 +
         b19*x1*x2*X1
     }
   } else {
     th   <- fit$coef
     spec <- fit$spec             # "standard_lstar" | "minimal_lstar" | "minimal_estar"
     zvar <- fit$z_var            # "x_lag1" | "combo"
     v    <- fit$v %||% c(1,0,0)  # length-3 if combo
     function(x1,x2,X1) {
       z <- if (identical(zvar, "x_lag1")) x1 else v[1]*x1 + v[2]*x2 + v[3]*X1
       if (spec == "standard_lstar") {
         a0 <- th["a0"]; a1 <- th["a1"]
         b10 <- th["b10"]; b11 <- th["b11"]
         b20 <- th["b20"]; b21 <- th["b21"]
         b30 <- th["b30"]; b31 <- th["b31"]
         eta <- th["eta"];  c  <- th["c"]
         G <- .G_logistic(z, eta, c)
         (a0 + a1*G) + (b10 + b11*G)*x1 + (b20 + b21*G)*x2 + (b30 + b31*G)*X1
       } else {
         # minimal LSTAR or minimal ESTAR
         a0 <- th["a0"]
         b10 <- th["b10"]; b20 <- th["b20"]; b30 <- th["b30"]
         d1  <- th["d1"];  d2  <- th["d2"];  d3  <- th["d3"]
         eta <- th["eta"]; c   <- th["c"]
         G <- if (spec == "minimal_lstar") .G_logistic(z, eta, c) else .G_estar(z, eta, c)
         a0 + (b10 + d1*G)*x1 + (b20 + d2*G)*x2 + (b30 + d3*G)*X1
       }
     }
   }
 }
 
 find_fixed_points <- function(fit, lo, hi) {
   F3 <- .build_F3(fit); g <- function(s) F3(s,s,s) - s
   grid <- seq(lo, hi, length.out=601)
   gv <- sapply(grid, g); ok <- which(is.finite(gv))
   grid <- grid[ok]; gv <- gv[ok]
   chg  <- which(diff(sign(gv)) != 0)
   roots <- c()
   for (k in chg) {
     a <- grid[k]; b <- grid[k+1]
     r <- try(uniroot(g, c(a,b))$root, silent=TRUE)
     if (!inherits(r,"try-error") && is.finite(r)) roots <- c(roots, r)
   }
   unique(round(roots, 10))
 }
 
 jacobian_at_ss <- function(fit, delta, sstar) {
   if (!requireNamespace("numDeriv", quietly=TRUE))
     stop("Please install.packages('numDeriv')")
   F3 <- .build_F3(fit)
   H  <- function(z) {
     x1<-z[1]; x2<-z[2]; X1<-z[3]
     h <- F3(x1,x2,X1)
     c(h, x1, (1-delta)*X1 + delta*h)
   }
   numDeriv::jacobian(H, c(sstar, sstar, sstar))
 }
 
 linear_type_at_ss <- function(fit, delta, sstar) {
   ev <- eigen(jacobian_at_ss(fit, delta, sstar), only.values=TRUE)$values
   rho <- max(Mod(ev))
   if (rho >= 1) return(list(type="unstable", ev=ev, rho=rho))
   if (any(abs(Im(ev)) > 1e-10)) return(list(type="focus", ev=ev, rho=rho))
   list(type="node", ev=ev, rho=rho)
 }
 
 tail_diagnostics <- function(h, burn_in) {
   if (length(h) <= burn_in + 10L) return(list(ok=FALSE))
   x  <- tail(h, length(h) - burn_in)
   
   f0 <- spectral_peak_frequency(x)
   if (!is.finite(f0) || f0 <= 0) return(list(ok=FALSE))
   xb <- try(bandpass_filter(x, f0*0.85, min(0.5, f0*1.15)), silent=TRUE)
   if (inherits(xb,"try-error")) xb <- x
   
   ex <- get_extrema(xb); idx_max <- ex$max; nmax <- length(idx_max)
   if (nmax < MIN_PEAKS_TAIL) return(list(ok=FALSE))
   per <- diff(idx_max); if (length(per) < 2L) return(list(ok=FALSE))
   Pbar <- mean(per); cvP <- sd(per)/Pbar
   
   if (cvP > PERIOD_CV_MAX) return(list(ok=FALSE))
   if (!spectral_gate(xb))  return(list(ok=FALSE))
   if (diff(range(tail(xb, max(1L, round(Pbar))))) < AMP_MIN_TAIL) return(list(ok=FALSE))
   
   Amax   <- pmax(abs(xb[idx_max]), 1e-12)
   lograt <- diff(log(Amax))
   med_lr <- median(lograt); iqr_lr <- IQR(lograt)
   wt <- suppressWarnings(wilcox.test(lograt, alternative="less", mu=0))
   p_less <- if (is.list(wt) && is.finite(wt$p.value)) wt$p.value else 1
   
   list(ok=TRUE, Pbar=Pbar, cvP=cvP, med_lr=med_lr, iqr_lr=iqr_lr, p_less=p_less)
 }
 
 classify_dynamics <- function(path_df, fit, delta, obs_per_year_real,
                               burn_in = BURN_IN_MULT * obs_per_year_real,
                               linearise_fallback = TRUE) {
   h <- path_df$h
   if (!finite_ok(h)) return("Explosive")
   
   td <- tail_diagnostics(h, burn_in)
   
   if (isTRUE(td$ok)) {
     if (abs(td$med_lr) <= LC_DAMP_ABS_MAX && td$iqr_lr <= 3*LC_DAMP_ABS_MAX) return("LC")
     if (td$med_lr <= -FOCUS_DAMP_MIN && td$p_less < FOCUS_SIGN_PVAL)         return("Focus")
     return("LC")
   }
   
   # FAST MODE: skip the (slow) linearisation entirely if asked
   if (!linearise_fallback) return("Direct")
   
   # fallback: linearisation at nearby fixed point (slow)
   lo <- quantile(h, 0.01, na.rm=TRUE) - 2*sd(h, na.rm=TRUE)
   hi <- quantile(h, 0.99, na.rm=TRUE) + 2*sd(h, na.rm=TRUE)
   roots <- find_fixed_points(fit, lo, hi)
   if (length(roots)) {
     sstar <- roots[ which.min(abs(mean(tail(h, 50)) - roots)) ]
     lt <- linear_type_at_ss(fit, delta, sstar)
     if (lt$type == "focus") return("Focus")
     if (lt$type == "node")  return("Direct")
   }
   "Direct"
 }
 
 
 # SIMULATORS
 simulate_core_recursive <- function(h_tm1, h_tm2, H_tm1, horizon, delta, law_h) {
   h_out <- H_out <- numeric(horizon)
   for (t in 1:horizon) {
     h_t <- law_h(h_tm1, h_tm2, H_tm1)            # model's law of motion
     H_t <- (1 - delta) * H_tm1 + delta * h_t     # EWMA recursion (∞ lag)
     h_out[t] <- h_t; H_out[t] <- H_t
     h_tm2 <- h_tm1; h_tm1 <- h_t; H_tm1 <- H_t   # roll lags
   }
   data.frame(Time = 1:horizon, h = h_out, H = H_out)
 }
 
 simulate_forecast_ar2 <- function(model, h_tm1, h_tm2, horizon) {
   co <- coef(model); b0 <- co["(Intercept)"] %||% 0; b1 <- co["x_lag1"] %||% 0; b2 <- co["x_lag2"] %||% 0
   h_out <- numeric(horizon)
   for (t in 1:horizon) {
     h_t <- b0 + b1*h_tm1 + b2*h_tm2
     h_out[t] <- h_t; h_tm2 <- h_tm1; h_tm1 <- h_t
   }
   data.frame(h = h_out)
 }
 
 simulate_forecast_linear <- function(model, h0, h_1, H_1, horizon, delta) {
   co <- coef(model); b0 <- co["(Intercept)"] %||% 0; b1 <- co["x_lag1"] %||% 0
   b2 <- co["x_lag2"] %||% 0; b3 <- co["X_lag1"] %||% 0
   simulate_core_recursive(h0, h_1, H_1, horizon, delta,
                           function(h1,h2,H1) b0 + b1*h1 + b2*h2 + b3*H1)
 }
 
 simulate_forecast_minimal <- function(model, h0, h_1, H_1, horizon, delta) {
   co <- coef(model); z <- function(nm) co[nm] %||% 0
   simulate_core_recursive(h0, h_1, H_1, horizon, delta,
                           function(h1,h2,H1) z("(Intercept)") + z("x_lag1")*h1 + z("x_lag2")*h2 +
                             z("X_lag1")*H1 + z("x_lag1_cub")*h1^3)
 }
 
 simulate_forecast_intermediate <- function(model, h0, h_1, H_1, horizon, delta) {
   co <- coef(model); z <- function(nm) co[nm] %||% 0
   simulate_core_recursive(h0, h_1, H_1, horizon, delta,
                           function(h1,h2,H1)
                             z("(Intercept)") + z("x_lag1")*h1 + z("x_lag2")*h2 + z("X_lag1")*H1 +
                             z("x_lag1_cub")*h1^3 + z("x_lag2_cub")*h2^3 + z("X_lag1_cub")*H1^3 +
                             z("x_lag1_sq_X_lag1")*(h1^2)*H1 + z("x_lag1_X_lag1_sq")*h1*(H1^2))
 }
 
 simulate_forecast_full <- function(model, h0, h_1, H_1, horizon, delta) {
   co <- coef(model); z <- function(nm) co[nm] %||% 0
   simulate_core_recursive(h0, h_1, H_1, horizon, delta,
                           function(h1,h2,H1)
                             z("(Intercept)") + z("x_lag1")*h1 + z("x_lag2")*h2 + z("X_lag1")*H1 +
                             z("x_lag1_sq")*h1^2 + z("x_lag2_sq")*h2^2 + z("X_lag1_sq")*H1^2 +
                             z("x_lag1_x_lag2")*h1*h2 + z("x_lag1_X_lag1")*h1*H1 + z("x_lag2_X_lag1")*h2*H1 +
                             z("x_lag1_cub")*h1^3 + z("x_lag2_cub")*h2^3 + z("X_lag1_cub")*H1^3 +
                             z("x_lag1_sq_x_lag2")*(h1^2)*h2 + z("x_lag1_x_lag2_sq")*h1*(h2^2) +
                             z("x_lag1_sq_X_lag1")*(h1^2)*H1 + z("x_lag2_sq_X_lag1")*(h2^2)*H1 +
                             z("x_lag1_X_lag1_sq")*h1*(H1^2) + z("x_lag2_X_lag1_sq")*h2*(H1^2) +
                             z("x_lag1_x_lag2_X_lag1")*h1*h2*H1)
 }
 
 # Generic path simulator: if fit is STAR, use its F3; else route by name
 simulate_path <- function(mod_name, fit, h0, h_1, H_1, horizon, delta) {
   if (is_star_fit(fit) || grepl("STAR", mod_name, fixed=TRUE)) {
     F3 <- .build_F3(fit)
     return(simulate_core_recursive(h0, h_1, H_1, horizon, delta,
                                    function(h1,h2,H1) F3(h1,h2,H1)))
   }
   switch(mod_name,
          AR2          = simulate_forecast_ar2(fit, h0, h_1, horizon),
          Linear       = simulate_forecast_linear(fit, h0, h_1, H_1, horizon, delta),
          Minimal      = simulate_forecast_minimal(fit, h0, h_1, H_1, horizon, delta),
          Intermediate = simulate_forecast_intermediate(fit, h0, h_1, H_1, horizon, delta),
          Full         = simulate_forecast_full(fit, h0, h_1, H_1, horizon, delta),
          stop("Unknown model: ", mod_name))
 }
 
 # BOOTSTRAP
 .add_features_local <- function(df) {
   df$x_lag1_sq <- df$x_lag1^2; df$x_lag2_sq <- df$x_lag2^2; df$X_lag1_sq <- df$X_lag1^2
   df$x_lag1_cub <- df$x_lag1^3; df$x_lag2_cub <- df$x_lag2^3; df$X_lag1_cub <- df$X_lag1^3
   df$x_lag1_x_lag2 <- df$x_lag1*df$x_lag2; df$x_lag1_X_lag1 <- df$x_lag1*df$X_lag1
   df$x_lag2_X_lag1 <- df$x_lag2*df$X_lag1
   df$x_lag1_sq_x_lag2 <- df$x_lag1_sq*df$x_lag2; df$x_lag1_x_lag2_sq <- df$x_lag1*df$x_lag2_sq
   df$x_lag1_sq_X_lag1 <- df$x_lag1_sq*df$X_lag1; df$x_lag2_sq_X_lag1 <- df$x_lag2_sq*df$X_lag1
   df$x_lag1_X_lag1_sq <- df$x_lag1*df$X_lag1_sq; df$x_lag2_X_lag1_sq <- df$x_lag2*df$X_lag1_sq
   df$x_lag1_x_lag2_X_lag1 <- df$x_lag1*df$x_lag2*df$X_lag1
   df
 }
 add_features_safe <- function(df) {
   if (exists("add_features", mode="function")) add_features(df) else .add_features_local(df)
 }
 coef_or0 <- function(co, nm) ifelse(!is.na(co[nm]), co[nm], 0)
 
 # Mean-law bootstrap generator for lm and STAR (uses residual resampling)
 gen_boot_df <- function(model_name, fit, n, init_h, init_hlag, init_H, delta,
                         scheme=c("iid","wild")) {
   scheme <- match.arg(scheme)
   res0 <- if (is_star_fit(fit)) fit$resid else resid(fit)
   eps  <- if (scheme=="iid") sample(res0, n, TRUE) else sample(c(-1,1), n, TRUE)*res0
   
   # conditional mean law
   law <- if (is_star_fit(fit) || grepl("STAR", model_name, fixed=TRUE)) {
     F3 <- .build_F3(fit)
     function(h1,h2,H1) F3(h1,h2,H1)
   } else {
     co <- coef(fit)
     switch(model_name,
            AR2 = function(h1,h2,H1) coef_or0(co,"(Intercept)") + coef_or0(co,"x_lag1")*h1 + coef_or0(co,"x_lag2")*h2,
            Linear = function(h1,h2,H1) coef_or0(co,"(Intercept)") + coef_or0(co,"x_lag1")*h1 + coef_or0(co,"x_lag2")*h2 + coef_or0(co,"X_lag1")*H1,
            Minimal = function(h1,h2,H1) coef_or0(co,"(Intercept)") + coef_or0(co,"x_lag1")*h1 + coef_or0(co,"x_lag2")*h2 +
              coef_or0(co,"X_lag1")*H1 + coef_or0(co,"x_lag1_cub")*h1^3,
            Intermediate = function(h1,h2,H1) coef_or0(co,"(Intercept)") + coef_or0(co,"x_lag1")*h1 + coef_or0(co,"x_lag2")*h2 +
              coef_or0(co,"X_lag1")*H1 + coef_or0(co,"x_lag1_cub")*h1^3 + coef_or0(co,"x_lag2_cub")*h2^3 +
              coef_or0(co,"X_lag1_cub")*H1^3 + coef_or0(co,"x_lag1_sq_X_lag1")*(h1^2)*H1 +
              coef_or0(co,"x_lag1_X_lag1_sq")*h1*(H1^2),
            Full = function(h1,h2,H1) coef_or0(co,"(Intercept)") + coef_or0(co,"x_lag1")*h1 + coef_or0(co,"x_lag2")*h2 + coef_or0(co,"X_lag1")*H1 +
              coef_or0(co,"x_lag1_sq")*h1^2 + coef_or0(co,"x_lag2_sq")*h2^2 + coef_or0(co,"X_lag1_sq")*H1^2 +
              coef_or0(co,"x_lag1_x_lag2")*h1*h2 + coef_or0(co,"x_lag1_X_lag1")*h1*H1 + coef_or0(co,"x_lag2_X_lag1")*h2*H1 +
              coef_or0(co,"x_lag1_cub")*h1^3 + coef_or0(co,"x_lag2_cub")*h2^3 + coef_or0(co,"X_lag1_cub")*H1^3 +
              coef_or0(co,"x_lag1_sq_x_lag2")*(h1^2)*h2 + coef_or0(co,"x_lag1_x_lag2_sq")*h1*(h2^2) +
              coef_or0(co,"x_lag1_sq_X_lag1")*(h1^2)*H1 + coef_or0(co,"x_lag2_sq_X_lag1")*(h2^2)*H1 +
              coef_or0(co,"x_lag1_X_lag1_sq")*h1*(H1^2) + coef_or0(co,"x_lag2_X_lag1_sq")*h2*(H1^2) +
              coef_or0(co,"x_lag1_x_lag2_X_lag1")*h1*h2*H1
     )
   }
   
   h <- H <- numeric(n); h1 <- init_h; h2 <- init_hlag; H1 <- init_H
   for (t in 1:n) {
     mu <- law(h1,h2,H1); ht <- mu + eps[t]; Ht <- (1-delta)*H1 + delta*ht
     h[t] <- ht; H[t] <- Ht; h2 <- h1; h1 <- ht; H1 <- Ht
   }
   dfb <- data.frame(
     h      = h,
     x_lag1 = c(init_h, head(h, -1)),
     x_lag2 = c(init_hlag, init_h, head(h, -2)),
     X_lag1 = c(init_H, head(H, -1))
   )
   add_features_safe(dfb)
 }
 
 make_star_RV <- function(fit, df) {
   Rr <- make_star_R(fit)
   V  <- vcov_star(fit, df)     # uses numDeriv once per STAR spec
   list(R = Rr$R, r = Rr$r, V = V)
 }
 
 
 # Wald: lm via car::linearHypothesis; STAR via user-supplied list(R,r)
 wald_reject_5 <- function(fit, hyps) {
   if (is.null(hyps)) return(TRUE)
   if (!is_star_fit(fit)) {
     if (!requireNamespace("car", quietly=TRUE)) stop("Please install.packages('car')")
     out <- try(car::linearHypothesis(fit, hyps), silent=TRUE)
     if (inherits(out,"try-error")) return(FALSE)
     p <- suppressWarnings(out[["Pr(>F)"]][2]); is.finite(p) && p < 0.05
   } else {
     stopifnot(is.list(hyps), !is.null(hyps$R))
     th <- fit$coef
     nm <- names(th)
     R  <- hyps$R[, nm, drop = FALSE]
     r  <- hyps$r %||% rep(0, nrow(R))
     V  <- hyps$V
     if (is.null(V)) {
       # last-resort fallback (not recommended): try diag(se^2)
       se <- fit$se
       if (all(is.na(se))) return(TRUE)  # don’t block on missing SE
       V <- diag(se^2); dimnames(V) <- list(nm,nm)
     }
     V <- V[nm, nm, drop = FALSE]        # align
     diff <- as.numeric(R %*% th - r)
     RVRT <- try(solve(R %*% V %*% t(R)), silent = TRUE)
     if (inherits(RVRT, "try-error")) return(TRUE)
     W <- drop(t(diff) %*% RVRT %*% diff)
     p <- pchisq(W, df = nrow(R), lower.tail = FALSE)
     is.finite(p) && p < 0.05
   }
 }
 
 limit_cycle_frequency <- function(models, hyps_list, B, obs_per_year_real,
                                   init_h, init_hlag, init_H, horizon, delta,
                                   bootstrap_scheme="iid", verbose=TRUE) {
   mods <- names(models)
   cats <- c("LC","Focus","Direct")
   out_counts <- setNames(replicate(length(mods), integer(3), simplify=FALSE), mods)
   lc_detect_counts <- setNames(integer(length(mods)), mods)
   out_mat <- matrix(NA_character_, nrow=B, ncol=length(mods), dimnames=list(NULL, mods))
   diag <- setNames(vector("list", length(mods)), mods)
   
   nT <- nrow(df_samp)  # match sample length
   
   for (m in mods) {
     fit0 <- models[[m]]
     kept <- 0L; draws <- 0L; wald_pass <- 0L
     while (kept < B && draws < 10*B) {
       draws <- draws + 1L
       dfb <- gen_boot_df(m, fit0, nT, init_h, init_hlag, init_H, delta, scheme=bootstrap_scheme)
       
       # re-estimate on bootstrap sample
       fit_b <- if (!is_star_fit(fit0)) {
         try(lm(formula(fit0), data = dfb), silent = TRUE)
       } else {
         try(estimate_star(
           dfb,
           spec  = fit0$spec,
           z_var = fit0$z_var,
           v     = fit0$v %||% c(1,0,0),
           starts_gamma = c(1, 5),      # tighter grid -> faster
           starts_c     = c("median"),  # single center -> faster
           trace        = FALSE,
           maxit        = 800           # smaller budget for bootstrap
         ), silent = TRUE)
       }

       if (inherits(fit_b,"try-error")) next
       
       path  <- try(simulate_path(m, fit_b, init_h, init_hlag, init_H, horizon, delta), silent=TRUE)
       if (inherits(path,"try-error")) next
       cls   <- classify_dynamics(path, fit_b, delta, obs_per_year_real)
       if (identical(cls,"Explosive")) next
       
       kept <- kept + 1L; out_mat[kept, m] <- cls
       
       hyps <- hyps_list[[m]]
       wok  <- if (is.null(hyps)) TRUE else wald_reject_5(fit_b, hyps)
       if (isTRUE(wok)) wald_pass <- wald_pass + 1L
       
       if (cls %in% cats) out_counts[[m]][ match(cls,cats) ] <- out_counts[[m]][ match(cls,cats) ] + 1L
       if (identical(cls,"LC") && isTRUE(wok)) lc_detect_counts[[m]] <- lc_detect_counts[[m]] + 1L
     }
     if (verbose) {
       tab <- out_counts[[m]]
       message(sprintf("%-14s kept=%4d draws=%5d Wald-pass(kept)=%4d  LC=%3d  Focus=%3d  Direct=%3d  LC&Wald=%3d",
                       m, kept, draws, wald_pass, tab[1], tab[2], tab[3], lc_detect_counts[[m]]))
     }
     diag[[m]] <- list(draws=draws, kept=kept, wald_pass_kept=wald_pass,
                       counts_path=setNames(out_counts[[m]], cats),
                       lc_detect=lc_detect_counts[[m]])
   }
   
   freq_table_path <- do.call(rbind, lapply(names(out_counts), function(m) {
     tot <- sum(out_counts[[m]]); pr <- if (tot>0) 100*out_counts[[m]]/tot else c(0,0,0)
     data.frame(Model=m, `%-LC(path)`=round(pr[1],1), `%-Focus(path)`=round(pr[2],1), `%-Direct(path)`=round(pr[3],1))
   }))
   rownames(freq_table_path) <- NULL
   
   table4_like <- do.call(rbind, lapply(names(out_counts), function(m) {
     kept_m <- sum(out_counts[[m]]); lc_wald <- lc_detect_counts[[m]]
     data.frame(Model=m, `%-LC (LC & Wald)`= if (kept_m>0) round(100*lc_wald/kept_m,1) else 0)
   }))
   rownames(table4_like) <- NULL
   
   list(freq_table_path=freq_table_path, table4_like=table4_like, diag=diag, raw=out_mat)
 }
 

 ## B runner: B=1000 for base models, B=100 for STAR models
 run_lcf_mixedB <- function(models, hyps_list,
                            B_base = 1000, B_star = 100,
                            obs_per_year_real, init_h, init_hlag, init_H,
                            horizon, delta, forecast_dates = NULL, horizon_data = NULL,
                            bootstrap_scheme = "iid", verbose = TRUE) {
   
   base_mods <- c("AR2","Linear","Minimal","Intermediate","Full")
   star_mods <- c("STAR_std","STAR_min","ESTAR_min")
   
   models_base <- models[intersect(names(models), base_mods)]
   models_star <- models[intersect(names(models), star_mods)]
   
   hyps_base <- hyps_list[intersect(names(hyps_list), names(models_base))]
   hyps_star <- hyps_list[intersect(names(hyps_list), names(models_star))]
   
   # Run each block with its own B
   res_base <- limit_cycle_frequency(
     models_base, hyps_base, B = B_base,
     obs_per_year_real = obs_per_year_real,
     init_h = init_h, init_hlag = init_hlag, init_H = init_H,
     horizon = horizon, delta = delta,
     bootstrap_scheme = bootstrap_scheme, verbose = verbose
   )
   
   res_star <- limit_cycle_frequency(
     models_star, hyps_star, B = B_star,
     obs_per_year_real = obs_per_year_real,
     init_h = init_h, init_hlag = init_hlag, init_H = init_H,
     horizon = horizon, delta = delta,
     bootstrap_scheme = bootstrap_scheme, verbose = verbose
   )
   
   # Combine to a single result object compatible with make_table4_summary()
   list(
     freq_table_path = rbind(res_base$freq_table_path, res_star$freq_table_path),
     table4_like     = rbind(res_base$table4_like,     res_star$table4_like),
     diag            = c(res_base$diag,                res_star$diag),
     raw             = list(base = res_base$raw, star = res_star$raw)  # keep both raws
   )
 }
 
 
 # Simulations 
 deterministic_paths <- function(models, init_h, init_hlag, init_H, horizon, delta,
                                 forecast_dates=NULL, horizon_data=NULL) {
   det <- lapply(names(models), function(m) {
     fit <- models[[m]]
     path <- simulate_path(m, fit, init_h, init_hlag, init_H, horizon, delta)
     n <- if (is.null(horizon_data)) nrow(path) else min(nrow(path), horizon_data)
     out <- path[seq_len(n), , drop=FALSE]; out$model <- m; out$value <- out$h
     out$date <- if (!is.null(forecast_dates) && length(forecast_dates) >= n)
       forecast_dates[seq_len(n)] else seq_len(n)
     out[, c("date","model","value")]
   })
   do.call(rbind, det)
 }
 
 sample_bootstrap_paths <- function(models, hyps_list, K, class_filter=c("any","LC","Focus","Direct","LC_Wald"),
                                    obs_per_year_real, init_h, init_hlag, init_H, horizon, delta,
                                    forecast_dates=NULL, horizon_data=NULL,
                                    bootstrap_scheme=c("iid","wild"), max_draws_per_model=5000, verbose=TRUE) {
   class_filter <- match.arg(class_filter); bootstrap_scheme <- match.arg(bootstrap_scheme)
   use_wald <- identical(class_filter, "LC_Wald")
   mods <- names(models); res <- vector("list", length(mods)); names(res) <- mods
   nT <- nrow(df_samp)
   
   for (m in mods) {
     fit0 <- models[[m]]
     kept <- 0L; draws <- 0L; bucket <- list()
     while (kept < K && draws < max_draws_per_model) {
       draws <- draws + 1L
       dfb <- gen_boot_df(m, fit0, nT, init_h, init_hlag, init_H, delta, scheme=bootstrap_scheme)
       
       fit_b <- if (!is_star_fit(fit0)) {
         try(lm(formula(fit0), data = dfb), silent = TRUE)
       } else {
         try(estimate_star(
           dfb,
           spec  = fit0$spec,
           z_var = fit0$z_var,
           v     = fit0$v %||% c(1,0,0),
           starts_gamma = c(1, 5),      # tighter grid -> faster
           starts_c     = c("median"),  # single center -> faster
           trace        = FALSE,
           maxit        = 800           # smaller budget for bootstrap
         ), silent = TRUE)
       }
       
       if (inherits(fit_b,"try-error")) next
       
       path  <- try(simulate_path(m, fit_b, init_h, init_hlag, init_H, horizon, delta), silent=TRUE)
       if (inherits(path,"try-error")) next
       if (!finite_ok(path$h)) next
       cls <- classify_dynamics(path, fit_b, delta, obs_per_year_real, )
       
       wald_ok <- TRUE
       if (use_wald) { hyps <- hyps_list[[m]]; wald_ok <- wald_reject_5(fit_b, hyps) }
       
       cls_ok <- (class_filter=="any") ||
         (class_filter=="LC"      && cls=="LC") ||
         (class_filter=="Focus"   && cls=="Focus") ||
         (class_filter=="Direct"  && cls=="Direct") ||
         (class_filter=="LC_Wald" && cls=="LC" && isTRUE(wald_ok))
       if (!cls_ok) next
       
       kept <- kept + 1L
       n <- if (is.null(horizon_data)) nrow(path) else min(nrow(path), horizon_data)
       out <- path[seq_len(n), , drop=FALSE]
       out$model <- m; out$sim_id <- kept; out$value <- out$h
       out$date <- if (!is.null(forecast_dates) && length(forecast_dates) >= n)
         forecast_dates[seq_len(n)] else seq_len(n)
       bucket[[kept]] <- out[, c("date","model","sim_id","value")]
     }
     if (verbose) {
       msg <- if (class_filter=="LC_Wald") "LC & Wald" else class_filter
       message(sprintf("[%-12s] kept=%d (draws=%d) for class='%s'", m, kept, draws, msg))
     }
     res[[m]] <- if (length(bucket)) do.call(rbind, bucket) else NULL
   }
   do.call(rbind, Filter(Negate(is.null), res))
 }
 
 snapshot_plot <- function(models, hyps_list, class_filter, K,
                           obs_per_year_real, init_h, init_hlag, init_H,
                           horizon, delta, forecast_dates, horizon_data) {
   if (!requireNamespace("ggplot2", quietly=TRUE)) stop("Please install.packages('ggplot2')")
   det_df  <- deterministic_paths(models, init_h, init_hlag, init_H, horizon, delta,
                                  forecast_dates, horizon_data)
   boot_df <- sample_bootstrap_paths(models, hyps_list, K, class_filter,
                                     obs_per_year_real, init_h, init_hlag, init_H,
                                     horizon, delta, forecast_dates, horizon_data,
                                     bootstrap_scheme="iid", verbose=TRUE)
   
   ggplot2::ggplot() +
     ggplot2::geom_line(data = boot_df,
                        ggplot2::aes(x = date, y = value, group = interaction(model, sim_id)),
                        color = "grey70", alpha = 0.6, linewidth = 0.4) +
     ggplot2::geom_line(data = det_df,
                        ggplot2::aes(x = date, y = value, color = model),
                        linewidth = 0.9) +
     ggplot2::facet_wrap(~ model, scales = "free_y", ncol = 2) +
     ggplot2::labs(
       title = "Deterministic path (color) with selected bootstrap deterministic paths (grey)",
       subtitle = sprintf("Class filter: %s • K per model: %d • Origin: %s",
                          class_filter, K,
                          if (!is.null(forecast_dates) && length(forecast_dates)>=1) format(forecast_dates[1], "%Y-Q%q") else "t=1"),
       x = "Date", y = "Cyclical component h_t") +
     ggplot2::theme_minimal(base_size = 12) +
     ggplot2::theme(legend.position = "none",
                    panel.grid.minor = ggplot2::element_blank())
 }
 

 # table printed
 make_table4_summary <- function(res, title = NULL, digits = 1) {
   di   <- res$diag
   mods <- names(di)
   
   pct   <- function(x, d) if (d > 0) 100 * x / d else 0
   fmtcp <- function(n, p) sprintf("%d (%.1f%%)", as.integer(n), round(p, digits))
   
   rows <- lapply(mods, function(m) {
     d        <- di[[m]]
     kept     <- as.integer(d$kept)
     draws    <- as.integer(d$draws)
     waldpass <- as.integer(d$wald_pass_kept)
     lc       <- as.integer(d$counts_path["LC"])
     foc      <- as.integer(d$counts_path["Focus"])
     dirc     <- as.integer(d$counts_path["Direct"])
     lc_wald  <- as.integer(d$lc_detect)
     
     c(
       Model             = m,
       Kept              = as.character(kept),
       Draws             = as.character(draws),
       `Wald-pass (kept)`= fmtcp(waldpass, pct(waldpass, kept)),
       LC                = fmtcp(lc,       pct(lc,       kept)),
       Focus             = fmtcp(foc,      pct(foc,      kept)),
       Direct            = fmtcp(dirc,     pct(dirc,     kept)),
       `LC & Wald`       = fmtcp(lc_wald,  pct(lc_wald,  kept))
     )
   })
   
   tab <- do.call(rbind, rows)
   tab <- as.data.frame(tab, stringsAsFactors = FALSE, check.names = FALSE)
   
   if (!is.null(title)) {
     cat("\n", title, "\n", strrep("=", nchar(title)), "\n", sep = "")
   }
   
   print_table <- function(df) {
     cols <- names(df)
     widths <- sapply(seq_along(df), function(j) max(nchar(c(cols[j], df[[j]]))))
     sep_line <- paste0("+", paste0(sapply(widths, function(w) strrep("-", w + 2)), collapse = "+"), "+")
     cat(sep_line, "\n")
     cat("|", paste0(sprintf(paste0(" %-", widths, "s "), cols), collapse = "|"), "|\n", sep = "")
     cat(sep_line, "\n")
     for (i in seq_len(nrow(df))) {
       cat("|", paste0(sprintf(paste0(" %-", widths, "s "), df[i, ]), collapse = "|"), "|\n", sep = "")
     }
     cat(sep_line, "\n")
   }
   
   print_table(tab)
   invisible(tab)
 }
 


 make_star_R <- function(fit) {
   stopifnot(is_star_fit(fit))
   nm   <- names(fit$coef)
   need <- switch(fit$spec,
                  "standard_lstar" = c("a1","b11","b21","b31"),
                  "minimal_lstar"  = c("d1","d2","d3"),
                  "minimal_estar"  = c("d1","d2","d3"),
                  stop("Unknown STAR spec: ", fit$spec))
   have <- intersect(need, nm)
   if (!length(have))
     stop("None of the target coefficients found in STAR fit. Coefs are: ",
          paste(nm, collapse = ", "))
   R <- matrix(0, nrow = length(have), ncol = length(nm),
               dimnames = list(have, nm))
   for (i in seq_along(have)) R[i, have[i]] <- 1
   list(R = R, r = numeric(nrow(R)))
 }
 
  models <- list(
    AR2 = model_ar2, Linear = model_linear, Minimal = model_minimal,
    Intermediate = model_intermediate, Full = model_full,
    STAR_std   = model_star_std,    
    STAR_min   = model_star_min,      
    ESTAR_min  = model_estar_min    
  )
  hyps_list <- list(
    AR2 = NULL,
    Linear = NULL,
    Minimal = "x_lag1_cub = 0",
    Intermediate = c("x_lag1_cub = 0","x_lag2_cub = 0","X_lag1_cub = 0",
                     "x_lag1_sq_X_lag1 = 0","x_lag1_X_lag1_sq = 0"),
    Full = c("x_lag1_sq = 0","x_lag2_sq = 0","X_lag1_sq = 0",
             "x_lag1_x_lag2 = 0","x_lag1_X_lag1 = 0","x_lag2_X_lag1 = 0",
             "x_lag1_cub = 0","x_lag2_cub = 0","X_lag1_cub = 0",
             "x_lag1_sq_x_lag2 = 0","x_lag1_x_lag2_sq = 0",
             "x_lag1_sq_X_lag1 = 0","x_lag2_sq_X_lag1 = 0",
             "x_lag1_X_lag1_sq = 0","x_lag2_X_lag1_sq = 0",
             "x_lag1_x_lag2_X_lag1 = 0")
  )
  
  model_star_std  <- estimate_star(df_samp, spec="standard_lstar", z_var="x_lag1")
  model_star_min  <- estimate_star(df_samp, spec="minimal_lstar",  z_var="x_lag1")
  model_estar_min <- estimate_star(df_samp, spec="minimal_estar",  z_var="x_lag1")
  
  # hyps
  if ("STAR_std"  %in% names(models))  hyps_list[["STAR_std"]]  <- make_star_RV(models[["STAR_std"]],  df_samp)
  if ("STAR_min"  %in% names(models))  hyps_list[["STAR_min"]]  <- make_star_RV(models[["STAR_min"]],  df_samp)
  if ("ESTAR_min" %in% names(models))  hyps_list[["ESTAR_min"]] <- make_star_RV(models[["ESTAR_min"]], df_samp)
  
  
  
###### CALL
  res <- run_lcf_mixedB(
    models, hyps_list,
    B_base = 1000, B_star = 100,
    obs_per_year_real = obs_per_year_real,
    init_h = initial_h, init_hlag = initial_h_lag, init_H = H_tm1_start,
    horizon = horizon, delta = delta,
    bootstrap_scheme = "iid", verbose = TRUE
  )
  
  make_table4_summary(res, title = "Limit-cycle classification — mixed B (base=1000, STAR=100)")
  
  make_table4_summary(res)
 
 
 
 
  
  
  
  
  
  
  
  
  
 

 ################################################################################
 # PLOT OF BOOTSTRAPPED CLASSIFIED TREAJECTORIES
 ################################################################################
 
 if (!exists("finite_ok", mode = "function")) {
   EXP_BOUND_LOCAL <- get0("EXP_BOUND", ifnotfound = 1e6)
   finite_ok <- function(x) all(is.finite(x)) && max(abs(x)) <= EXP_BOUND_LOCAL
 }
 
 # deterministic (baseline) paths (explicit args) 
 deterministic_paths <- function(models,
                                 init_h, init_hlag, init_H,
                                 horizon, delta,
                                 forecast_dates = NULL,
                                 horizon_data   = NULL) {
   det <- lapply(names(models), function(m) {
     fit  <- models[[m]]
     path <- simulate_path(m, fit, init_h, init_hlag, init_H, horizon, delta)
     stopifnot(is.data.frame(path), "h" %in% names(path))
     n <- if (is.null(horizon_data)) nrow(path) else min(nrow(path), horizon_data)
     
     out <- path[seq_len(n), , drop = FALSE]
     out$model <- m
     out$value <- out$h
     
     # stamp dates only if available & matching
     out$date <- if (!is.null(forecast_dates) && length(forecast_dates) >= n) {
       forecast_dates[seq_len(n)]
     } else {
       seq_len(n)  # numeric fallback
     }
     
     out[, c("date","model","value")]
   })
   do.call(rbind, det)
 }
 
 # sample bootstrap paths by class (optional Wald for LC_Wald) 
  sample_bootstrap_paths <- function(models, hyps_list,
                                     K = 10,
                                     class_filter = c("any","LC","Focus","Direct","LC_Wald"),
                                     obs_per_year_real,
                                     init_h, init_hlag, init_H,
                                     horizon, delta,
                                     forecast_dates = NULL,
                                     horizon_data   = NULL,
                                     bootstrap_scheme = c("iid","wild"),
                                     max_draws_per_model = 300,
                                     verbose = TRUE) {
    class_filter     <- match.arg(class_filter)
    bootstrap_scheme <- match.arg(bootstrap_scheme)
    use_wald <- identical(class_filter, "LC_Wald")
    
    mods <- names(models)
    res  <- vector("list", length(mods)); names(res) <- mods
    nT   <- nrow(df_samp)   
    
    for (m in mods) {
      fit0 <- models[[m]]
      kept <- 0L; draws <- 0L
      bucket <- list()
      
      while (kept < K && draws < max_draws_per_model) {
        draws <- draws + 1L
        
        # 1) bootstrap pseudo-sample (works for lm and STAR)
        dfb <- gen_boot_df(m, fit0, nT, init_h, init_hlag, init_H, delta,
                           scheme = bootstrap_scheme)
        
        # 2) re-estimate on bootstrap sample
        fit_b <- if (!is_star_fit(fit0)) {
          try(lm(formula(fit0), data = dfb), silent = TRUE)
        } else {
          try(estimate_star(
            dfb,
            spec  = fit0$spec,
            z_var = fit0$z_var,
            v     = fit0$v %||% c(1,0,0),
            starts_gamma = c(1, 5),      # tighter grid -> faster
            starts_c     = c("median"),  # single center -> faster
            trace        = FALSE,
            maxit        = 800           # smaller budget for bootstrap
          ), silent = TRUE)
        }
        
        if (inherits(fit_b, "try-error")) next
        
        # 3) simulate deterministic skeleton; discard ONLY if explosive
        path <- try(simulate_path(m, fit_b, init_h, init_hlag, init_H, horizon, delta),
                    silent = TRUE)
        if (inherits(path, "try-error")) next
        if (!finite_ok(path$h)) next
        
        # 4) classify (works for lm and STAR)
        cls <- classify_dynamics(
          path_df = path,
          fit     = fit_b,
          delta   = delta,
          obs_per_year_real = obs_per_year_real,
          linearise_fallback = TRUE   # match limit_cycle_frequency()
        )
        
        
        # 5) optional Wald filter (ONLY when class_filter == "LC_Wald")
        #    For STAR, hyps_list[[m]] should be list(R=..., r=..., V=...) precomputed once
        #    (so we never need SEs here).
        wald_ok <- TRUE
        if (use_wald) {
          hyps <- hyps_list[[m]]
          wald_ok <- wald_reject_5(fit_b, hyps)
        }
        
        # 6) class filter
        keep_this <-
          (class_filter == "any") ||
          (class_filter == "LC"       && cls == "LC") ||
          (class_filter == "Focus"    && cls == "Focus") ||
          (class_filter == "Direct"   && cls == "Direct") ||
          (class_filter == "LC_Wald"  && cls == "LC" && isTRUE(wald_ok))
        if (!keep_this) next
        
        # 7) keep this path for plotting
        kept <- kept + 1L
        n <- if (is.null(horizon_data)) nrow(path) else min(nrow(path), horizon_data)
        out <- path[seq_len(n), , drop = FALSE]
        out$model  <- m
        out$sim_id <- kept
        out$value  <- out$h
        out$cls    <- cls
        
        out$date <- if (!is.null(forecast_dates) && length(forecast_dates) >= n) {
          forecast_dates[seq_len(n)]
        } else {
          seq_len(n)
        }
        
        bucket[[kept]] <- out[, c("date","model","sim_id","value","cls")]
      }
      
      if (verbose) {
        msg_cf <- if (class_filter == "LC_Wald") "LC & Wald" else class_filter
        message(sprintf("[%-12s] kept=%d (draws=%d) for class='%s'", m, kept, draws, msg_cf))
      }
      res[[m]] <- if (length(bucket)) do.call(rbind, bucket) else NULL
    }
    
    do.call(rbind, Filter(Negate(is.null), res))
  }
 
 # ready-to-run example (build data + plot) 
 # Choose what to show
 K_plot      <- 10
 class_filter <- "any"  # "any" | "LC" | "Focus" | "Direct" | "LC_Wald"
 
 # 1) deterministic colored paths
 det_df <- deterministic_paths(
   models,
   init_h         = initial_h,
   init_hlag      = initial_h_lag,
   init_H         = H_tm1_start,
   horizon        = horizon,
   delta          = delta,
   forecast_dates = get0("forecast_dates", ifnotfound = NULL),
   horizon_data   = get0("horizon_data",   ifnotfound = nrow(df_samp))
 )
 
 # 2) pick bootstrap class and sample paths
 boot_df <- sample_bootstrap_paths(
   models, hyps_list,
   K                  = K_plot,
   class_filter       = class_filter,
   obs_per_year_real  = obs_per_year_real,
   init_h             = initial_h,
   init_hlag          = initial_h_lag,
   init_H             = H_tm1_start,
   horizon            = horizon,
   delta              = delta,
   forecast_dates     = get0("forecast_dates", ifnotfound = NULL),
   horizon_data       = get0("horizon_data",   ifnotfound = nrow(df_samp)),
   bootstrap_scheme   = "iid",
   verbose            = TRUE
 )
 
 # 3) fallback if none matched the filter
 if (is.null(boot_df) || !nrow(boot_df)) {
   message("No bootstrap paths matched the filter; plotting deterministic paths only.")
   boot_df <- det_df; boot_df$sim_id <- 0L
 }
 
 # 4) subtitle helpers
 origin_str <- if (!is.null(get0("forecast_dates")) && length(forecast_dates) >= 1)
   format(forecast_dates[1], "%Y-Q%q") else "t=1"
 
 sub_str <- sprintf("Class filter: %s • K per model: %d • Origin: %s",
                    if (identical(class_filter, "LC_Wald")) "LC (Wald)" else class_filter,
                    K_plot, origin_str)
 
 # 5) plot
 ####### EACH SELECTED TRAJECTORY WITH A DISTINCTIVE COLOUR

 # 1) Trajectory class colours 
 class_cols <- c(
   LC     = "forestgreen",
   Focus  = "royalblue",
   Direct = "firebrick"
 )
 
 # 2) Reuse  existing model palette (from above)
 #    Add aliases so names match what snapshot code uses
 model_cols_snapshot <- c(
   model_cols,                       
   "AR2"        = model_cols["AR(2)"],
   "STAR_std"   = model_cols["STAR"],
   "STAR_min"   = model_cols["STAR_min"],
   "ESTAR_min"  = model_cols["ESTAR"]
 )
 model_cols_snapshot <- model_cols_snapshot[!is.na(model_cols_snapshot)]
 
 # 3) One unified palette for BOTH aesthetics
 pal_all <- c(class_cols, model_cols_snapshot)
 
 # 4) Plot
 ggplot2::ggplot() +
   ggplot2::geom_line(
     data = boot_df,
     ggplot2::aes(x = date, y = value,
                  group = interaction(model, sim_id),
                  color = cls),
     alpha = 0.3, linewidth = 0.5
   ) +
   ggplot2::geom_line(
     data = det_df,
     ggplot2::aes(x = date, y = value, color = model),
     linewidth = 0.9
   ) +
   ggplot2::facet_wrap(~ model, scales = "free_y", ncol = 2) +
   ggplot2::scale_color_manual(values = pal_all, guide = "none") +
   ggplot2::labs(
     title = "Deterministic path (bold) and bootstrap paths (transparent by class)",
     subtitle = sub_str, x = "Date", y = "Cyclical component h[t]"
   ) +
   ggplot2::theme_minimal(base_size = 12) +
   ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
 












 
 ################################################################
 ## CYCLE-LENGTH ESTIMATION & COMPARISON (EMPIRICAL vs. MODELS) 
 ################################################################
 `%||%` <- function(a,b) if (length(a)==0L || isTRUE(is.na(a))) b else a
 
 if (!exists("is_star_fit", mode="function")) {
   is_star_fit <- function(x) inherits(x, "star_fit") ||
     (is.list(x) && all(c("spec","z_var","coef") %in% names(x)))
 }
 
 # Small guard for finite_ok
 if (!exists("finite_ok", mode = "function")) {
   EXP_BOUND_LOCAL <- get0("EXP_BOUND", ifnotfound = 1e6)
   finite_ok <- function(x) all(is.finite(x)) && max(abs(x)) <= EXP_BOUND_LOCAL
 }
 
 # Helper: compute G_t for a given df and baseline STAR fit (uses eta,c,v)
 .compute_G <- function(dfb, fit0) {
   stopifnot(is.list(fit0), !is.null(fit0$spec), !is.null(fit0$coef))
   z <- if (identical(fit0$z_var, "x_lag1")) {
     dfb$x_lag1
   } else {
     v <- fit0$v %||% c(1,0,0)
     v[1]*dfb$x_lag1 + v[2]*dfb$x_lag2 + v[3]*dfb$X_lag1
   }
   eta <- unname(fit0$coef["eta"])
   c0  <- unname(fit0$coef["c"])
   if (fit0$spec == "standard_lstar" || fit0$spec == "minimal_lstar") {
     1 / (1 + exp(-exp(eta) * (z - c0)))
   } else { # minimal_estar
     1 - exp(-exp(eta) * (z - c0)^2)
   }
 }
 
 # FAST STAR RE-FIT: OLS-only update of linear parameters; keep (eta,c,v)
 fast_refit_star <- function(dfb, fit0) {
   stopifnot(is.list(fit0), !is.null(fit0$spec))
   G <- .compute_G(dfb, fit0)
   y <- dfb$h; x1 <- dfb$x_lag1; x2 <- dfb$x_lag2; X1 <- dfb$X_lag1
   
   coef_new <- fit0$coef
   
   if (fit0$spec == "standard_lstar") {
     # h = (a0 + a1*G) + (b10 + b11*G)x1 + (b20 + b21*G)x2 + (b30 + b31*G)X1
     X <- cbind(1, G, x1, G*x1, x2, G*x2, X1, G*X1)
     colnames(X) <- c("a0","a1","b10","b11","b20","b21","b30","b31")
   } else { # minimal_lstar / minimal_estar
     # h = a0 + (b10 + d1*G)x1 + (b20 + d2*G)x2 + (b30 + d3*G)X1
     X <- cbind(1, x1, G*x1, x2, G*x2, X1, G*X1)
     colnames(X) <- c("a0","b10","d1","b20","d2","b30","d3")
   }
   
   keep <- intersect(colnames(X), names(coef_new))
   fit  <- stats::.lm.fit(X[, keep, drop=FALSE], y)
   coef_new[keep] <- fit$coefficients
   resid <- y - drop(X[, keep, drop=FALSE] %*% fit$coefficients)
   
   # Return a STAR-like list sufficient for .build_F3/simulate_path/Wald
   list(
     spec  = fit0$spec,
     z_var = fit0$z_var,
     v     = fit0$v %||% c(1,0,0),
     coef  = coef_new,
     resid = resid
   )
 }
 
 # Safe call into classifier (skip fallback when available)
 .classify_fast <- function(path_df, fit, delta, opy, burn_in) {
   fmls <- try(names(formals(classify_dynamics)), silent = TRUE)
   if (!inherits(fmls, "try-error") && "linearise_fallback" %in% fmls) {
     classify_dynamics(path_df = path_df, fit = fit, delta = delta,
                       obs_per_year_real = opy, burn_in = burn_in,
                       linearise_fallback = FALSE)
   } else {
     classify_dynamics(path_df = path_df, fit = fit, delta = delta,
                       obs_per_year_real = opy, burn_in = burn_in)
   }
 }
 
 # Uniform spectral peak estimator (Daniell via spec.pgram) 
 .parabolic_peak <- function(x, y, i) {
   n <- length(x)
   if (i <= 1L || i >= n) return(x[i])
   x1 <- x[i-1]; x2 <- x[i]; x3 <- x[i+1]
   y1 <- y[i-1]; y2 <- y[i]; y3 <- y[i+1]
   denom <- (y1 - 2*y2 + y3)
   if (!is.finite(denom) || abs(denom) < .Machine$double.eps) return(x2)
   delta <- 0.5 * (y1 - y3) / denom
   delta <- max(-1, min(1, delta))
   x2 + delta * (x3 - x2)
 }
 
 .pad_pow2_to_at_least <- function(n, N_target) {
   Nt <- suppressWarnings(as.numeric(N_target))
   if (!is.finite(Nt) || Nt < n) Nt <- 2^ceiling(log2(n))
   p  <- ceiling(log2(Nt)) - ceiling(log2(n))
   max(0L, as.integer(p))
 }
 
 dominant_freq_spec <- function(x, opy,
                                band_quarters = c(6L, 32L),
                                N_plot = 2048,
                                spans_plot = 13,
                                taper_plot = 0.10,
                                detrend = TRUE) {
   x <- as.numeric(x); x <- x[is.finite(x)]
   if (length(x) < 40L)
     return(list(f_star=NA_real_, T_obs=NA_real_, T_quarters=NA_real_,
                 T_years=NA_real_, idx=NA_integer_))
   pad <- .pad_pow2_to_at_least(length(x), N_plot)
   sp  <- spec.pgram(x, spans = spans_plot, pad = pad,
                     taper = taper_plot, detrend = detrend, plot = FALSE)
   f <- as.numeric(sp$freq); S <- as.numeric(sp$spec)
   ok <- is.finite(f) & is.finite(S) & f > 0 & f < 0.5
   f <- f[ok]; S <- S[ok]
   if (length(f) < 5L)
     return(list(f_star=NA_real_, T_obs=NA_real_, T_quarters=NA_real_,
                 T_years=NA_real_, idx=NA_integer_))
   lo_q <- band_quarters[1]; hi_q <- band_quarters[2]
   stopifnot(is.finite(lo_q), is.finite(hi_q), lo_q > 0, hi_q > lo_q)
   f_lo <- 4/(opy*hi_q); f_hi <- 4/(opy*lo_q)
   keep <- (f >= f_lo) & (f <= f_hi)
   if (!any(keep))
     return(list(f_star=NA_real_, T_obs=NA_real_, T_quarters=NA_real_,
                 T_years=NA_real_, idx=NA_integer_))
   fk <- f[keep]; Sk <- S[keep]
   imax <- which.max(Sk)
   f_star <- .parabolic_peak(fk, Sk, imax)
   if (!is.finite(f_star)) f_star <- fk[imax]
   f_star <- max(min(f_star, max(fk)), min(fk))
   T_obs   <- 1 / f_star
   T_quart <- (1 / f_star) * (4/opy)
   T_years <- T_obs / opy
   list(f_star=f_star, T_obs=T_obs, T_quarters=T_quart, T_years=T_years, idx=imax)
 }
 
 # Tail extractor (small burn-in for spectra)
 tail_after_burn <- function(path_df, burn_in_small, tail_len = 800L, min_len = 40L) {
   x <- as.numeric(path_df$h); n <- length(x)
   if (!is.finite(burn_in_small)) burn_in_small <- 0L
   start <- min(n, max(1L, burn_in_small + 1L))
   xt <- x[start:n]; xt <- xt[is.finite(xt)]
   if (length(xt) > tail_len) xt <- utils::tail(xt, tail_len)
   if (length(xt) < min_len) return(numeric(0))
   xt
 }
 
 # Empirical spectral period (data)
 empirical_spectral_period <- function(h, opy,
                                       band_quarters = c(6L,32L),
                                       N_plot=2048, spans_plot=13, taper_plot=0.10) {
   out <- dominant_freq_spec(h, opy, band_quarters, N_plot, spans_plot, taper_plot, detrend=TRUE)
   list(point_obs = out$T_obs,
        point_years = out$T_years,
        details = out)
 }
 
 # Deterministic spectral period per nonlinear model
 deterministic_spectral_by_model <- function(models,
                                             init_h, init_hlag, init_H,
                                             horizon, delta,
                                             opy,
                                             band_quarters = c(6L,32L),
                                             burn_in_small = 4L * opy,
                                             tail_len = 800L,
                                             N_plot=2048, spans_plot=13, taper_plot=0.10) {
   keep <- setdiff(names(models), c("AR2","Linear"))
   if (!length(keep)) stop("No nonlinear models in `models`.")
   purrr::map_dfr(keep, function(m) {
     fit  <- models[[m]]
     path <- simulate_path(m, fit, init_h, init_hlag, init_H, horizon, delta)
     xtail <- tail_after_burn(path, burn_in_small, tail_len)
     if (length(xtail) == 0L)
       return(tibble::tibble(Model = m, Det_T_obs = NA_real_, Det_T_years = NA_real_))
     est <- dominant_freq_spec(xtail, opy, band_quarters, N_plot, spans_plot, taper_plot, detrend=TRUE)
     tibble::tibble(Model = m, Det_T_obs = est$T_obs, Det_T_years = est$T_years)
   })
 }
 
 # LC + Wald–gated bootstrap spectral periods (FAST STAR)
 bootstrap_spectral_periods <- function(models, hyps_list,
                                        B = 250L,
                                        bootstrap_scheme = c("iid","wild"),
                                        init_h, init_hlag, init_H,
                                        delta, horizon, opy,
                                        band_quarters = c(6L,32L),
                                        burn_in_class = NULL,
                                        burn_in_small = 10L * opy,
                                        tail_len = 800L,
                                        N_plot = 2048, spans_plot = 13, taper_plot = 0.10,
                                        verbose = TRUE,
                                        max_draws_per_model = 10L * B,
                                        accept_explosive = FALSE,
                                        refit_star = c("linear","none","full")) {
   bootstrap_scheme <- match.arg(bootstrap_scheme)
   refit_star <- match.arg(refit_star)
   if (is.null(burn_in_class)) burn_in_class <- get0("burn_in") %||% (20L*opy)
   
   mods <- setdiff(names(models), c("AR2","Linear"))
   out  <- setNames(vector("list", length(mods)), mods)
   diag <- setNames(vector("list", length(mods)), mods)
   
   for (mod in mods) {
     fit0 <- models[[mod]]
     
     nT <- if (!is_star_fit(fit0)) {
       nrow(model.frame(fit0))
     } else {
       nrow(get0("df_samp"))
     }
     if (!is.finite(nT)) stop("Cannot infer nT; set df_samp or pass model.frame for ", mod)
     
     Ts <- numeric(0)
     kept <- 0L; draws <- 0L; lc_kept <- 0L; wald_pass <- 0L
     hyps <- hyps_list[[mod]]
     
     if (verbose) message(sprintf("[BOOT %s] target kept B=%d", mod, B))
     pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
     
     while (kept < B && draws < max_draws_per_model) {
       draws <- draws + 1L
       
       # 1) Bootstrap sample (mean-law, residual resampling)
       dfb <- try(gen_boot_df(mod, fit0, nT, init_h, init_hlag, init_H, delta,
                              scheme = bootstrap_scheme), silent = TRUE)
       if (inherits(dfb, "try-error")) { utils::setTxtProgressBar(pb, kept); next }
       
       # 2) Re-estimate (LM: OLS; STAR: fast OLS-only update or reuse baseline)
       fit_b <- if (!is_star_fit(fit0)) {
         form <- formula(fit0)
         try(stats::lm(form, data = dfb), silent = TRUE)
       } else if (refit_star == "none") {
         fit0
       } else if (refit_star == "linear") {
         try(fast_refit_star(dfb, fit0), silent = TRUE)
       } else { # "full" (slow; avoid unless necessary)
         try(estimate_star(
           dfb,
           spec  = fit0$spec,
           z_var = fit0$z_var,
           v     = fit0$v %||% c(1,0,0),
           starts_gamma = c(1,5),
           starts_c     = "median",
           trace        = FALSE,
           maxit        = 800
         ), silent = TRUE)
       }
       if (inherits(fit_b, "try-error")) { utils::setTxtProgressBar(pb, kept); next }
       
       # 3) Wald gate first (use precomputed V from hyps_list for STAR)
       wald_ok <- if (is.null(hyps)) TRUE else wald_reject_5(fit_b, hyps)
       if (!isTRUE(wald_ok)) { utils::setTxtProgressBar(pb, kept); next }
       wald_pass <- wald_pass + 1L
       
       # 4) Simulate path; optionally discard explosive
       path_b <- try(simulate_path(mod, fit_b, init_h, init_hlag, init_H, horizon, delta),
                     silent = TRUE)
       if (inherits(path_b, "try-error")) { utils::setTxtProgressBar(pb, kept); next }
       if (!accept_explosive && !finite_ok(path_b$h)) { utils::setTxtProgressBar(pb, kept); next }
       
       # 5) Fast classification (skip linearisation fallback)
       cls_b <- .classify_fast(path_b, fit_b, delta, opy, burn_in = burn_in_class)
       if (!identical(cls_b, "LC")) { utils::setTxtProgressBar(pb, kept); next }
       lc_kept <- lc_kept + 1L
       
       # 6) Tail spectrum
       xtail <- tail_after_burn(path_b, burn_in_small, tail_len)
       if (length(xtail) == 0L) { utils::setTxtProgressBar(pb, kept); next }
       
       est <- dominant_freq_spec(xtail, opy, band_quarters, N_plot, spans_plot,
                                 taper_plot, detrend = TRUE)
       Ts  <- c(Ts, est$T_obs)
       kept <- kept + 1L
       utils::setTxtProgressBar(pb, kept)
       
       if (verbose && (kept %% 25L == 0L))
         message(sprintf("%s: kept=%d (draws=%d; Wald=%d; LC=%d)",
                         mod, kept, draws, wald_pass, lc_kept))
     }
     close(pb)
     if (verbose) message(sprintf("%s: FINAL kept=%d (draws=%d; Wald=%d; LC=%d)",
                                  mod, kept, draws, wald_pass, lc_kept))
     
     out[[mod]] <- Ts
     diag[[mod]] <- list(kept=kept, draws=draws, wald_pass=wald_pass, lc_kept=lc_kept)
   }
   
   list(T_obs = out, diag = diag)
 }
 
 summarize_boot_T <- function(Ts_obs, opy) {
   if (!length(Ts_obs)) return(list(
     q = setNames(rep(NA_real_,5), c("16%","50%","84%","2.5%","97.5%")),
     q_years = setNames(rep(NA_real_,5), c("16%","50%","84%","2.5%","97.5%"))
   ))
   qs <- stats::quantile(Ts_obs, probs=c(.16,.5,.84,.025,.975), type=8, names=TRUE)
   list(q = qs, q_years = qs / opy)
 }
 
 # Glue: run everything and build a compact table
 compute_cycle_lengths_spectral <- function(models, hyps_list,
                                            data_h,
                                            init_h, init_hlag, init_H,
                                            horizon, delta, opy,
                                            band_quarters = c(6L,32L),
                                            burn_in_small = 4L * opy,
                                            burn_in_class = NULL,
                                            tail_len = 800L,
                                            B_boot = 100L,
                                            bootstrap_scheme = "iid",
                                            N_plot=2048, spans_plot=13, taper_plot=0.10,
                                            verbose = TRUE,
                                            refit_star = c("linear","none","full")) {
   
   refit_star <- match.arg(refit_star)
   
   # 0) Empirical
   emp <- empirical_spectral_period(data_h, opy, band_quarters, N_plot, spans_plot, taper_plot)
   
   # 1) Deterministic (nonlinear only)
   det_tbl <- deterministic_spectral_by_model(models, init_h, init_hlag, init_H,
                                              horizon, delta, opy,
                                              band_quarters, burn_in_small, tail_len,
                                              N_plot, spans_plot, taper_plot)
   
   # 2) Bootstrap LC+Wald spectral periods (FAST STAR)
   boot <- bootstrap_spectral_periods(
     models, hyps_list, B = B_boot,
     bootstrap_scheme = bootstrap_scheme,
     init_h = init_h, init_hlag = init_hlag, init_H = init_H,
     delta = delta, horizon = horizon, opy = opy,
     band_quarters = band_quarters,
     burn_in_class = burn_in_class,
     burn_in_small = burn_in_small, tail_len = tail_len,
     N_plot = N_plot, spans_plot = spans_plot, taper_plot = taper_plot,
     verbose = verbose,
     max_draws_per_model = 10L * B_boot,
     accept_explosive = FALSE,
     refit_star = refit_star
   )
   
   # Summaries per model
   mods <- setdiff(names(models), c("AR2","Linear"))
   boot_summ <- purrr::map_dfr(mods, function(m){
     s <- summarize_boot_T(boot$T_obs[[m]], opy)
     tibble::tibble(
       Model = m,
       Boot_T_q16   = unname(s$q["16%"]),
       Boot_T_med   = unname(s$q["50%"]),
       Boot_T_q84   = unname(s$q["84%"]),
       Boot_T_q025  = unname(s$q["2.5%"]),
       Boot_T_q975  = unname(s$q["97.5%"]),
       Boot_T_q16_y = unname(s$q_years["16%"]),
       Boot_T_med_y = unname(s$q_years["50%"]),
       Boot_T_q84_y = unname(s$q_years["84%"]),
       Boot_T_q025_y= unname(s$q_years["2.5%"]),
       Boot_T_q975_y= unname(s$q_years["97.5%"])
     )
   })
   
   # 3) Final table
   to_yrs <- function(Tobs) ifelse(is.finite(Tobs), Tobs/opy, NA_real_)
   res_tbl <- det_tbl %>% dplyr::left_join(boot_summ, by="Model") %>%
     dplyr::mutate(Det_T_years = to_yrs(Det_T_obs))
   
   emp_row <- tibble::tibble(
     Model = "Data",
     Det_T_obs = NA_real_, Det_T_years = NA_real_,
     Boot_T_q16 = NA_real_, Boot_T_med = NA_real_, Boot_T_q84 = NA_real_,
     Boot_T_q025 = NA_real_, Boot_T_q975 = NA_real_,
     Boot_T_q16_y = NA_real_, Boot_T_med_y = NA_real_, Boot_T_q84_y = NA_real_,
     Boot_T_q025_y = NA_real_, Boot_T_q975_y = NA_real_
   )
   tbl <- dplyr::bind_rows(emp_row, res_tbl)
   
   fmt2 <- function(x) ifelse(is.finite(x), sprintf("%.2f", x), "")
   rng2 <- function(a,b) if (is.finite(a) && is.finite(b)) sprintf("[%.2f, %.2f]", a,b) else ""
   
   final_tbl <- tibble::tibble(
     Model                 = tbl$Model,
     `Empirical spec (yrs)`= c(sprintf("%.2f", emp$point_years), rep("", nrow(tbl)-1L)),
     `Deterministic (obs)` = vapply(tbl$Det_T_obs,   fmt2, character(1)),
     `Deterministic (yrs)` = vapply(tbl$Det_T_years, fmt2, character(1)),
     `Boot med (obs)`      = vapply(tbl$Boot_T_med,  fmt2, character(1)),
     `Boot 68% (obs)`      = vapply(seq_len(nrow(tbl)), function(i)
       rng2(tbl$Boot_T_q16[i], tbl$Boot_T_q84[i]), character(1)),
     `Boot med (yrs)`      = vapply(tbl$Boot_T_med_y, fmt2, character(1)),
     `Boot 68% (yrs)`      = vapply(seq_len(nrow(tbl)), function(i)
       rng2(tbl$Boot_T_q16_y[i], tbl$Boot_T_q84_y[i]), character(1))
   )
   
   list(empirical = emp,
        deterministic = det_tbl,
        bootstrap = boot,
        table = final_tbl)
 }
 
 
 ##### CALL 
 res <- compute_cycle_lengths_spectral(
   models  = models,
   hyps_list = hyps_list,
   data_h  = data_merged_real$h,
   init_h  = initial_h,
   init_hlag = initial_h_lag,
   init_H  = H_tm1_start,
   horizon = horizon,
   delta   = delta,
   opy     = obs_per_year_real,
   band_quarters = c(6L, 60L),
   burn_in_small = 20 * obs_per_year_real,
   burn_in_class = 20 * obs_per_year_real,
   tail_len = 800L,
   B_boot  = 100L,
   bootstrap_scheme = "iid",
   N_plot  = 2048,
   spans_plot = 13,
   taper_plot = 0.10,
   verbose = TRUE,
   refit_star = "linear"   # <<< FAST STAR (reuse gamma/c, OLS for linear params)
 )
 
 cat("\n=== Cycle Lengths (Spectral) ===\n")
 print(knitr::kable(res$table, align = "lrrrrrrr"))
 
 # BAR SNAPSHOT: Cycle lengths (years) 
 plot_cycle_length_snapshot <- function(res, order_by = c("deterministic","closeness_to_empirical")){
   order_by <- match.arg(order_by)
   
   tbl <- res$table %>%
     dplyr::filter(Model != "Data") %>%
     dplyr::transmute(
       Model,
       Det_Y = suppressWarnings(as.numeric(`Deterministic (yrs)`)),
       Boot_med_Y  = suppressWarnings(as.numeric(`Boot med (yrs)`)),
       Boot_low_Y  = suppressWarnings(as.numeric(sub("\\[(.*),.*\\]","\\1", `Boot 68% (yrs)`))),
       Boot_high_Y = suppressWarnings(as.numeric(sub("\\[.*,(.*)\\]","\\1", `Boot 68% (yrs)`)))
     )
   
   emp_years <- res$empirical$point_years
   
   if (order_by == "closeness_to_empirical" && is.finite(emp_years)) {
     tbl <- tbl %>% dplyr::mutate(closeness = abs(Det_Y - emp_years)) %>%
       dplyr::arrange(closeness, Det_Y)
   } else {
     tbl <- tbl %>% dplyr::arrange(Det_Y)
   }
   
   tbl <- tbl %>% dplyr::mutate(Model = forcats::fct_inorder(Model))
   
   ggplot(tbl, aes(x = Det_Y, y = Model)) +
     geom_col(width = 0.6, alpha = 0.45, na.rm = TRUE) +
     geom_errorbarh(aes(xmin = Boot_low_Y, xmax = Boot_high_Y), height = 0.24, size = 0.7, na.rm = TRUE) +
     geom_point(aes(x = Boot_med_Y), size = 2.2, na.rm = TRUE) +
     { if (is.finite(emp_years)) geom_vline(xintercept = emp_years, linetype = 2) else NULL } +
     labs(
       x = "Cycle length (years)",
       y = NULL,
       title = "Cycle Lengths: Deterministic (bars), Bootstrap median±68% (points/whiskers)",
       subtitle = if (is.finite(emp_years))
         sprintf("Dashed line = empirical spectral period (%.2f years)", emp_years) else
           "Empirical spectral period unavailable"
     ) +
     theme_minimal(base_size = 12) +
     theme(panel.grid.major.y = element_blank(),
           plot.title.position = "plot")
 }
 

 p <- plot_cycle_length_snapshot(res, order_by = "closeness_to_empirical")
 print(p)
 

  
  
  









  
  
  
  


#################################################################################################
# ROBUSTNESS CHECKS #############################################################################
#################################################################################################

# helpers (steady state + eigenvalues) 

# 0) One palette to rule them all 
model_pal <- c(
  "Linear"       = "#1f77b4",  # blue
  "Minimal"      = "#d62728",  # red
  "Intermediate" = "#2ca02c",  # green
  "Full"         = "#9467bd"   # purple
)

# small helper: enforce common levels so colors map consistently
as_model_factor <- function(x) factor(x, levels = names(model_pal))


`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L || (length(a) == 1L && is.na(a))) b else a
}

# robust root finder for 1D equations on [L,U]:
uniroot_all_grid <- function(f, L=-10, U=10, m=400) {
  xs <- seq(L, U, length.out = m+1)
  vals <- sapply(xs, f)
  okL  <- is.finite(head(vals, -1))
  okR  <- is.finite(tail(vals, -1))
  cross <- head(vals, -1) * tail(vals, -1) <= 0
  idx <- which(cross & okL & okR)
  roots <- numeric(0)   # <-- numeric, not NULL
  for (i in idx) {
    a <- xs[i]; b <- xs[i+1]
    fa <- vals[i]; fb <- vals[i+1]
    if (!is.finite(fa) || !is.finite(fb) || fa == fb) next
    r <- try(uniroot(f, lower=a, upper=b)$root, silent=TRUE)
    if (!inherits(r, "try-error")) roots <- c(roots, as.numeric(r))
  }
  if (!length(roots)) return(numeric(0))  # <-- skip rounding on empty
  unique(round(roots, 10))
}

# build X_{t} with **truncated** kernel (N quarters)
build_X_trunc <- function(h, delta, N=40) {
  w <- delta * (1 - delta)^(0:(N-1))
  X <- rep(NA_real_, length(h))
  for (t in seq_along(h)) {
    k <- min(t, N)
    X[t] <- sum(w[1:k] * rev(h[(t-k+1):t]))
  }
  X
}

# law of motion h_t = f(h_{t-1}, h_{t-2}, X_{t-1}) (+ innovation added outside)
law_h_eval <- function(coefs, h1, h2, X1, model=c("Linear","Minimal","Intermediate")) {
  model <- match.arg(model)
  b <- function(nm) coefs[nm] %||% 0
  out <- b("(Intercept)") + b("x_lag1")*h1 + b("x_lag2")*h2 + b("X_lag1")*X1
  if (model != "Linear") out <- out + b("x_lag1_cub")*h1^3
  if (model == "Intermediate") {
    out <- out + (b("x_lag2_cub")*h2^3 + b("X_lag1_cub")*X1^3 +
                    b("x_lag1_sq_X_lag1")*(h1^2)*X1 + b("x_lag1_X_lag1_sq")*h1*(X1^2))
  }
  out
}

# max-modulus eigenvalue at the **steady state** of the estimated model
rho_ss <- function(coefs, delta, model=c("Linear","Minimal","Intermediate"), mean_h_reference=0) {
  model <- match.arg(model)
  b <- function(nm) coefs[nm] %||% 0
  
  # steady state: h* = f(h*, h*, h*)
  f_minus_id <- function(h) law_h_eval(coefs, h, h, h, model) - h
  roots <- uniroot_all_grid(f_minus_id, L=-10, U=10, m=600)
  # choose the root closest to the sample mean of h
  hbar <- if (length(roots)) roots[which.min(abs(roots - mean_h_reference))] else 0
  
  # Jacobian entries at (h*, h*, X*) with X*=h*
  a11 <- b("x_lag1")
  a12 <- b("x_lag2")
  a13 <- b("X_lag1")
  if (model != "Linear") a11 <- a11 + 3*(b("x_lag1_cub"))*(hbar^2)
  # Intermediate adds **no extra** first derivatives w.r.t h1 at steady state beyond the h1-cubic
  
  J <- matrix(c(
    a11,       a12,       a13,
    1,         0,         0,
    delta*a11, delta*a12, (1-delta)+delta*a13
  ), 3, 3, byrow=TRUE)
  
  max(Mod(eigen(J, only.values=TRUE)$values))
}

# simulate one pseudo-sample for a **given model** with bootstrapped residuals
simulate_nonlinear_sample <- function(coefs, T, delta, N=40, resid_pool,
                                      model=c("Linear","Minimal","Intermediate"),
                                      burnin=300, cap=Inf) {
  model <- match.arg(model)
  e_star <- sample(resid_pool - mean(resid_pool, na.rm=TRUE), size=T+burnin, replace=TRUE)
  
  # initial conditions (washed out by burn-in)
  h <- numeric(T+burnin)
  X <- numeric(T+burnin)
  # start at zero; 
  h[1:2] <- 0
  X[1] <- 0
  
  for (t in 3:(T+burnin)) {
    # X_{t-1} with **truncation N** over the last N realized h's
    Xlag <- sum(delta * (1 - delta)^(0:(min(t-1, N)-1)) * rev(h[(t-1 - (min(t-1,N)-1)):(t-1)]))
    mu <- law_h_eval(coefs, h[t-1], h[t-2], Xlag, model)
    h[t] <- mu + e_star[t]
    if (abs(h[t]) > cap || !is.finite(h[t])) stop("explosive")
  }
  # drop burn-in
  h <- h[(burnin+1):(T+burnin)]
  h
}

# from h -> features used in estimation (lags and truncated X_{t-1})
make_features_from_h <- function(h, delta, N=40) {
  X <- build_X_trunc(h, delta, N)
  df <- data.frame(
    h       = h,
    x_lag1  = dplyr::lag(h, 1),
    x_lag2  = dplyr::lag(h, 2),
    X_lag1  = dplyr::lag(X, 1)
  )
  df <- tidyr::drop_na(df)
  dplyr::mutate(df,
                x_lag1_sq            = x_lag1^2,
                x_lag2_sq            = x_lag2^2,
                X_lag1_sq            = X_lag1^2,
                x_lag1_cub           = x_lag1^3,
                x_lag2_cub           = x_lag2^3,
                X_lag1_cub           = X_lag1^3,
                x_lag1_x_lag2        = x_lag1 * x_lag2,
                x_lag1_X_lag1        = x_lag1 * X_lag1,
                x_lag2_X_lag1        = x_lag2 * X_lag1,
                x_lag1_sq_x_lag2     = (x_lag1^2) * x_lag2,
                x_lag1_x_lag2_sq     = x_lag1 * (x_lag2^2),
                x_lag1_sq_X_lag1     = (x_lag1^2) * X_lag1,
                x_lag2_sq_X_lag1     = (x_lag2^2) * X_lag1,
                x_lag1_X_lag1_sq     = x_lag1 * (X_lag1^2),
                x_lag2_X_lag1_sq     = x_lag2 * (X_lag1^2),
                x_lag1_x_lag2_X_lag1 = x_lag1 * x_lag2 * X_lag1
  )
}


# Only defined if missing; relies on earlier .G_logistic() and compute_G_from_df()
if (!exists(".get_v", mode="function")) {
  .get_v <- function(star_fit) star_fit$v %||% c(1,0,0)
}
if (!exists("predict_star_min_mu", mode="function")) {
  # mean law: h_t = a0 + (b10 + d1*G)*x1 + (b20 + d2*G)*x2 + (b30 + d3*G)*X1
  predict_star_min_mu <- function(df, beta, star_fit) {
    need <- c("a0","b10","d1","b20","d2","b30","d3")
    stopifnot(all(need %in% names(beta)))
    G <- compute_G_from_df(df, star_fit)
    with(df, beta["a0"] +
           (beta["b10"] + beta["d1"]*G)*x_lag1 +
           (beta["b20"] + beta["d2"]*G)*x_lag2 +
           (beta["b30"] + beta["d3"]*G)*X_lag1)
  }
}
if (!exists("simulate_star_min_sample", mode="function")) {
  simulate_star_min_sample <- function(star_fit, T, delta, resid_pool,
                                       burnin=300, cap=1e6) {
    th <- star_fit$coef
    beta <- c(a0 = th["a0"], b10 = th["b10"], d1 = th["d1"],
              b20 = th["b20"], d2 = th["d2"], b30 = th["b30"], d3 = th["d3"])
    v <- .get_v(star_fit)
    e_star <- sample(resid_pool - mean(resid_pool, na.rm=TRUE), T+burnin, TRUE)
    h <- numeric(T+burnin); H <- numeric(T+burnin)
    # start at zero
    h[1:2] <- 0; H[1] <- 0
    for (t in 3:(T+burnin)) {
      h1 <- h[t-1]; h2 <- h[t-2]; H1 <- H[t-1]
      z  <- if (identical(star_fit$z_var,"x_lag1")) h1 else v[1]*h1 + v[2]*h2 + v[3]*H1
      G  <- .G_logistic(z, eta = unname(th["eta"]), c = unname(th["c"]))
      mu <- beta["a0"] +
        (beta["b10"] + beta["d1"]*G)*h1 +
        (beta["b20"] + beta["d2"]*G)*h2 +
        (beta["b30"] + beta["d3"]*G)*H1
      ht <- mu + e_star[t]
      Ht <- (1-delta)*H1 + delta*ht
      if (!is.finite(ht) || abs(ht) > cap) stop("explosive")
      h[t] <- ht; H[t] <- Ht
    }
    h[(burnin+1):(T+burnin)]
  }
}
if (!exists("rho_ss_star_min", mode="function")) {
  # Analytic Jacobian at steady state (fast; no numDeriv)
  rho_ss_star_min <- function(beta, star_fit, delta, mean_h_reference = 0) {
    need <- c("a0","b10","d1","b20","d2","b30","d3")
    stopifnot(all(need %in% names(beta)))
    eta <- unname(star_fit$coef["eta"]); c0 <- unname(star_fit$coef["c"])
    v   <- .get_v(star_fit)
    
    # steady state s*: solve F(s,s,s)=s
    f3 <- function(s) {
      z <- if (identical(star_fit$z_var,"x_lag1")) s else (v[1]+v[2]+v[3])*s
      G <- .G_logistic(z, eta, c0)
      beta["a0"] +
        (beta["b10"] + beta["d1"]*G)*s +
        (beta["b20"] + beta["d2"]*G)*s +
        (beta["b30"] + beta["d3"]*G)*s
    }
    g  <- function(s) f3(s) - s
    roots <- uniroot_all_grid(g, L = -10, U = 10, m = 600)
    sstar <- if (length(roots)) roots[ which.min(abs(roots - mean_h_reference)) ] else 0
    
    # derivatives at s*
    zstar <- if (identical(star_fit$z_var,"x_lag1")) sstar else (v[1]+v[2]+v[3])*sstar
    Gs    <- .G_logistic(zstar, eta, c0)
    k     <- exp(eta)                   # slope parameter
    gprime <- k * Gs * (1 - Gs)
    
    dz_dh1 <- if (identical(star_fit$z_var,"x_lag1")) 1 else v[1]
    dz_dh2 <- if (identical(star_fit$z_var,"x_lag1")) 0 else v[2]
    dz_dH1 <- if (identical(star_fit$z_var,"x_lag1")) 0 else v[3]
    S      <- sstar
    Sdot   <- (beta["d1"]*S + beta["d2"]*S + beta["d3"]*S)  # d1*h1 + d2*h2 + d3*H1 at ss
    
    a11 <- (beta["b10"] + beta["d1"]*Gs) + gprime*dz_dh1*Sdot
    a12 <- (beta["b20"] + beta["d2"]*Gs) + gprime*dz_dh2*Sdot
    a13 <- (beta["b30"] + beta["d3"]*Gs) + gprime*dz_dH1*Sdot
    
    J <- matrix(c(
      a11, a12, a13,
      1,   0,   0,
      delta*a11, delta*a12, (1-delta)+delta*a13
    ), 3, 3, byrow = TRUE)
    
    max(Mod(eigen(J, only.values = TRUE)$values))
  }
}






##############################################################################
# Confidence Interval and Posterior Draws
##############################################################################

# 1) Shared setup: initial state, forecast dates, δ, number of draws
# last in-sample date (origin), strictly before the first forecast date
origin_row <- max(which(data_merged_real$date < from))
stopifnot(is.finite(origin_row))

initial_state <- with(data_merged_real[origin_row, ], list(
  x_tm1 = h,        # h_t0
  x_tm2 = x_lag1,   # h_{t0-1}
  X_tm1 = X_lag1    # X_t0
))
# dates for the simulated path remain the same:
horizon_ci <- if (exists("forecast_dates")) {
  length(forecast_dates)
} else if (exists("horizon_data")) {
  horizon_data
} else {
  100L  # fallback demo length
}
# first forecast point corresponds to 'from'
forecast_dates_ci <- seq(from = from, by = obs_unit_real, length.out = horizon_ci)

delta   <- 1 - (1 - 0.05)^(4 / obs_per_year_real)
N_draws <- 1000

# 2) Core deterministic simulation
simulate_det <- function(beta_vec) {
  law_h <- function(h1, h2, H1) {
    regs <- list(
      `(Intercept)` = 1,
      x_lag1        = h1,
      x_lag2        = h2,
      X_lag1        = H1
    )
    # optional terms — only canonical names
    if ("x_lag1_sq"            %in% names(beta_vec)) regs$x_lag1_sq            <- h1^2
    if ("x_lag2_sq"            %in% names(beta_vec)) regs$x_lag2_sq            <- h2^2
    if ("X_lag1_sq"            %in% names(beta_vec)) regs$X_lag1_sq            <- H1^2
    if ("x_lag1_cub"           %in% names(beta_vec)) regs$x_lag1_cub           <- h1^3
    if ("x_lag2_cub"           %in% names(beta_vec)) regs$x_lag2_cub           <- h2^3
    if ("X_lag1_cub"           %in% names(beta_vec)) regs$X_lag1_cub           <- H1^3
    if ("x_lag1_x_lag2"        %in% names(beta_vec)) regs$x_lag1_x_lag2        <- h1 * h2
    if ("x_lag1_X_lag1"        %in% names(beta_vec)) regs$x_lag1_X_lag1        <- h1 * H1
    if ("x_lag2_X_lag1"        %in% names(beta_vec)) regs$x_lag2_X_lag1        <- h2 * H1
    if ("x_lag1_sq_x_lag2"     %in% names(beta_vec)) regs$x_lag1_sq_x_lag2     <- h1^2 * h2
    if ("x_lag1_x_lag2_sq"     %in% names(beta_vec)) regs$x_lag1_x_lag2_sq     <- h1 * h2^2
    if ("x_lag1_sq_X_lag1"     %in% names(beta_vec)) regs$x_lag1_sq_X_lag1     <- h1^2 * H1
    if ("x_lag2_sq_X_lag1"     %in% names(beta_vec)) regs$x_lag2_sq_X_lag1     <- h2^2 * H1
    if ("x_lag1_X_lag1_sq"     %in% names(beta_vec)) regs$x_lag1_X_lag1_sq     <- h1 * H1^2
    if ("x_lag2_X_lag1_sq"     %in% names(beta_vec)) regs$x_lag2_X_lag1_sq     <- h2 * H1^2
    if ("x_lag1_x_lag2_X_lag1" %in% names(beta_vec)) regs$x_lag1_x_lag2_X_lag1 <- h1 * h2 * H1
    
    use <- intersect(names(beta_vec), names(regs))
    sum(unlist(regs[use]) * beta_vec[use])
  }
  
  simulate_core_recursive(
    h_tm1 = initial_state$x_tm1,
    h_tm2 = initial_state$x_tm2,
    H_tm1 = initial_state$X_tm1,
    horizon = horizon_ci,
    delta = delta,
    law_h = law_h
  )$h
}

# 4) CI‑plot helper using posterior draws‑plot helper using posterior draws
make_ci_plot <- function(model_obj, name) {
  beta_hat <- coef(model_obj)
  cov_beta <- vcov(model_obj)
  
  # draw from posterior
  post_draws <- MASS::mvrnorm(n = N_draws, mu = beta_hat, Sigma = cov_beta)
  
  # simulate each draw
  sims <- sapply(seq_len(N_draws), function(i) simulate_det(post_draws[i, ]))
  
  # build 95% bands
  lo95 <- apply(sims, 1, quantile, probs = 0.025, na.rm = TRUE)
  hi95 <- apply(sims, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  conf_df <- data.frame(
    date = forecast_dates_ci,
    lo95 = lo95,
    hi95 = hi95
  )
  
  # deterministic path
  det_path <- simulate_det(beta_hat)
  det_df <- data.frame(
    date = forecast_dates_ci,
    det  = det_path
  )
  
  # real historical series up to last forecast date
  real_df <- subset(
    data_merged_real,
    date <= max(forecast_dates_ci),
    select = c(date, h)
  )
  names(real_df)[2] <- "value"
  
  # plot
  ggplot() +
    geom_line(data = real_df, aes(date, value), color = "darkgrey", size = 0.8) +
    geom_ribbon(data = conf_df, aes(date, ymin = lo95, ymax = hi95), fill = "grey70", alpha = 0.4) +
    geom_line(data = det_df, aes(date, det), color = "red", size = 1) +
    labs(
      title = paste0(name, " Model – 95% Posterior Band"),
      x = "Date",
      y = "Cyclical Component (h)"
    ) +
    theme_minimal()
}


# STAR_min (FAST posterior bands)
# logistic link used by minimal LSTAR
.G_logistic <- function(z, eta, c) 1/(1 + exp(-exp(eta) * (z - c)))

# small safe getter for v (avoid %||% on vectors)
.get_v <- function(fit) {
  v <- fit$v
  if (is.null(v) || length(v) != 3 || anyNA(v) || !all(is.finite(v))) c(1,0,0) else v
}

# Build G_{t-1} for a dataframe given baseline STAR fit (uses z_var & v if 'combo')
compute_G_from_df <- function(df, star_fit) {
  stopifnot(is.list(star_fit), !is.null(star_fit$spec), star_fit$spec %in% c("minimal_lstar"))
  z <- if (identical(star_fit$z_var, "x_lag1")) {
    df$x_lag1
  } else {
    v <- .get_v(star_fit)
    v[1]*df$x_lag1 + v[2]*df$x_lag2 + v[3]*df$X_lag1
  }
  .G_logistic(z, eta = unname(star_fit$coef["eta"]), c = unname(star_fit$coef["c"]))
}

# OLS “linear part” re-estimation given fixed (eta,c,v); returns beta & cov
ols_linear_part_star_min <- function(df, star_fit) {
  G <- compute_G_from_df(df, star_fit)
  X <- cbind(
    a0  = 1,
    b10 = df$x_lag1, d1 = G*df$x_lag1,
    b20 = df$x_lag2, d2 = G*df$x_lag2,
    b30 = df$X_lag1, d3 = G*df$X_lag1
  )
  y <- df$h
  ok <- stats::complete.cases(X, y)
  X <- X[ok, , drop=FALSE]; y <- y[ok]
  
  fit <- stats::.lm.fit(X, y)
  
  # name the coefficients
  beta_hat <- setNames(drop(fit$coefficients), colnames(X))
  
  resid  <- y - drop(X %*% beta_hat)
  sigma2 <- sum(resid^2) / (nrow(X) - ncol(X))
  V <- sigma2 * solve(crossprod(X))
  dimnames(V) <- list(colnames(X), colnames(X))   # keep names consistent
  list(beta_hat = beta_hat, V = V)
}

simulate_det_STAR_min <- function(beta_vec, star_fit, init_state, horizon, delta) {
  # enforce expected order; error out if names missing
  need <- c("a0","b10","d1","b20","d2","b30","d3")
  if (is.null(names(beta_vec)) || !all(need %in% names(beta_vec))) {
    stop("simulate_det_STAR_min: named coefficients required: ", paste(need, collapse=", "))
  }
  beta_vec <- beta_vec[need]
  
  eta <- unname(star_fit$coef["eta"]); c0 <- unname(star_fit$coef["c"])
  v   <- .get_v(star_fit)
  
  law_h <- function(h1, h2, H1) {
    z <- if (identical(star_fit$z_var, "x_lag1")) h1 else v[1]*h1 + v[2]*h2 + v[3]*H1
    G <- .G_logistic(z, eta, c0)
    beta_vec["a0"] +
      (beta_vec["b10"] + beta_vec["d1"]*G)*h1 +
      (beta_vec["b20"] + beta_vec["d2"]*G)*h2 +
      (beta_vec["b30"] + beta_vec["d3"]*G)*H1
  }
  
  simulate_core_recursive(
    h_tm1 = init_state$x_tm1,
    h_tm2 = init_state$x_tm2,
    H_tm1 = init_state$X_tm1,
    horizon = horizon,
    delta   = delta,
    law_h   = law_h
  )$h
}

make_ci_plot_STAR_min <- function(
    star_fit, name = "STAR_min",
    df_design      = get0("df_samp"),
    init_state     = get0("initial_state"),
    horizon        = get0("horizon_ci"),
    delta          = get0("delta"),
    N_draws        = get0("N_draws"),
    forecast_dates = get0("forecast_dates_ci")
) {
  if (is.null(df_design) || is.null(init_state) || is.null(horizon) ||
      is.null(delta) || is.null(N_draws) || is.null(forecast_dates)) {
    stop("make_ci_plot_STAR_min: missing globals.")
  }
  
  ols <- ols_linear_part_star_min(df_design, star_fit)
  beta_hat <- ols$beta_hat
  V        <- ols$V
  
  det_path <- simulate_det_STAR_min(beta_hat, star_fit, init_state, horizon, delta)
  
  # draws (screened)
  raw_draws <- MASS::mvrnorm(n = N_draws, mu = beta_hat, Sigma = V)
  colnames(raw_draws) <- names(beta_hat)
  
  # keep only locally stable draws (ρ < 1 at the steady state)
  mean_h_ref <- mean(df_design$h, na.rm = TRUE)
  keep <- apply(raw_draws, 1, function(row) {
    b <- setNames(row, names(beta_hat))
    tryCatch(
      rho_ss_star_min(beta = b, star_fit = star_fit, delta = delta,
                      mean_h_reference = mean_h_ref) < 1,
      error = function(...) FALSE
    )
  })
  draws <- raw_draws[keep, , drop = FALSE]
  
  if (nrow(draws) < 50) warning("Few stable STAR_min draws passed the screen: ", nrow(draws))
  
  # simulate each kept draw, with explosion guard
  cap <- max(10, 8 * stats::sd(df_design$h, na.rm = TRUE))  # sensible scale cap (pct pts)
  sim_one <- function(b) {
    s <- try(
      simulate_det_STAR_min(setNames(b, names(beta_hat)), star_fit,
                            init_state, horizon, delta),
      silent = TRUE
    )
    if (inherits(s, "try-error")) return(rep(NA_real_, horizon))
    if (any(!is.finite(s)) || max(abs(s)) > cap) return(rep(NA_real_, horizon))
    s
  }
  sims <- t(apply(draws, 1, sim_one))      # n_draws_kept x horizon
  if (!is.matrix(sims)) sims <- matrix(sims, nrow = 1)
  
  # drop bad/NA paths
  ok_cols <- apply(sims, 1, function(v) all(is.finite(v)))
  sims <- sims[ok_cols, , drop = FALSE]
  
  if (!nrow(sims)) stop("All STAR_min simulations were filtered out (unstable or exploded).")
  
  # pointwise bands 
  lo95 <- apply(sims, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  hi95 <- apply(sims, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  
  conf_df <- data.frame(date = forecast_dates, lo95 = lo95, hi95 = hi95)
  
  # deterministic path (already stable)
  det_path <- simulate_det_STAR_min(beta_hat, star_fit, init_state, horizon, delta)
  det_df  <- data.frame(date = forecast_dates, det = det_path)
  
  real_df <- subset(get0("data_merged_real"), date <= max(forecast_dates), select = c(date, h))
  names(real_df)[2] <- "value"
  
  ggplot2::ggplot() +
    ggplot2::geom_line(data = real_df, ggplot2::aes(date, value), color = "darkgrey", linewidth = 0.8) +
    ggplot2::geom_ribbon(data = conf_df, ggplot2::aes(date, ymin = lo95, ymax = hi95),
                         fill = "grey70", alpha = 0.4) +
    ggplot2::geom_line(data = det_df, ggplot2::aes(date, det), color = "red", linewidth = 1.0) +  # <- red for consistency
    ggplot2::labs(
      title = paste0("Minimal LSTAR) - 95% band)"),
      subtitle = "", # "η, c (and v) held fixed from baseline fit; draws for a0,b10,b20,b30,d1,d2,d3"
      x = "Date", y = "Cyclical Component (h)"
    ) +
    ggplot2::theme_minimal()
}


# 5) Generate and arrange plots
# ADD STAR_min 
plots <- list(
  Minimal       = make_ci_plot(model_minimal,      "Minimal"),
  Intermediate  = make_ci_plot(model_intermediate, "Intermediate"),
  STAR_min      = make_ci_plot_STAR_min(
    model_star_min,
    df_design = df_samp,
    init_state = initial_state,
    horizon = horizon_ci,
    delta = delta,
    N_draws = N_draws,
    forecast_dates = forecast_dates_ci
  )
)

gridExtra::grid.arrange(plots$Minimal, plots$Intermediate, plots$STAR_min, ncol = 1)







##############################################################################
# SIGNIFICANCE TEST: PDF & DISTRIBUTION OF MAX EIGENVALUES
# BOOTSTRAP DISTRIBUTION OF λ_max AT ZERO EQUILIBRIUM
##############################################################################
# Residual bootstrap with simulation

B    <- (1000)
Tlen <- nrow(df_samp)
mean_h_reference <- mean(df_samp$h, na.rm=TRUE)
stopifnot(Tlen > 100)

# residual pools for each fitted model
res_lin <- residuals(model_linear)
res_min <- residuals(model_minimal)
res_int <- residuals(model_intermediate)

# storage
max_eig_lin <- max_eig_min <- max_eig_int <- numeric(B)

b <- 1
while (b <= B) {
  # --- Linear
  h_star <- try(simulate_nonlinear_sample(coef(model_linear), T=Tlen,
                                          delta=delta, N=40, resid_pool=res_lin,
                                          model="Linear", burnin=300, cap=1e6),
                silent=TRUE)
  if (inherits(h_star, "try-error")) next
  
  df_star <- make_features_from_h(h_star, delta, N=40)
  fit_lin <- lm(h ~ x_lag1 + x_lag2 + X_lag1, data=df_star)
  max_eig_lin[b] <- rho_ss(coef(fit_lin), delta, "Linear", mean_h_reference)
  
  # --- Minimal
  h_star <- try(simulate_nonlinear_sample(coef(model_minimal), T=Tlen,
                                          delta=delta, N=40, resid_pool=res_min,
                                          model="Minimal", burnin=300, cap=1e6),
                silent=TRUE)
  if (inherits(h_star, "try-error")) next
  
  df_star <- make_features_from_h(h_star, delta, N=40)
  fit_min <- lm(h ~ x_lag1 + x_lag2 + X_lag1 + x_lag1_cub, data=df_star)
  max_eig_min[b] <- rho_ss(coef(fit_min), delta, "Minimal", mean_h_reference)
  
  # --- Intermediate
  h_star <- try(simulate_nonlinear_sample(coef(model_intermediate), T=Tlen,
                                          delta=delta, N=40, resid_pool=res_int,
                                          model="Intermediate", burnin=300, cap=1e6),
                silent=TRUE)
  if (inherits(h_star, "try-error")) next
  
  df_star <- make_features_from_h(h_star, delta, N=40)
  fit_int <- lm(h ~ x_lag1 + x_lag2 + X_lag1 +
                  x_lag1_cub + x_lag2_cub + X_lag1_cub +
                  x_lag1_sq_X_lag1 + x_lag1_X_lag1_sq,
                data=df_star)
  max_eig_int[b] <- rho_ss(coef(fit_int), delta, "Intermediate", mean_h_reference)
  
  b <- b + 1
}

df_eig <- data.frame(
  evalue = c(max_eig_lin, max_eig_min, max_eig_int),
  Model  = as_model_factor(rep(c("Linear","Minimal","Intermediate"), each = B))
)

ggplot(df_eig, aes(x = evalue, fill = Model)) +
  geom_density(alpha = 0.5, na.rm = TRUE) +
  labs(title = "Bootstrap distribution of |λ|max at steady state (residual bootstrap)",
       x = expression("|λ|max"), y = "Density", fill = "Model") +
  coord_cartesian(xlim = c(0.9, 1.2)) +
  scale_fill_manual(values = model_pal, drop = TRUE) +
  theme_minimal() + theme(legend.position = "top")







##############################################################################
# P‑VALUE TRUE EIGENVALUES >1 UNDER STABLE AR(2) DGP
##############################################################################
B_param <- 2000
B_data  <- 1000
Tlen    <- nrow(df_samp)
burnin  <- 500
mean_h_reference <- mean(df_samp$h, na.rm=TRUE)

# Normalize STAR_min linear parameters to the names rho_ss_star_min() expects
extract_star_min_linear_beta <- function(star_fit, warn = TRUE) {
  cf <- star_fit$coef
  nm <- names(cf)
  
  pick <- function(cands) {
    hit <- intersect(cands, nm)
    if (length(hit)) unname(cf[hit[1]]) else 0
  }
  
  missing_any <- function(cands) !any(cands %in% nm)
  
  out <- c(
    a0  = pick(c("a0","const","(Intercept)")),
    b10 = pick(c("b10","b1","phi1")),
    d1  = pick(c("d1","g1","gamma1")),
    b20 = pick(c("b20","b2","phi2")),
    d2  = pick(c("d2","g2","gamma2")),
    b30 = pick(c("b30","b3","phi3","X_lag1")),
    d3  = pick(c("d3","g3","gamma3"))
  )
  
  if (warn) {
    miss <- c(
      if (missing_any(c("a0","const","(Intercept)"))) "a0",
      if (missing_any(c("b10","b1","phi1"))) "b10",
      if (missing_any(c("d1","g1","gamma1"))) "d1",
      if (missing_any(c("b20","b2","phi2"))) "b20",
      if (missing_any(c("d2","g2","gamma2"))) "d2",
      if (missing_any(c("b30","b3","phi3","X_lag1"))) "b30",
      if (missing_any(c("d3","g3","gamma3"))) "d3"
    )
    if (length(miss)) message("STAR_min: filled 0 for missing linear coef(s): ", paste(miss, collapse=", "))
  }
  out
}

# 1) AR(2) fit to the cycle
ar2_fit <- lm(h ~ x_lag1 + x_lag2, data=df_samp)
u_hat   <- residuals(ar2_fit)
co_hat  <- coef(ar2_fit)

sim_ar2 <- function(beta, Tlen, u_pool, burnin=500, cap=1e6) {
  b0 <- beta["(Intercept)"] %||% 0
  b1 <- beta["x_lag1"] %||% 0
  b2 <- beta["x_lag2"] %||% 0
  r <- polyroot(c(1, -b1, -b2)); if (any(Mod(r) <= 1)) stop("nonstationary")
  e_star <- sample(u_pool - mean(u_pool), size=Tlen+burnin, replace=TRUE)
  x <- numeric(Tlen+burnin); x[1:2] <- 0
  for (t in 3:(Tlen+burnin)) {
    x[t] <- b0 + b1*x[t-1] + b2*x[t-2] + e_star[t]
    if (!is.finite(x[t]) || abs(x[t]) > cap) stop("explosive")
  }
  x[(burnin+1):(Tlen+burnin)]
}

# 2) Parameter bootstrap under AR(2) null
param_draws <- matrix(NA_real_, nrow=B_param, ncol=3,
                      dimnames=list(NULL, c("(Intercept)", "x_lag1", "x_lag2")))
i <- 1
while (i <= B_param) {
  x_star <- try(sim_ar2(co_hat, Tlen, u_hat, burnin), silent=TRUE)
  if (inherits(x_star,"try-error")) next
  df_tmp <- data.frame(h=x_star, x_lag1=dplyr::lag(x_star,1), x_lag2=dplyr::lag(x_star,2))
  df_tmp <- tidyr::drop_na(df_tmp)
  fit    <- lm(h ~ x_lag1 + x_lag2, data=df_tmp)
  beta_i <- coef(fit)
  rr <- polyroot(c(1, -beta_i["x_lag1"], -beta_i["x_lag2"]))
  if (any(Mod(rr) <= 1)) next
  param_draws[i,] <- beta_i[c("(Intercept)", "x_lag1", "x_lag2")]
  i <- i + 1
}

# 3) Observed |lambda|max for each model (include STAR_min)
obs_min  <- rho_ss(coef(model_minimal),      delta, "Minimal",      mean_h_reference)
obs_int  <- rho_ss(coef(model_intermediate), delta, "Intermediate", mean_h_reference)
obs_full <- rho_ss(coef(model_full),         delta, "Intermediate", mean_h_reference)

# STAR_min observed (use baseline beta from the STAR fit)
th <- model_star_min$coef
beta_star_obs <- extract_star_min_linear_beta(model_star_min)
obs_star <- rho_ss_star_min(beta_star_obs, model_star_min, delta, mean_h_reference)


# 4) Under AR(2) null: datasets -> fit models -> |lambda|max
boot_min  <- numeric(B_data)
boot_int  <- numeric(B_data)
boot_full <- numeric(B_data)
boot_star <- numeric(B_data)

set.seed(321)
b <- 1
while (b <= B_data) {
  beta <- param_draws[sample.int(nrow(param_draws), 1), ]
  x_star <- try(sim_ar2(beta, Tlen, u_hat, burnin), silent=TRUE)
  if (inherits(x_star,"try-error")) next
  
  df_sim <- make_features_from_h(x_star, delta, N=40)
  
  # polynomial models
  m_min  <- try(lm(h ~ x_lag1 + x_lag2 + X_lag1 + x_lag1_cub, data=df_sim), silent=TRUE)
  m_int  <- try(lm(h ~ x_lag1 + x_lag2 + X_lag1 +
                     x_lag1_cub + x_lag2_cub + X_lag1_cub +
                     x_lag1_sq_X_lag1 + x_lag1_X_lag1_sq, data=df_sim), silent=TRUE)
  m_full <- try(lm(h ~ x_lag1 + x_lag2 + X_lag1 +
                     x_lag1_sq + x_lag2_sq + X_lag1_sq +
                     x_lag1_x_lag2 + x_lag1_X_lag1 + x_lag2_X_lag1 +
                     x_lag1_cub + x_lag2_cub + X_lag1_cub +
                     x_lag1_sq_x_lag2 + x_lag1_x_lag2_sq +
                     x_lag1_sq_X_lag1 + x_lag2_sq_X_lag1 +
                     x_lag1_X_lag1_sq + x_lag2_X_lag1_sq +
                     x_lag1_x_lag2_X_lag1, data=df_sim), silent=TRUE)
  if (inherits(m_min,"try-error") || inherits(m_int,"try-error") || inherits(m_full,"try-error")) next
  
  boot_min[b]  <- rho_ss(coef(m_min),  delta, "Minimal",      mean_h_reference)
  boot_int[b]  <- rho_ss(coef(m_int),  delta, "Intermediate", mean_h_reference)
  boot_full[b] <- rho_ss(coef(m_full), delta, "Intermediate", mean_h_reference)
  
  # STAR_min on pseudo data (fast: eta,c,v fixed from baseline)
  ols_star <- try(ols_linear_part_star_min(df_sim, model_star_min), silent=TRUE)
  if (inherits(ols_star,"try-error")) next
  boot_star[b] <- rho_ss_star_min(ols_star$beta_hat, model_star_min, delta, mean_h_reference)
  
  b <- b + 1
}

# 5) p-values (STAR_min included)
p_min   <- mean(boot_min   >= obs_min,   na.rm=TRUE)
p_int   <- mean(boot_int   >= obs_int,   na.rm=TRUE)
p_full  <- mean(boot_full  >= obs_full,  na.rm=TRUE)
p_star  <- mean(boot_star  >= obs_star,  na.rm=TRUE)

cat("Minimal p-value (AR(2) null):     ", round(p_min,4),  "\n")
cat("Intermediate p-value (AR(2) null):", round(p_int,4),  "\n")
cat("Full p-value (AR(2) null):        ", round(p_full,4), "\n")
cat("STAR_min p-value (AR(2) null):    ", round(p_star,4), "\n")

# density plot including STAR_min
model_pal_ext <- c(model_pal, "STAR_min" = "#ff7f0e")
as_model_factor_ext <- function(x) factor(x, levels = names(model_pal_ext))

df_boot <- dplyr::tibble(
  evalue = c(boot_min, boot_int, boot_full, boot_star),
  Model  = as_model_factor_ext(rep(c("Minimal","Intermediate","Full","STAR_min"), each = B_data))
)

ggplot(df_boot, aes(x = evalue, fill = Model)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(0.85, 1.15)) +
  labs(title = "Distribution of |λ|max under AR(2) null (parameter+innovation bootstrap)",
       x = expression("|λ|max"), y = "Density", fill = "Model") +
  scale_fill_manual(values = model_pal_ext, drop = TRUE) +
  theme_minimal() + theme(legend.position = "top")






















##############################################################################
# CHAOS FALSE-POSITIVE TEST
##############################################################################

# A) Simulate Day’s Chaotic Map
# Chaotic Logistic Map:
# The map x_{t+1} = r * x_t * (1 - x_t) (with 3.57 < r ≤ 4) is the canonical
# logistic map exhibiting deterministic chaos.  It admits two rigorous micro‐
# founded readings in economics:
#  Day’s “Irregular Growth Cycles” Model (One‐Sector Growth)

# 1. same date vector
dates   <- data_merged_real$date
n_obs   <- length(dates)
N <= 40
# 2. Simulate Logistic‐cobweb (Day’s “Pollution” map)
simulate_day <- function(n, B = 3.9, m = 1, y = 1, a = 1, X = 0, k0 = 0.2) {
  k <- numeric(n)
  k[1] <- k0
  for(t in 2:n) {
    k[t] <- (a * B * k[t-1] * (m - k[t-1])) / (1 + X)
  }
  k
}
h_sim <- simulate_day(n = n_obs, B = 3.9, m = 1, y = 1, a = 1, X = 0, k0 = 0.2)

# 3. Attach to data frame
data_merged_real <- data_merged_real %>%
  mutate(
    h      = h_sim,      # chaotic “cyclical” component
    date   = dates
  )

ggplot(data_merged_real, aes(x = date, y = h)) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Simulated Logistic–Cobweb (Day’s “Pollution” Map)",
    x = "Date",
    y = "h (chaotic component)"
  ) +
  theme_minimal()


##  utilities
`%||%` <- function(a,b) if (is.null(a) || length(a)==0L || (length(a)==1L && is.na(a))) b else a

infer_opy <- function(dates) {
  # infer observations-per-year from the median spacing in days
  dt <- median(as.numeric(diff(as.Date(dates))), na.rm = TRUE)
  if (!is.finite(dt) || dt <= 0) return(4L)
  as.integer(round(365.25 / dt))
}

# truncated EWMA accumulator X_t = sum_{k=0}^{N-1} delta*(1-delta)^k * h_{t-k}
calc_accumulation_truncated <- function(h, delta, N = 40L) {
  n <- length(h); out <- rep(NA_real_, n); w <- delta * (1 - delta)^(0:(N-1))
  for (t in seq_len(n)) {
    k <- min(t, N)
    out[t] <- sum(w[1:k] * rev(h[(t-k+1):t]))
  }
  out
}

# Feature builder 
.add_features_local <- function(df) {
  transform(df,
            x_lag1_sq = x_lag1^2, x_lag2_sq = x_lag2^2, X_lag1_sq = X_lag1^2,
            x_lag1_cub = x_lag1^3, x_lag2_cub = x_lag2^3, X_lag1_cub = X_lag1^3,
            x_lag1_x_lag2 = x_lag1*x_lag2, x_lag1_X_lag1 = x_lag1*X_lag1, x_lag2_X_lag1 = x_lag2*X_lag1,
            x_lag1_sq_x_lag2 = (x_lag1^2)*x_lag2, x_lag1_x_lag2_sq = x_lag1*(x_lag2^2),
            x_lag1_sq_X_lag1 = (x_lag1^2)*X_lag1, x_lag2_sq_X_lag1 = (x_lag2^2)*X_lag1,
            x_lag1_X_lag1_sq = x_lag1*(X_lag1^2), x_lag2_X_lag1_sq = x_lag2*(X_lag1^2),
            x_lag1_x_lag2_X_lag1 = x_lag1*x_lag2*X_lag1
  )
}
add_features_safe <- function(df) if (exists("add_features", mode="function")) add_features(df) else .add_features_local(df)

# Core recursive simulator (shared by all models)
simulate_core_recursive <- function(h_tm1, h_tm2, H_tm1, horizon, delta, law_h) {
  h_out <- H_out <- numeric(horizon)
  for (t in 1:horizon) {
    h_t <- law_h(h_tm1, h_tm2, H_tm1)
    H_t <- (1 - delta) * H_tm1 + delta * h_t
    h_out[t] <- h_t; H_out[t] <- H_t
    h_tm2 <- h_tm1; h_tm1 <- h_t; H_tm1 <- H_t
  }
  data.frame(Time = seq_len(horizon), h = h_out, H = H_out)
}

## STAR helpers (estimation + simulation) 
if (!exists("estimate_star", mode="function")) stop("estimate_star() not found.")
if (!exists("vcov_star",     mode="function")) stop("vcov_star() not found.")

is_star_fit <- function(x) inherits(x, "star_fit") ||
  (is.list(x) && all(c("spec","z_var","coef") %in% names(x)))

.G_logistic <- function(z, eta, c) 1/(1 + exp(-exp(eta) * (z - c)))
.G_estar    <- function(z, eta, c) 1 - exp(-exp(eta) * (z - c)^2)

# Build F(h_{t-1}, h_{t-2}, X_{t-1}) from a STAR fit
.build_F3 <- function(fit) {
  stopifnot(is_star_fit(fit))
  th <- fit$coef; spec <- fit$spec; zvar <- fit$z_var; v <- fit$v %||% c(1,0,0)
  function(x1,x2,X1){
    z <- if (identical(zvar,"x_lag1")) x1 else v[1]*x1 + v[2]*x2 + v[3]*X1
    if (spec == "standard_lstar") {
      a0 <- th["a0"]; a1 <- th["a1"]
      b10 <- th["b10"]; b11 <- th["b11"]
      b20 <- th["b20"]; b21 <- th["b21"]
      b30 <- th["b30"]; b31 <- th["b31"]
      eta <- th["eta"];  c  <- th["c"]
      G <- .G_logistic(z, eta, c)
      (a0 + a1*G) + (b10 + b11*G)*x1 + (b20 + b21*G)*x2 + (b30 + b31*G)*X1
    } else {
      a0 <- th["a0"]
      b10 <- th["b10"]; b20 <- th["b20"]; b30 <- th["b30"]
      d1  <- th["d1"];  d2  <- th["d2"];  d3  <- th["d3"]
      eta <- th["eta"]; c   <- th["c"]
      G <- if (spec == "minimal_lstar") .G_logistic(z, eta, c) else .G_estar(z, eta, c)
      a0 + (b10 + d1*G)*x1 + (b20 + d2*G)*x2 + (b30 + d3*G)*X1
    }
  }
}

simulate_forecast_star <- function(star_fit, h0, h_1, H_1, horizon, delta) {
  F3 <- .build_F3(star_fit)
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, function(h1,h2,H1) F3(h1,h2,H1))
}




# B) Build Accumulator X and All Lags / Nonlinear Terms

##  (B) Build accumulator X and features 
data_merged_real <- data_merged_real %>%
  arrange(date) %>%
  mutate(
    X      = calc_accumulation_truncated(h, delta, N),
    X_lag1 = dplyr::lag(X, 1),
    x_lag1 = dplyr::lag(h, 1),
    x_lag2 = dplyr::lag(h, 2)
  ) %>%
  slice(-(1:2)) %>%        # drop first two rows (NA lags)
  add_features_safe()

## (C) Estimate nested polynomial models 
model_ar2 <- lm(h ~ x_lag1 + x_lag2, data = data_merged_real)

model_linear <- lm(h ~ x_lag1 + x_lag2 + X_lag1, data = data_merged_real)

model_minimal <- lm(h ~ x_lag1 + x_lag2 + X_lag1 + x_lag1_cub, data = data_merged_real)

model_intermediate <- lm(
  h ~ x_lag1 + x_lag2 + X_lag1 +
    x_lag1_cub + x_lag2_cub + X_lag1_cub +
    x_lag1_sq_X_lag1 + x_lag1_X_lag1_sq,
  data = data_merged_real
)

model_full <- lm(
  h ~ x_lag1 + x_lag2 + X_lag1 +
    x_lag1_sq + x_lag2_sq + X_lag1_sq +
    x_lag1_x_lag2 + x_lag1_X_lag1 + x_lag2_X_lag1 +
    x_lag1_cub + x_lag2_cub + X_lag1_cub +
    x_lag1_sq_x_lag2 + x_lag1_x_lag2_sq +
    x_lag1_sq_X_lag1 + x_lag2_sq_X_lag1 +
    x_lag1_X_lag1_sq + x_lag2_X_lag1_sq +
    x_lag1_x_lag2_X_lag1,
  data = data_merged_real
)
model_full$coefficients[is.na(model_full$coefficients)] <- 0

## LASSO (CV) 
mf_full <- model.frame(model_full)
tt_full <- terms(model_full)
X_raw   <- model.matrix(tt_full, data = mf_full)[, -1, drop = FALSE]
y_full  <- model.response(mf_full)
stopifnot(nrow(X_raw) == length(y_full))

scaler  <- pmax(apply(X_raw, 2, sd), .Machine$double.eps)
X_full  <- sweep(X_raw, 2, scaler, "/")

set.seed(123)
cvfit     <- glmnet::cv.glmnet(x = X_full, y = y_full, alpha = 1, standardize = FALSE, nfolds = 10)
lambda_cv <- cvfit$lambda.min
beta_hat  <- as.matrix(stats::coef(cvfit, s = lambda_cv))

survivors <- setdiff(rownames(beta_hat)[beta_hat[,1] != 0], "(Intercept)")
if (length(survivors) == 0L) {
  model_lasso <- lm(h ~ 1, data = mf_full)
} else {
  lasso_form  <- reformulate(survivors, response = "h")
  model_lasso <- lm(lasso_form, data = mf_full)
}
model_lasso$coefficients[is.na(model_lasso$coefficients)] <- 0
coef_lasso_vec <- coef(model_lasso)

simulate_forecast_lasso <- function(beta, h0, h_1, H_1, horizon, delta) {
  z <- function(nm) ifelse(!is.na(beta[nm]), beta[nm], 0)
  law_h <- function(h1,h2,H1) {
    z("(Intercept)") +
      z("x_lag1")*h1 + z("x_lag2")*h2 + z("X_lag1")*H1 +
      z("x_lag1_sq")*h1^2 + z("x_lag2_sq")*h2^2 + z("X_lag1_sq")*H1^2 +
      z("x_lag1_x_lag2")*h1*h2 + z("x_lag1_X_lag1")*h1*H1 + z("x_lag2_X_lag1")*h2*H1 +
      z("x_lag1_cub")*h1^3 + z("x_lag2_cub")*h2^3 + z("X_lag1_cub")*H1^3 +
      z("x_lag1_sq_x_lag2")*(h1^2)*h2 + z("x_lag1_x_lag2_sq")*h1*(h2^2) +
      z("x_lag1_sq_X_lag1")*(h1^2)*H1 + z("x_lag2_sq_X_lag1")*(h2^2)*H1 +
      z("x_lag1_X_lag1_sq")*h1*(H1^2) + z("x_lag2_X_lag1_sq")*h2*(H1^2) +
      z("x_lag1_x_lag2_X_lag1")*h1*h2*H1
  }
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, law_h)
}

## STAR models: FULL estimation + SEs 
# Thorough search; SEs via vcov_star()
model_star_std  <- estimate_star(
  data_merged_real, spec = "standard_lstar", z_var = "x_lag1",
  starts_gamma = c(0.5,1,2,5,10), starts_c = c("median","mean"),
  trace = TRUE, maxit = 4000
)
model_star_min  <- estimate_star(
  data_merged_real, spec = "minimal_lstar",  z_var = "x_lag1",
  starts_gamma = c(0.5,1,2,5,10), starts_c = c("median","mean"),
  trace = TRUE, maxit = 4000
)
model_estar_min <- estimate_star(
  data_merged_real, spec = "minimal_estar",  z_var = "x_lag1",
  starts_gamma = c(0.5,1,2,5,10), starts_c = c("median","mean"),
  trace = TRUE, maxit = 4000
)

for (m in c("model_star_std","model_star_min","model_estar_min")) {
  fit <- get(m)
  V   <- vcov_star(fit, data_merged_real)
  fit$vcov <- V
  fit$se   <- sqrt(diag(V))
  assign(m, fit, inherits = TRUE)
}

## (D) Forecast setup 
# horizon & origin
horizon    <- get0("horizon") %||% as.integer(1000 * opy / 4)
inc_days   <- median(as.numeric(diff(dates)), na.rm = TRUE)   # days step
if (!is.finite(inc_days) || inc_days <= 0) inc_days <- 90
inc        <- inc_days                                       # used with Date arithmetic

N_needed   <- N
start_idx  <- N_needed + 3L            # need at least two lags + N history
start_date <- data_merged_real$date[start_idx]

# initial state from chaotic data just before origin
h_hist <- data_merged_real$h[1:(start_idx-1)]
H_hist <- calc_accumulation_truncated(h_hist, delta, N)
h0  <- h_hist[length(h_hist)]
h_1 <- h_hist[length(h_hist)-1L]
H_1 <- H_hist[length(H_hist)]


## (E) Deterministic forecasts for all models
# AR2 and polynomial simulators (explicit laws)
simulate_forecast_ar2 <- function(model, h0, h_1, horizon) {
  co <- coef(model); b0 <- co["(Intercept)"] %||% 0; b1 <- co["x_lag1"] %||% 0; b2 <- co["x_lag2"] %||% 0
  h_out <- numeric(horizon); h_tm1 <- h0; h_tm2 <- h_1
  for (t in 1:horizon) { h_t <- b0 + b1*h_tm1 + b2*h_tm2; h_out[t] <- h_t; h_tm2 <- h_tm1; h_tm1 <- h_t }
  data.frame(Time = seq_len(horizon), h = h_out)
}
simulate_forecast_linear <- function(model, h0, h_1, H_1, horizon, delta) {
  co <- coef(model); z <- function(nm) co[nm] %||% 0
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, function(h1,h2,H1) z("(Intercept)") + z("x_lag1")*h1 + z("x_lag2")*h2 + z("X_lag1")*H1)
}
simulate_forecast_minimal <- function(model, h0, h_1, H_1, horizon, delta) {
  co <- coef(model); z <- function(nm) co[nm] %||% 0
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, function(h1,h2,H1)
    z("(Intercept)") + z("x_lag1")*h1 + z("x_lag2")*h2 + z("X_lag1")*H1 + z("x_lag1_cub")*h1^3)
}
simulate_forecast_intermediate <- function(model, h0, h_1, H_1, horizon, delta) {
  co <- coef(model); z <- function(nm) co[nm] %||% 0
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, function(h1,h2,H1)
    z("(Intercept)") + z("x_lag1")*h1 + z("x_lag2")*h2 + z("X_lag1")*H1 +
      z("x_lag1_cub")*h1^3 + z("x_lag2_cub")*h2^3 + z("X_lag1_cub")*H1^3 +
      z("x_lag1_sq_X_lag1")*(h1^2)*H1 + z("x_lag1_X_lag1_sq")*h1*(H1^2))
}
simulate_forecast_full <- function(model, h0, h_1, H_1, horizon, delta) {
  co <- coef(model); z <- function(nm) co[nm] %||% 0
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, function(h1,h2,H1)
    z("(Intercept)") + z("x_lag1")*h1 + z("x_lag2")*h2 + z("X_lag1")*H1 +
      z("x_lag1_sq")*h1^2 + z("x_lag2_sq")*h2^2 + z("X_lag1_sq")*H1^2 +
      z("x_lag1_x_lag2")*h1*h2 + z("x_lag1_X_lag1")*h1*H1 + z("x_lag2_X_lag1")*h2*H1 +
      z("x_lag1_cub")*h1^3 + z("x_lag2_cub")*h2^3 + z("X_lag1_cub")*H1^3 +
      z("x_lag1_sq_x_lag2")*(h1^2)*h2 + z("x_lag1_x_lag2_sq")*h1*(h2^2) +
      z("x_lag1_sq_X_lag1")*(h1^2)*H1 + z("x_lag2_sq_X_lag1")*(h2^2)*H1 +
      z("x_lag1_X_lag1_sq")*h1*(H1^2) + z("x_lag2_X_lag1_sq")*h2*(H1^2) +
      z("x_lag1_x_lag2_X_lag1")*h1*h2*H1)
}

# Forecasts
fc_ar2        <- simulate_forecast_ar2         (model_ar2,        h0, h_1,              horizon)
fc_lin        <- simulate_forecast_linear      (model_linear,     h0, h_1, H_1, horizon, delta)
fc_min        <- simulate_forecast_minimal     (model_minimal,    h0, h_1, H_1, horizon, delta)
fc_int        <- simulate_forecast_intermediate(model_intermediate,h0, h_1, H_1, horizon, delta)
fc_full       <- simulate_forecast_full        (model_full,       h0, h_1, H_1, horizon, delta)
fc_las        <- simulate_forecast_lasso       (coef_lasso_vec,   h0, h_1, H_1, horizon, delta)

# STAR (all three)
fc_star_std   <- simulate_forecast_star(model_star_std,  h0, h_1, H_1, horizon, delta)
fc_star_min   <- simulate_forecast_star(model_star_min,  h0, h_1, H_1, horizon, delta)
fc_estar_min  <- simulate_forecast_star(model_estar_min, h0, h_1, H_1, horizon, delta)

# Collect
fc_list <- list(
  AR2           = fc_ar2,
  Linear        = fc_lin,
  Minimal       = fc_min,
  Intermediate  = fc_int,
  Full          = fc_full,
  LASSO         = fc_las,
  STAR        = fc_star_std,
  STAR_min      = fc_star_min,
  ESTAR        = fc_estar_min
)



## (F) Phase portraits
# Empirical chaotic segment start (one step after start_date)
empirical <- data_merged_real %>%
  dplyr::filter(date >= start_date + inc) %>%
  dplyr::transmute(t = dplyr::row_number(), h = h, h_lag = dplyr::lag(h)) %>%
  dplyr::filter(!is.na(h_lag)) %>%
  dplyr::mutate(model = "TrueChaos")

phase_dfs <- lapply(names(fc_list), function(nm) {
  fc <- fc_list[[nm]]
  tibble(
    t     = fc$Time,
    h     = fc$h,
    h_lag = dplyr::lag(fc$h),
    model = nm
  ) %>% dplyr::filter(!is.na(h_lag))
})

phase_all <- dplyr::bind_rows(empirical, dplyr::bind_rows(phase_dfs))

ggplot(phase_all, aes(x = h_lag, y = h, color = model)) +
  geom_path(alpha = 0.8) +
  facet_wrap(~ model, ncol = 2, scales = "free") +
  scale_color_manual(values = c(
    "TrueChaos"   = "black",
    "AR2"         = "#1f78b4",
    "Linear"      = "#fb9a99",
    "Minimal"     = "#33a02c",
    "Intermediate"= "#e31a1c",
    "Full"        = "#ff7f00",
    "LASSO"       = "#6a3d9a",
    "STAR"    = "#1b9e77",
    "STAR_min"    = "#d95f02",
    "ESTAR"   = "#7570b3"
  )) +
  labs(title = "",
       x = expression(h[t-1]), y = expression(h[t]), color = "Model") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "bottom")

## (G) Time-series small multiples 
true_start <- start_date + inc
true_end   <- max(data_merged_real$date, na.rm = TRUE)
fc_dates   <- seq(from = as.Date(true_start), by = inc_days, length.out = horizon)

ts_dfs <- lapply(names(fc_list), function(nm) {
  fc <- fc_list[[nm]]
  tibble(date = fc_dates, value = fc$h, model = nm) %>%
    dplyr::filter(date <= true_end)
})
ts_true <- data_merged_real %>%
  dplyr::filter(date >= true_start) %>%
  dplyr::transmute(date, value = h, model = "TrueChaos")

ts_all <- dplyr::bind_rows(ts_true, dplyr::bind_rows(ts_dfs)) %>%
  dplyr::mutate(panel = model)

ggplot(ts_all, aes(x = date, y = value, color = model)) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  facet_wrap(~ panel, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c(
    "TrueChaos"   = "black",
    "AR2"         = "#1f78b4",
    "Linear"      = "#fb9a99",
    "Minimal"     = "#33a02c",
    "Intermediate"= "#e31a1c",
    "Full"        = "#ff7f00",
    "LASSO"       = "#6a3d9a",
    "STAR"    = "#1b9e77",
    "STAR_min"    = "#d95f02",
    "ESTAR"   = "#7570b3"
  )) +
  scale_x_date(limits = c(as.Date(true_start), as.Date(true_end))) +
  labs(title = "",
       x = "Date", y = expression(h[t]), color = "Series") +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))



# 1) Accumulation series H_{t-1}
data_merged_real <- data_merged_real %>%
  arrange(date) %>%
  mutate(
    X      = calc_accumulation_truncated(h, delta, N),
    X_lag1 = lag(X),
    x_lag1 = lag(h, 1),
    x_lag2 = lag(h, 2)
  ) %>%
  slice(-c(1,2))  # drop the first two rows with NA lags

# 2) Nonlinear transformations
# Build all polynomial & interaction features with canonical names
data_merged_real <- add_features(data_merged_real)



# C) Estimate the Nested Models on the Chaotic Sample

# 1) AR(2)
model_ar2 <- lm(h ~ x_lag1 + x_lag2,
                data = data_merged_real)

# 2) Linear (adds X_lag1)
model_linear <- lm(h ~ x_lag1 + x_lag2 + X_lag1,
                   data = data_merged_real)

# 3) Minimal (adds x_lag1_cub)
model_minimal <- lm(h ~ x_lag1 + x_lag2 + X_lag1 + x_lag1_cub,
                    data = data_merged_real)

# 4) Intermediate
model_intermediate <- lm(
  h ~ x_lag1 + x_lag2 + X_lag1 +
    x_lag1_cub + x_lag2_cub + X_lag1_cub +
    x_lag1_sq_X_lag1 + x_lag1_X_lag1_sq,
  data = data_merged_real)

# 5) Full
model_full <- lm(
  h ~ x_lag1 + x_lag2 + X_lag1 +
    x_lag1_sq + x_lag2_sq + X_lag1_sq +
    x_lag1_x_lag2 + x_lag1_X_lag1 + x_lag2_X_lag1 +
    x_lag1_cub + x_lag2_cub + X_lag1_cub +
    x_lag1_sq_x_lag2 + x_lag1_x_lag2_sq +
    x_lag1_sq_X_lag1 + x_lag2_sq_X_lag1 +
    x_lag1_X_lag1_sq + x_lag2_X_lag1_sq +
    x_lag1_x_lag2_X_lag1,
  data = data_merged_real
)
# after fitting
model_full$coefficients[ is.na(model_full$coefficients) ] <- 0

# 6) LASSO (CV‐based)
## Build X and y from the *actual* data used by model_full ---
mf_full <- model.frame(model_full)                            # rows kept by lm (na.action applied)
tt_full <- terms(model_full)                                  # same terms as the fitted model
X_raw   <- model.matrix(tt_full, data = mf_full)[, -1, drop = FALSE]  # drop intercept
y_full  <- model.response(mf_full)

## sanity check
stopifnot(nrow(X_raw) == length(y_full))

## Column scaling (no zero-division) 
scaler <- pmax(apply(X_raw, 2, sd), .Machine$double.eps)
X_full <- sweep(X_raw, 2, scaler, "/")

## CV LASSO (note: use `s=` when extracting coefs) 
set.seed(123)  # for reproducibility
cvfit     <- cv.glmnet(x = X_full, y = y_full, alpha = 1, standardize = FALSE, nfolds = 10)
lambda_cv <- cvfit$lambda.min
beta_hat  <- as.matrix(coef(cvfit, s = lambda_cv))

## survivors & refit a clean OLS on those (optional but nice)
survivors  <- setdiff(rownames(beta_hat)[beta_hat[,1] != 0], "(Intercept)")
if (length(survivors) == 0L) {
  # fallback: intercept-only model
  model_lasso <- lm(h ~ 1, data = mf_full)
} else {
  lasso_form  <- reformulate(survivors, response = "h")
  model_lasso <- lm(lasso_form, data = mf_full)
}
model_lasso$coefficients[is.na(model_lasso$coefficients)] <- 0
coef_lasso_vec <- coef(model_lasso)
## LASSO simulator 
simulate_forecast_lasso <- function(beta, h0, h_1, H_1, horizon, delta) {
  z <- function(nm) ifelse(!is.na(beta[nm]), beta[nm], 0)
  law_h <- function(h1, h2, H1) {
    z("(Intercept)") +
      z("x_lag1")*h1 + z("x_lag2")*h2 + z("X_lag1")*H1 +
      z("x_lag1_sq")*h1^2 + z("x_lag2_sq")*h2^2 + z("X_lag1_sq")*H1^2 +
      z("x_lag1_x_lag2")*h1*h2 + z("x_lag1_X_lag1")*h1*H1 + z("x_lag2_X_lag1")*h2*H1 +
      z("x_lag1_cub")*h1^3 + z("x_lag2_cub")*h2^3 + z("X_lag1_cub")*H1^3 +
      z("x_lag1_sq_x_lag2")*(h1^2)*h2 + z("x_lag1_x_lag2_sq")*h1*(h2^2) +
      z("x_lag1_sq_X_lag1")*(h1^2)*H1 + z("x_lag2_sq_X_lag1")*(h2^2)*H1 +
      z("x_lag1_X_lag1_sq")*h1*(H1^2) + z("x_lag2_X_lag1_sq")*h2*(H1^2) +
      z("x_lag1_x_lag2_X_lag1")*h1*h2*H1
  }
  simulate_core_recursive(h0, h_1, H_1, horizon, delta, law_h)
}


# STAR fits (fast) 
# Guards
`%||%` <- function(a,b) if (is.null(a) || length(a)==0L || (length(a)==1L && is.na(a))) b else a
if (!exists("is_star_fit", mode = "function")) {
  is_star_fit <- function(x) inherits(x, "star_fit") ||
    (is.list(x) && all(c("spec","z_var","coef") %in% names(x)))
}
# Minimal STAR builder if .build_F3 isn't already in session
if (!exists(".build_F3", mode="function")) {
  .G_logistic <- function(z, eta, c) 1/(1 + exp(-exp(eta) * (z - c)))
  .G_estar    <- function(z, eta, c) 1 - exp(-exp(eta) * (z - c)^2)
  .build_F3 <- function(fit) {
    stopifnot(is_star_fit(fit))
    th   <- fit$coef
    spec <- fit$spec
    zvar <- fit$z_var
    v    <- fit$v %||% c(1,0,0)
    function(x1,x2,X1) {
      z <- if (identical(zvar, "x_lag1")) x1 else v[1]*x1 + v[2]*x2 + v[3]*X1
      if (spec == "standard_lstar") {
        a0 <- th["a0"]; a1 <- th["a1"]
        b10 <- th["b10"]; b11 <- th["b11"]
        b20 <- th["b20"]; b21 <- th["b21"]
        b30 <- th["b30"]; b31 <- th["b31"]
        eta <- th["eta"];  c  <- th["c"]
        G <- .G_logistic(z, eta, c)
        (a0 + a1*G) + (b10 + b11*G)*x1 + (b20 + b21*G)*x2 + (b30 + b31*G)*X1
      } else {
        # minimal_lstar or minimal_estar
        a0 <- th["a0"]
        b10 <- th["b10"]; b20 <- th["b20"]; b30 <- th["b30"]
        d1  <- th["d1"];  d2  <- th["d2"];  d3  <- th["d3"]
        eta <- th["eta"]; c   <- th["c"]
        G <- if (spec == "minimal_lstar") .G_logistic(z, eta, c) else .G_estar(z, eta, c)
        a0 + (b10 + d1*G)*x1 + (b20 + d2*G)*x2 + (b30 + d3*G)*X1
      }
    }
  }
}

# Quick STAR estimator calls (no SEs, tight starts) on the chaotic sample
stopifnot(exists("estimate_star"))
model_star_std  <- estimate_star(
  data_merged_real, spec = "standard_lstar", z_var = "x_lag1",
  starts_gamma = c(1,5), starts_c = "median", trace = FALSE, maxit = 800
)
model_star_min  <- estimate_star(
  data_merged_real, spec = "minimal_lstar",  z_var = "x_lag1",
  starts_gamma = c(1,5), starts_c = "median", trace = FALSE, maxit = 800
)
model_estar_min <- estimate_star(
  data_merged_real, spec = "minimal_estar",  z_var = "x_lag1",
  starts_gamma = c(1,5), starts_c = "median", trace = FALSE, maxit = 800
)

# Deterministic simulator for any STAR fit 
simulate_forecast_star <- function(star_fit, h0, h_1, H_1, horizon, delta) {
  F3 <- .build_F3(star_fit)
  simulate_core_recursive(
    h_tm1 = h0, h_tm2 = h_1, H_tm1 = H_1,
    horizon = horizon, delta = delta,
    law_h = function(h1,h2,H1) F3(h1,h2,H1)
  )
}




##### DETERMINISTIC SIMULATIONS (CHAOS)

# 0) Setup: horizon & starting state from chaotic data (N quarters + 2 lags of history)

horizon    <- 1000 * obs_per_year_real / 4     # same as before
N_needed   <- N                    # e.g. 10*obs_per_year_real = 40 quarters
warmup_obs <- max(N_needed, 2)     # also need 2 obs for x_lag2
# (b) index into date vector
#     take the (warmup_obs + 1)-th observation as the forecast origin
start_idx  <- warmup_obs + 1
# (c) pull the actual date
start_date <- data_merged_real$date[start_idx]
# (d) the increment between observations
inc        <- if(freq_real=="q") months(3) else if(freq_real=="m") months(1) else weeks(1)
# (e) confirm
cat("Forecasting from ", start_date,
    " (obs #", start_idx, 
    ") with ", warmup_obs, " prior obs available\n")
# build initial state (h_{t-1}, h_{t-2}, H_{t-1}) from the chaotic data
init_objs <- build_initial_objects(
  data_merged_real, start_date)
h0    <- init_objs$initial_state["h"]
h_1   <- init_objs$initial_state["h_lag"]
H_1   <- init_objs$H_start


# 1) Deterministic forecast for each model
fc_ar2  <- simulate_forecast_ar2         (model_ar2,  h0, h_1, horizon)
fc_lin  <- simulate_forecast_linear      (model_linear,  h0, h_1, H_1, horizon, delta)
fc_min  <- simulate_forecast_minimal     (model_minimal, h0, h_1, H_1, horizon, delta)
fc_int  <- simulate_forecast_intermediate(model_intermediate, h0, h_1, H_1, horizon, delta)
fc_full <- simulate_forecast_full        (model_full, h0, h_1, H_1, horizon, delta)
fc_las  <- simulate_forecast_lasso       (coef_lasso_vec, h0, h_1, H_1, horizon, delta)

head(fc_ar2)
head(fc_lin)
head(fc_min)
head(fc_int)
head(fc_full)
head(fc_las)

# STAR forecasts
fc_star_std  <- simulate_forecast_star(model_star_std,  h0, h_1, H_1, horizon, delta)
fc_star_min  <- simulate_forecast_star(model_star_min,  h0, h_1, H_1, horizon, delta)
fc_estar_min <- simulate_forecast_star(model_estar_min, h0, h_1, H_1, horizon, delta)

# Optionally peek
head(fc_star_std); head(fc_star_min); head(fc_estar_min)

# Extend the list
fc_list <- list(
  AR2           = fc_ar2,
  Linear        = fc_lin,
  Minimal       = fc_min,
  Intermediate  = fc_int,
  Full          = fc_full,
  LASSO         = fc_las,
  STAR_std      = fc_star_std,
  STAR_min      = fc_star_min,
  ESTAR_min     = fc_estar_min
)

# palette addition
star_cols <- c(
  "STAR"  = "#1b9e77",
  "STAR_min"  = "#d95f02",
  "ESTAR_min" = "#7570b3"
)

# for phase-portrait plot
scale_color_manual(
  values = c(
    "TrueChaos"   = "black",
    "AR2"         = "#1f78b4",
    "Linear"      = "#fb9a99",
    "Minimal"     = "#33a02c",
    "Intermediate"= "#e31a1c",
    "Full"        = "#ff7f00",
    "LASSO"       = "#6a3d9a",
    star_cols
  )
)




# 2) Phase–portrait: s_{t-1} → s_t for each model + true chaos
# first build a data‐frame of the true chaotic path (shortly after start_date)
empirical <- data_merged_real %>%
  dplyr::filter(date >= start_date + inc) %>%   # explicitly from dplyr
  dplyr::transmute(
    t     = row_number(),
    h     = h,
    h_lag = dplyr::lag(h)
  ) %>%
  dplyr::filter(!is.na(h_lag)) %>%              # explicitly from dplyr
  dplyr::mutate(model = "TrueChaos")

# now build the forecast phase‐data
phase_dfs <- lapply(names(fc_list), function(nm) {
  df <- fc_list[[nm]]
  df2 <- df %>%
    dplyr::transmute(
      t     = Time,
      h     = h,
      h_lag = dplyr::lag(h)
    ) %>%
    dplyr::filter(!is.na(h_lag)) %>%
    dplyr::mutate(model = nm)
  df2
})

phase_all <- bind_rows(
  empirical,      # already built with dplyr::filter above
  bind_rows(phase_dfs)
)

ggplot(phase_all, aes(x = h_lag, y = h, color = model)) +
  geom_path(alpha = 0.8) +
  facet_wrap(~ model, ncol = 2, scales = "free") +
  scale_color_manual(
    values = c(
      "TrueChaos"    = "black",
      "AR2"           = "#1f78b4",
      "Linear"        = "#fb9a99",
      "Minimal"       = "#33a02c",
      "Intermediate"  = "#e31a1c",
      "Full"          = "#ff7f00",
      "LASSO"         = "#6a3d9a"
    )
  ) +
  labs(
    title = "",
    x     = expression(h[t-1]),
    y     = expression(h[t]),
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    strip.text      = element_text(face = "bold"),
    panel.spacing   = unit(0.5, "lines"),
    legend.position = "bottom"
  )



# 3) Small‑multiples of time‑series: TrueChaos + each model, restricted to
#    the empirical chaotic sample span

# determine the date range of the true chaotic series
true_start <- start_date + inc
true_end   <- max(data_merged_real$date)

# build forecast dates
fc_dates <- seq(
  from   = true_start,
  by     = obs_unit_real,
  length = horizon
)

# convert each forecast to a long tibble, then restrict to true_end
ts_dfs <- lapply(names(fc_list), function(nm) {
  fc <- fc_list[[nm]]
  tibble(
    date  = fc_dates,
    value = fc$h,
    model = nm
  ) %>%
    dplyr::filter(date <= true_end)
})

# true chaotic series
ts_true <- data_merged_real %>%
  dplyr::filter(date >= true_start) %>%
  dplyr::transmute(date, value = h, model = "TrueChaos")

# combine all
ts_all <- bind_rows(ts_true, bind_rows(ts_dfs)) %>%
  dplyr::mutate(panel = model)

# plot with free y‑scales and x limited to [true_start, true_end]
ggplot(ts_all, aes(x = date, y = value, color = model)) +
  geom_line(size = 0.7, alpha = 0.9) +
  facet_wrap(~panel, ncol = 2, scales = "free_y") +
  scale_color_manual(
    values = c(
      "TrueChaos"    = "black",
      "AR2"           = "#1f78b4",
      "Linear"        = "#fb9a99",
      "Minimal"       = "#33a02c",
      "Intermediate"  = "#e31a1c",
      "Full"          = "#ff7f00",
      "LASSO"         = "#6a3d9a"
    )
  ) +
  scale_x_date(limits = c(true_start, true_end)) +
  labs(
    title = "",
    x     = "Date",
    y     = expression(h[t]),
    color = "Series"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )






