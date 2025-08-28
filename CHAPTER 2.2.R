################################################################################
# CHAPTER 2 - BIVARIATE ANALYSIS Real + Financial (Section 2)
################################################################################

# Camille Souffron - MASTER THESIS APE (PSE & ENS)

#NB: all packages in the DATA_LOADING file


# 0) Global frequency (target)
LOG_h <- 100 # 100 if in log deviations %
LOG_s <- 1
freq_real <- "q"                            # "q","m","w"
obs_per_year_real <- switch(freq_real, q=4, m=12, w=52)
obs_unit_real     <- switch(freq_real, q="quarter", m="month", w="week")


# 1) Helpers: frequency + dates
# Floor to start-of-period and aggregate to target frequency.
coerce_to_target_freq <- function(df, date_col, value_col, target_unit = obs_unit_real,
                                  how = c("first","mean","sum")) {
  how <- match.arg(how)
  stopifnot(date_col %in% names(df), value_col %in% names(df))
  df <- df %>%
    dplyr::mutate(.date_raw = .data[[date_col]],
                  .date     = lubridate::floor_date(.data[[date_col]], unit = target_unit)) %>%
    dplyr::arrange(.date, .date_raw)
  
  out <- switch(
    how,
    first = df %>% dplyr::group_by(.date) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup(),
    mean  = df %>% dplyr::group_by(.date) %>%
      dplyr::summarise(!!value_col := mean(.data[[value_col]], na.rm=TRUE), .groups="drop"),
    sum   = df %>% dplyr::group_by(.date) %>%
      dplyr::summarise(!!value_col := sum(.data[[value_col]], na.rm=TRUE),  .groups="drop")
  )
  out %>%
    dplyr::transmute(date = .date, !!value_col := .data[[value_col]]) %>%
    dplyr::distinct(date, .keep_all = TRUE)
}

.as_num <- function(x) as.numeric(x)


# 2) Unit-root tests 
run_unit_root_tests <- function(x, model = c("level","trend")) {
  model <- match.arg(model)
  x_num <- stats::na.omit(as.numeric(x))
  if (length(x_num) < 20) stop("Series too short for unit-root tests.")
  adf_type <- if (model=="level") "drift" else "trend"
  adf_res  <- urca::ur.df(x_num, type=adf_type, selectlags="AIC")
  pp_model <- if (model=="level") "constant" else "trend"
  pp_res   <- urca::ur.pp(x_num, type="Z-alpha", model=pp_model, lags="short")
  kpss_null<- if (model=="level") "Level" else "Trend"
  kpss_res <- tseries::kpss.test(x_num, null=kpss_null)
  dfgls_res<- urca::ur.ers(x_num, type="DF-GLS", model=pp_model,
                           lag.max = trunc(4 * (length(x_num)/100)^(1/4)))
  list(ADF=summary(adf_res), PP=summary(pp_res), KPSS=kpss_res, DFGLS=summary(dfgls_res))
}


# 3) Detrending filters 
ideal_highpass_filter <- function(x, cutoff_freq = 1 / (20 * obs_per_year_real)) {
  N <- length(x); X <- stats::fft(x); k <- 0:(N-1)
  freq_k <- ifelse(k <= N/2, k/N, (k-N)/N)
  keep <- abs(freq_k) >= cutoff_freq
  X[!keep] <- 0
  Re(stats::fft(X, inverse=TRUE) / N)
}

.require_Rssa <- function() {
  if (!requireNamespace("Rssa", quietly = TRUE))
    stop("Package 'Rssa' is required for SSA. Do: install.packages('Rssa')")
}

.ssa_full_realtime_levels <- function(x, L, max_K, trend_min_period_obs, trend_lowshare_min,
                                      pair_wcor_min, pair_centroid_tol) {
  if (!exists(".ssa_full_realtime_levels_user", mode="function")) {
    stop("Define `.ssa_full_realtime_levels_user(...)` or use method='ssa_full'.")
  }
  .ssa_full_realtime_levels_user(
    x, L=L, max_K=max_K,
    trend_min_period_obs=trend_min_period_obs,
    trend_lowshare_min=trend_lowshare_min,
    pair_wcor_min=pair_wcor_min,
    pair_centroid_tol=pair_centroid_tol
  )
}

get_cycle_trend <- function(x, method) {
  stopifnot(is.numeric(x), length(x) > 5L)
  as_cycle_with_residual_trend <- function(cyc, x)
    list(cycle = as.numeric(cyc), trend = as.numeric(x - cyc), meta=list(source="residual_trend"))
  switch(method,
         ideal20y = { cutoff <- 1 / (20 * obs_per_year_real); cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         ideal100y= { cutoff <- 1 / (100* obs_per_year_real); cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         ideal90y = { cutoff <- 1 / (90 * obs_per_year_real); cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         ideal70y = { cutoff <- 1 / (70 * obs_per_year_real); cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         ideal60y = { cutoff <- 1 / (60 * obs_per_year_real); cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         ideal50y = { cutoff <- 1 / (50 * obs_per_year_real); cyc <- ideal_highpass_filter(x, cutoff)
         as_cycle_with_residual_trend(cyc, x) },
         bandpass = { bf <- signal::butter(2, bp_bounds, type="pass")
         cyc <- as.numeric(signal::filtfilt(bf, as.numeric(x)))
         as_cycle_with_residual_trend(cyc, x) },
         hp1600   = { out <- mFilter::hpfilter(stats::ts(x, freq = obs_per_year_real), freq=1600, type="lambda")
         list(cycle=as.numeric(out$cycle), trend=as.numeric(out$trend), meta=list(source="hpfilter")) },
         poly2    = { t <- seq_along(x); cyc <- stats::resid(stats::lm(x ~ t + I(t^2)))
         as_cycle_with_residual_trend(cyc, x) },
         poly3    = { t <- seq_along(x); cyc <- stats::resid(stats::lm(x ~ t + I(t^2)+I(t^3)))
         as_cycle_with_residual_trend(cyc, x) },
         poly4    = { t <- seq_along(x); cyc <- stats::resid(stats::lm(x ~ t + I(t^2)+I(t^3)+I(t^4)))
         as_cycle_with_residual_trend(cyc, x) },
         poly5    = { t <- seq_along(x); cyc <- stats::resid(stats::lm(x ~ t + I(t^2)+I(t^3)+I(t^4)+I(t^5)))
         as_cycle_with_residual_trend(cyc, x) },
         linear   = { t <- seq_along(x); cyc <- stats::resid(stats::lm(x ~ t))
         as_cycle_with_residual_trend(cyc, x) },
         log      = { t <- seq_along(x); cyc <- stats::resid(stats::lm(x ~ log(t)))
         as_cycle_with_residual_trend(cyc, x) },
         ssa_full = { .require_Rssa()
           dec <- .ssa_full_decompose_levels(
             x, L = cfg$ssa_full_L, max_K = cfg$ssa_full_max_components,
             trend_min_period_obs = cfg$ssa_trend_min_period_obs,
             trend_lowshare_min   = cfg$ssa_trend_lowfreq_share_min,
             pair_wcor_min        = cfg$ssa_pair_wcor_min,
             pair_centroid_tol    = cfg$ssa_pair_centroid_tol
           )
           list(cycle=as.numeric(dec$cycle), trend=as.numeric(dec$trend), meta=list(source="ssa_full")) },
         ssa_full_rt = { .require_Rssa()
           dec <- .ssa_full_realtime_levels(
             x, L = cfg$ssa_full_L, max_K = cfg$ssa_full_max_components,
             trend_min_period_obs = cfg$ssa_trend_min_period_obs,
             trend_lowshare_min   = cfg$ssa_trend_lowfreq_share_min,
             pair_wcor_min        = cfg$ssa_pair_wcor_min,
             pair_centroid_tol    = cfg$ssa_pair_centroid_tol
           )
           list(cycle=as.numeric(dec$cycle), trend=as.numeric(dec$trend), meta=list(source="ssa_full_rt")) },
         none     = { list(cycle = as.numeric(x), trend = rep(0, length(x)), meta=list(source="none")) },
         stop("Unknown detrending method: ", method)
  )
}

detrend_method <- "ideal20y"
yrs_low  <- 1.5   # 6Q
yrs_high <- 8     # 32Q
bp_bounds <- c(1/(yrs_high*obs_per_year_real), 1/(yrs_low*obs_per_year_real)) * 2


# 4) Harmonize Real + Financial
value_col_real <- "log_hours"
value_col_fin  <- "finance"

real_aligned <- coerce_to_target_freq(data_merged_real, "date", value_col_real,
                                      target_unit = obs_unit_real, how = "first")
fin_aligned  <- coerce_to_target_freq(data_merged_fin,  "date", value_col_fin,
                                      target_unit = obs_unit_real, how = "first")

common_dates <- intersect(real_aligned$date, fin_aligned$date)
real_aligned <- real_aligned %>% dplyr::filter(date %in% common_dates) %>% dplyr::arrange(date)
fin_aligned  <- fin_aligned  %>% dplyr::filter(date %in% common_dates) %>% dplyr::arrange(date)

stopifnot(nrow(real_aligned) == nrow(fin_aligned),
          all(real_aligned$date == fin_aligned$date))


# 5) Cycles + trends for both
res_ct_real <- get_cycle_trend(real_aligned[[value_col_real]], detrend_method)
data_merged_real <- real_aligned %>%
  dplyr::mutate(
    h          = LOG_h * res_ct_real$cycle,
    trend_real = res_ct_real$trend
  )

res_ct_fin <- get_cycle_trend(fin_aligned[[value_col_fin]], detrend_method)
data_merged_fin <- fin_aligned %>%
  dplyr::mutate(
    s          = LOG_s * res_ct_fin$cycle,
    trend_fin  = res_ct_fin$trend
  )

cat("Method =", detrend_method,
    "→ corr(h, trend_real) =", round(stats::cor(data_merged_real$h, data_merged_real$trend_real), 4),
    "  |  corr(s, trend_fin) =", round(stats::cor(data_merged_fin$s, data_merged_fin$trend_fin), 4), "\n")


p1r <- ggplot2::ggplot(data_merged_real, ggplot2::aes(date)) +
  ggplot2::geom_line(ggplot2::aes(y = .data[[value_col_real]]), color="black", size=0.8) +
  ggplot2::geom_line(ggplot2::aes(y = trend_real), color="red", size=0.9) +
  ggplot2::labs(title="(a) Real — Level & Trend", x="Date", y=value_col_real) +
  ggplot2::theme_minimal()
p2r <- ggplot2::ggplot(data_merged_real, ggplot2::aes(date, h)) +
  ggplot2::geom_line(color="black", size=0.8) +
  ggplot2::labs(title="(b) Real — Cyclical Component (h)", x="Date", y="Pct. Deviation") +
  ggplot2::theme_minimal()
p1f <- ggplot2::ggplot(data_merged_fin, ggplot2::aes(date)) +
  ggplot2::geom_line(ggplot2::aes(y = .data[[value_col_fin]]), color="black", size=0.8) +
  ggplot2::geom_line(ggplot2::aes(y = trend_fin), color="red", size=0.9) +
  ggplot2::labs(title="(c) Financial — Level & Trend", x="Date", y=value_col_fin) +
  ggplot2::theme_minimal()
p2f <- ggplot2::ggplot(data_merged_fin, ggplot2::aes(date, s)) +
  ggplot2::geom_line(color="black", size=0.8) +
  ggplot2::labs(title="(d) Financial — Cyclical Component (s)", x="Date", y="Pct. Deviation") +
  ggplot2::theme_minimal()
gridExtra::grid.arrange(p1r, p2r, p1f, p2f, ncol=2)



# 6) Unit-root tests on cycles
cat("=== Real cycle h: unit-root tests ===\n")
print(run_unit_root_tests(data_merged_real$h, model = "level"))
cat("\n=== Financial cycle s: unit-root tests ===\n")
print(run_unit_root_tests(data_merged_fin$s, model = "level"))



# 7) Accumulators + lags (both)
delta <- 1 - (1 - 0.05)^(4 / obs_per_year_real)
N     <- 10 * obs_per_year_real

calc_accumulation_truncated <- function(x, delta, N) {
  w <- delta * (1 - delta)^(0:(N - 1))
  as.numeric(stats::filter(x, w, method = "convolution", sides = 1))
}

data_merged_real <- data_merged_real %>%
  dplyr::arrange(date) %>%
  dplyr::mutate(
    X      = calc_accumulation_truncated(h, delta, N),
    X_lag1 = dplyr::lag(X, 1),
    x_lag1 = dplyr::lag(h, 1),
    x_lag2 = dplyr::lag(h, 2)
  )

data_merged_fin <- data_merged_fin %>%
  dplyr::arrange(date) %>%
  dplyr::mutate(
    Y      = calc_accumulation_truncated(s, delta, N),
    Y_lag1 = dplyr::lag(Y, 1),
    y_lag1 = dplyr::lag(s, 1),
    y_lag2 = dplyr::lag(s, 2)
  )



# 8) Build joint bivariate sample
df_full_bivar <- data_merged_real %>%
  dplyr::select(date, h, X, X_lag1, x_lag1, x_lag2) %>%
  dplyr::inner_join(
    data_merged_fin %>% dplyr::select(date, s, Y, Y_lag1, y_lag1, y_lag2),
    by = "date"
  ) %>%
  dplyr::arrange(date)

prepare_window_bivar <- function(df_full,
                                 start_date = NULL,
                                 end_date   = NULL,
                                 N_required = 10 * obs_per_year_real) {
  df_full <- df_full %>% dplyr::arrange(date)
  warmup  <- max(N_required, 2)
  if (is.null(start_date)) start_date <- df_full$date[warmup + 1L]
  if (is.null(end_date))   end_date   <- max(df_full$date, na.rm = TRUE)
  df_win <- subset(df_full, date >= start_date & date <= end_date)
  df_win <- tidyr::drop_na(df_win, X_lag1, x_lag1, x_lag2, Y_lag1, y_lag1, y_lag2)
  list(
    df_samp_bi  = df_win,
    start_date_f= min(df_win$date),
    end_date_f  = max(df_win$date),
    horizon_data= nrow(df_win)
  )
}

win_bi <- prepare_window_bivar(df_full_bivar)
df_samp_bi   <- win_bi$df_samp_bi
start_date_f <- win_bi$start_date_f
horizon_data <- win_bi$horizon_data

forecast_starter <- 3 * obs_per_year_real
if (1 + forecast_starter > nrow(df_samp_bi))
  stop("Not enough data after harmonization for the requested forecast_starter.")
forecast_origin <- df_samp_bi$date[1 + forecast_starter]








# 9. (Bi)variate nonlinear features
make_row <- function(x_lag1, x_lag2, X_lag1, include_intercept = TRUE) {
  stopifnot(length(x_lag1)==1L, length(x_lag2)==1L, length(X_lag1)==1L)
  stopifnot(all(is.finite(c(x_lag1, x_lag2, X_lag1))))
  xl1 <- as.numeric(x_lag1); xl2 <- as.numeric(x_lag2); XL <- as.numeric(X_lag1)
  out <- data.frame(
    x_lag1 = xl1, x_lag2 = xl2, X_lag1 = XL,
    x_lag1_sq = xl1^2, x_lag2_sq = xl2^2, X_lag1_sq = XL^2,
    x_lag1_cub = xl1^3, x_lag2_cub = xl2^3, X_lag1_cub = XL^3,
    x_lag1_x_lag2 = xl1 * xl2,
    x_lag1_X_lag1 = xl1 * XL,
    x_lag2_X_lag1 = xl2 * XL,
    x_lag1_sq_x_lag2 = (xl1^2) * xl2,
    x_lag1_x_lag2_sq = xl1 * (xl2^2),
    x_lag1_sq_X_lag1 = (xl1^2) * XL,
    x_lag2_sq_X_lag1 = (xl2^2) * XL,
    x_lag1_X_lag1_sq = xl1 * (XL^2),
    x_lag2_X_lag1_sq = xl2 * (XL^2),
    x_lag1_x_lag2_X_lag1 = xl1 * xl2 * XL,
    check.names = FALSE
  )
  if (isTRUE(include_intercept)) out <- cbind("(Intercept)" = 1, out, check.names = FALSE)
  out
}

add_features_bivar <- function(df) {
  df %>%
    dplyr::mutate(
      # x-block (real) 
      x_lag1_sq = x_lag1^2, x_lag2_sq = x_lag2^2, X_lag1_sq = X_lag1^2,
      x_lag1_cub = x_lag1^3, x_lag2_cub = x_lag2^3, X_lag1_cub = X_lag1^3,
      x_lag1_x_lag2 = x_lag1 * x_lag2,
      x_lag1_X_lag1 = x_lag1 * X_lag1,
      x_lag2_X_lag1 = x_lag2 * X_lag1,
      x_lag1_sq_x_lag2 = x_lag1_sq * x_lag2,
      x_lag1_x_lag2_sq = x_lag1 * x_lag2_sq,
      x_lag1_sq_X_lag1 = x_lag1_sq * X_lag1,
      x_lag2_sq_X_lag1 = x_lag2_sq * X_lag1,
      x_lag1_X_lag1_sq = x_lag1 * X_lag1_sq,
      x_lag2_X_lag1_sq = x_lag2 * X_lag1_sq,
      x_lag1_x_lag2_X_lag1 = x_lag1 * x_lag2 * X_lag1,
      
      # y-block (financial) 
      y_lag1_sq = y_lag1^2, y_lag2_sq = y_lag2^2, Y_lag1_sq = Y_lag1^2,
      y_lag1_cub = y_lag1^3, y_lag2_cub = y_lag2^3, Y_lag1_cub = Y_lag1^3,
      y_lag1_y_lag2 = y_lag1 * y_lag2,
      y_lag1_Y_lag1 = y_lag1 * Y_lag1,
      y_lag2_Y_lag1 = y_lag2 * Y_lag1,
      y_lag1_sq_y_lag2 = y_lag1_sq * y_lag2,
      y_lag1_y_lag2_sq = y_lag1 * y_lag2_sq,
      y_lag1_sq_Y_lag1 = y_lag1_sq * Y_lag1,
      y_lag2_sq_Y_lag1 = y_lag2_sq * Y_lag1,
      y_lag1_Y_lag1_sq = y_lag1 * Y_lag1_sq,
      y_lag2_Y_lag1_sq = y_lag2 * Y_lag1_sq,
      y_lag1_y_lag2_Y_lag1 = y_lag1 * y_lag2 * Y_lag1,
      
      # cross-block (new)
      # bilinear cross-effects
      x_lag1_y_lag1 = x_lag1 * y_lag1,
      x_lag1_Y_lag1 = x_lag1 * Y_lag1,
      X_lag1_y_lag1 = X_lag1 * y_lag1,
      
      # selected 3rd-order cross interactions (sparse & economically interpretable)
      x_lag1_sq_y_lag1 = (x_lag1^2) * y_lag1,
      x_lag1_y_lag1_sq = x_lag1 * (y_lag1^2),
      
      x_lag1_sq_Y_lag1 = (x_lag1^2) * Y_lag1,
      x_lag1_Y_lag1_sq = x_lag1 * (Y_lag1^2),
      
      X_lag1_sq_y_lag1 = (X_lag1^2) * y_lag1,
      X_lag1_y_lag1_sq = X_lag1 * (y_lag1^2),
      
      # symmetric terms used by the s-equation
      y_lag1_x_lag1 = y_lag1 * x_lag1,  # alias; same numeric column, explicit label ok
      y_lag1_X_lag1 = y_lag1 * X_lag1,
      Y_lag1_x_lag1 = Y_lag1 * x_lag1,
      
      y_lag1_sq_x_lag1 = (y_lag1^2) * x_lag1,
      y_lag1_x_lag1_sq = y_lag1 * (x_lag1^2),
      
      y_lag1_sq_X_lag1 = (y_lag1^2) * X_lag1,
      y_lag1_X_lag1_sq = y_lag1 * (X_lag1^2),
      
      Y_lag1_sq_x_lag1 = (Y_lag1^2) * x_lag1,
      Y_lag1_x_lag1_sq = Y_lag1 * (x_lag1^2)
    )
}

df_samp_bi <- add_features_bivar(df_samp_bi)

# Guards: explicit lists 
need_x <- c("x_lag1","x_lag2","X_lag1",
            "x_lag1_sq","x_lag2_sq","X_lag1_sq",
            "x_lag1_cub","x_lag2_cub","X_lag1_cub",
            "x_lag1_x_lag2","x_lag1_X_lag1","x_lag2_X_lag1",
            "x_lag1_sq_x_lag2","x_lag1_x_lag2_sq",
            "x_lag1_sq_X_lag1","x_lag2_sq_X_lag1",
            "x_lag1_X_lag1_sq","x_lag2_X_lag1_sq",
            "x_lag1_x_lag2_X_lag1")
need_y <- c("y_lag1","y_lag2","Y_lag1",
            "y_lag1_sq","y_lag2_sq","Y_lag1_sq",
            "y_lag1_cub","y_lag2_cub","Y_lag1_cub",
            "y_lag1_y_lag2","y_lag1_Y_lag1","y_lag2_Y_lag1",
            "y_lag1_sq_y_lag2","y_lag1_y_lag2_sq",
            "y_lag1_sq_Y_lag1","y_lag2_sq_Y_lag1",
            "y_lag1_Y_lag1_sq","y_lag2_Y_lag1_sq",
            "y_lag1_y_lag2_Y_lag1")
need_cross <- c(
  "x_lag1_y_lag1","x_lag1_Y_lag1","X_lag1_y_lag1",
  "x_lag1_sq_y_lag1","x_lag1_y_lag1_sq",
  "x_lag1_sq_Y_lag1","x_lag1_Y_lag1_sq",
  "X_lag1_sq_y_lag1","X_lag1_y_lag1_sq",
  "y_lag1_x_lag1","y_lag1_X_lag1","Y_lag1_x_lag1",
  "y_lag1_sq_x_lag1","y_lag1_x_lag1_sq",
  "y_lag1_sq_X_lag1","y_lag1_X_lag1_sq",
  "Y_lag1_sq_x_lag1","Y_lag1_x_lag1_sq"
)
stopifnot(all(need_x %in% names(df_samp_bi)), all(need_y %in% names(df_samp_bi)))


# 10. Inspect and plot processed series
str(df_samp_bi)
summary(df_samp_bi[, c("h","X","x_lag1","x_lag2","X_lag1","s","Y","y_lag1","y_lag2","Y_lag1")])

pHx <- ggplot2::ggplot(df_samp_bi, ggplot2::aes(x = date)) +
  ggplot2::geom_line(ggplot2::aes(y = h), color = "red",  size = 1) +
  ggplot2::geom_line(ggplot2::aes(y = X), color = "blue", size = 1) +
  ggplot2::labs(title = "Real block: h (red) vs X (blue)", x="Date", y="Value") +
  ggplot2::theme_minimal()
pSy <- ggplot2::ggplot(df_samp_bi, ggplot2::aes(x = date)) +
  ggplot2::geom_line(ggplot2::aes(y = s), color = "red",  size = 1) +
  ggplot2::geom_line(ggplot2::aes(y = Y), color = "blue", size = 1) +
  ggplot2::labs(title = "Financial block: s (red) vs Y (blue)", x="Date", y="Value") +
  ggplot2::theme_minimal()
gridExtra::grid.arrange(pHx, pSy, ncol = 1)



# STAR / VSTAR helpers
G_logistic <- function(z, eta, c) {           # eta = log(gamma)
  gamma <- exp(eta)
  # plogis is stable even for large arguments
  plogis(gamma * (z - c))
}
G_estar <- function(z, eta, c) {              # ESTAR
  gamma <- exp(eta)
  # safe even for big gamma*(z-c)^2 ⇒ tends to 1 smoothly
  1 - exp(-gamma * (z - c)^2)
}

make_z_bi <- function(df, z_var = "x_lag1", v = c(1,0,0, 0,0,0)) {
  if (identical(z_var, "x_lag1")) return(df$x_lag1)
  if (identical(z_var, "y_lag1")) return(df$y_lag1)
  if (identical(z_var, "combo"))  {
    stopifnot(length(v) == 6L)
    return(v[1]*df$x_lag1 + v[2]*df$x_lag2 + v[3]*df$X_lag1 +
             v[4]*df$y_lag1 + v[5]*df$y_lag2 + v[6]*df$Y_lag1)
  }
  stop("Unknown z_var for bivariate STAR: use 'x_lag1', 'y_lag1', or 'combo'.")
}

.vstar_default_regs <- c("x_lag1","x_lag2","X_lag1","y_lag1","y_lag2","Y_lag1")

vstar_eval <- function(theta, df, regs = .vstar_default_regs,
                       spec = c("standard_lstar","minimal_lstar","minimal_estar"),
                       z) {
  spec <- match.arg(spec)
  R <- length(regs)
  if (spec == "standard_lstar") {
    idx <- 0L
    a0_h <- theta[idx <- idx+1]; a1_h <- theta[idx <- idx+1]
    b0_h <- theta[(idx+1):(idx+R)]; idx <- idx+R
    b1_h <- theta[(idx+1):(idx+R)]; idx <- idx+R
    a0_s <- theta[idx <- idx+1]; a1_s <- theta[idx <- idx+1]
    b0_s <- theta[(idx+1):(idx+R)]; idx <- idx+R
    b1_s <- theta[(idx+1):(idx+R)]; idx <- idx+R
    eta  <- theta[idx <- idx+1];    c <- theta[idx <- idx+1]
    G <- G_logistic(z, eta, c)
    Xmat <- as.matrix(df[, regs, drop=FALSE])
    mu_h <- a0_h + a1_h*G
    mu_s <- a0_s + a1_s*G
    beta_h <- matrix(rep(b0_h, each = nrow(df)), nrow=nrow(df)) +
      matrix(rep(b1_h, each = nrow(df)), nrow=nrow(df)) * G
    beta_s <- matrix(rep(b0_s, each = nrow(df)), nrow=nrow(df)) +
      matrix(rep(b1_s, each = nrow(df)), nrow=nrow(df)) * G
    yhat_h <- as.numeric(mu_h + rowSums(Xmat * beta_h))
    yhat_s <- as.numeric(mu_s + rowSums(Xmat * beta_s))
  } else {
    idx <- 0L
    a0_h <- theta[idx <- idx+1]
    b0_h <- theta[(idx+1):(idx+R)]; idx <- idx+R
    d_h  <- theta[(idx+1):(idx+R)]; idx <- idx+R
    a0_s <- theta[idx <- idx+1]
    b0_s <- theta[(idx+1):(idx+R)]; idx <- idx+R
    d_s  <- theta[(idx+1):(idx+R)]; idx <- idx+R
    eta  <- theta[idx <- idx+1];    c <- theta[idx <- idx+1]
    G <- if (spec == "minimal_lstar") G_logistic(z, eta, c) else G_estar(z, eta, c)
    Xmat <- as.matrix(df[, regs, drop=FALSE])
    beta_h <- matrix(rep(b0_h, each = nrow(df)), nrow=nrow(df)) +
      matrix(rep(d_h,  each = nrow(df)), nrow=nrow(df)) * G
    beta_s <- matrix(rep(b0_s, each = nrow(df)), nrow=nrow(df)) +
      matrix(rep(d_s,  each = nrow(df)), nrow=nrow(df)) * G
    yhat_h <- as.numeric(a0_h + rowSums(Xmat * beta_h))
    yhat_s <- as.numeric(a0_s + rowSums(Xmat * beta_s))
  }
  list(yhat_h = yhat_h, yhat_s = yhat_s)
}

estimate_vstar <- function(df,
                           regs = .vstar_default_regs,
                           spec = c("standard_lstar","minimal_lstar","minimal_estar"),
                           z_var = "x_lag1", v = c(1,0,0, 0,0,0),
                           starts_gamma = c(0.25, 0.5, 1, 2, 5, 10),
                           starts_c = c("median","q10","q25","q50","q75","q90"),
                           trace = TRUE, maxit = 2000, compute_se = FALSE) {
  spec <- match.arg(spec)
  
  stopifnot(is.data.frame(df))
  stopifnot(all(c("h","s") %in% names(df)))
  
  # keep only regressors that exist
  regs <- intersect(regs, names(df))
  if (length(regs) < 1L) stop("No valid regressors found in df for regs=...")
  Z <- df[, regs, drop = FALSE]
  
  # build transition (then standardize for numerical stability)
  z_raw <- make_z_bi(df, z_var = z_var, v = v)
  z_mu  <- stats::median(z_raw, na.rm = TRUE)
  z_sd  <- stats::sd(z_raw, na.rm = TRUE); if (!is.finite(z_sd) || z_sd==0) z_sd <- stats::mad(z_raw, na.rm = TRUE)
  z     <- as.numeric((z_raw - z_mu) / z_sd)
  
  # sample mask
  ok <- is.finite(df$h) & is.finite(df$s) & is.finite(z) & apply(Z, 1, function(r) all(is.finite(r)))
  y_h <- df$h[ok]; y_s <- df$s[ok]; Z <- Z[ok, , drop = FALSE]; z <- z[ok]
  n   <- length(y_h); R <- ncol(Z)
  stopifnot(n >= 50, R >= 1)
  
  # Linear warm starts (SUR if available; else OLS per eq) 
  a0h <- 0; a0s <- 0
  b0h <- rep(0, R); b0s <- rep(0, R)
  
  if (requireNamespace("systemfit", quietly = TRUE)) {
    # equation labels MUST NOT have underscores/spaces
    forms <- list(
      hEq = stats::as.formula(paste0("h ~ ", paste(regs, collapse = " + "))),
      sEq = stats::as.formula(paste0("s ~ ", paste(regs, collapse = " + ")))
    )
    df_ok <- cbind.data.frame(h = y_h, s = y_s, Z)
    fit0  <- try(systemfit::systemfit(forms, data = df_ok, method = "SUR"), silent = TRUE)
    if (!inherits(fit0, "try-error")) {
      bh <- stats::coef(fit0$eq$hEq)
      bs <- stats::coef(fit0$eq$sEq)
      a0h <- unname(bh[1]); a0s <- unname(bs[1])
      b0h <- unname(bh[-1]); b0s <- unname(bs[-1])
    }
  }
  # If SUR failed or not available, fall back to OLS
  if (all(!is.finite(b0h)) || all(b0h==0)) {
    lm_h <- stats::lm(y_h ~ ., data = as.data.frame(Z))
    lm_s <- stats::lm(y_s ~ ., data = as.data.frame(Z))
    a0h  <- unname(stats::coef(lm_h)[1]); a0s <- unname(stats::coef(lm_s)[1])
    b0h  <- unname(stats::coef(lm_h)[-1]); b0s <- unname(stats::coef(lm_s)[-1])
  }
  a0h[!is.finite(a0h)] <- 0; a0s[!is.finite(a0s)] <- 0
  b0h[!is.finite(b0h)] <- 0; b0s[!is.finite(b0s)] <- 0
  
  # starting grid for (eta,c), on standardized z
  qfun <- function(z, key)
    switch(key,
           median = stats::median(z, na.rm = TRUE),
           q10    = stats::quantile(z, 0.10, na.rm = TRUE),
           q25    = stats::quantile(z, 0.25, na.rm = TRUE),
           q50    = stats::quantile(z, 0.50, na.rm = TRUE),
           q75    = stats::quantile(z, 0.75, na.rm = TRUE),
           q90    = stats::quantile(z, 0.90, na.rm = TRUE),
           stats::median(z, na.rm = TRUE))
  c_candidates <- unique(vapply(starts_c, function(k) qfun(z, k), numeric(1)))
  g_candidates <- starts_gamma
  
  # parameter vector factories 
  par_linear_standard <- function() {
    c(a0_h = a0h, a1_h = 0,
      b0_h = b0h, b1_h = rep(0, R),
      a0_s = a0s, a1_s = 0,
      b0_s = b0s, b1_s = rep(0, R),
      eta  = 0,    # gamma=1
      c    = 0)    # centered z ⇒ c≈0
  }
  par_linear_minimal <- function() {
    c(a0_h = a0h, b0_h = b0h, d_h = rep(0, R),
      a0_s = a0s, b0_s = b0s, d_s = rep(0, R),
      eta  = 0, c = 0)
  }
  
  # objective (joint SSE)
  obj <- function(th) {
    df_ok <- cbind.data.frame(Z)            # vstar_eval expects the X-matrix columns
    pr <- vstar_eval(th, df_ok, regs = colnames(Z), spec = spec, z = z)
    sum((y_h - pr$yhat_h)^2) + sum((y_s - pr$yhat_s)^2)
  }
  
  # initialize with a safe linear baseline (ALWAYS finite)
  st0 <- if (spec == "standard_lstar") par_linear_standard() else par_linear_minimal()
  best <- list(val = obj(st0), par = st0, conv = 0L)
  if (isTRUE(trace)) message(sprintf("VSTAR(%s) baseline SSE = %.6f", spec, best$val))
  
  # grid search over (gamma,c) + two optimizers per start 
  make_start <- function(g0, c0) {
    if (spec == "standard_lstar") {
      c(a0_h = a0h, a1_h = 0,
        b0_h = b0h, b1_h = rep(0, R),
        a0_s = a0s, a1_s = 0,
        b0_s = b0s, b1_s = rep(0, R),
        eta  = log(g0), c = c0)
    } else {
      c(a0_h = a0h, b0_h = b0h, d_h = rep(0, R),
        a0_s = a0s, b0_s = b0s, d_s = rep(0, R),
        eta  = log(g0), c = c0)
    }
  }
  
  for (g0 in g_candidates) for (c0 in c_candidates) {
    st <- make_start(g0, c0)
    
    # Try BFGS
    opt1 <- try(stats::optim(st, obj, method = "BFGS",
                             control = list(maxit = maxit, reltol = 1e-10)),
                silent = TRUE)
    cand <- NULL
    if (!inherits(opt1, "try-error") && is.finite(opt1$value)) cand <- opt1
    # If BFGS bad, try Nelder–Mead
    if (is.null(cand)) {
      opt2 <- try(stats::optim(st, obj, method = "Nelder-Mead",
                               control = list(maxit = maxit, reltol = 1e-8)),
                  silent = TRUE)
      if (!inherits(opt2, "try-error") && is.finite(opt2$value)) cand <- opt2
    }
    
    # Keep the best feasible candidate
    if (!is.null(cand) && cand$value < best$val) {
      best <- list(val = cand$value, par = cand$par, conv = cand$convergence)
      if (isTRUE(trace)) message(sprintf("  improved: SSE = %.6f  (g0=%.3g, c0=%.3f)", best$val, g0, c0))
    }
  }
  
  theta <- best$par
  
  # fitteds, residuals, criteria (on original scale) 
  df_ok <- cbind.data.frame(Z)
  pr    <- vstar_eval(theta, df_ok, regs = colnames(Z), spec = spec, z = z)
  res_h <- y_h - pr$yhat_h
  res_s <- y_s - pr$yhat_s
  s2    <- (mean(res_h^2) + mean(res_s^2)) / 2
  k     <- length(theta)
  ll    <- -0.5 * (2*n) * (log(2*pi*s2) + 1)
  AIC   <- (2*n) * log(s2) + 2 * k
  BIC   <- (2*n) * log(s2) + k * log(2*n)
  
  # optional SEs
  se <- rep(NA_real_, k)
  if (compute_se && requireNamespace("numDeriv", quietly = TRUE)) {
    f <- function(th) {
      pr <- vstar_eval(th, df_ok, regs = colnames(Z), spec = spec, z = z)
      c(pr$yhat_h, pr$yhat_s)
    }
    G <- try(numDeriv::jacobian(f, theta), silent = TRUE)  # (2n x k)
    if (!inherits(G, "try-error")) {
      V <- try(chol2inv(chol(crossprod(G))) * s2, silent = TRUE)
      if (!inherits(V, "try-error")) se <- sqrt(pmax(diag(V), 0))
    }
  }
  
  nm <- function(side) if (spec == "standard_lstar") {
    c(paste0("a0_", side), paste0("a1_", side),
      as.vector(rbind(paste0("b0_", regs, "_", side),
                      paste0("b1_", regs, "_", side))))
  } else {
    c(paste0("a0_", side),
      paste0("b0_", regs, "_", side),
      paste0("d_",  regs, "_", side))
  }
  names(theta) <- c(nm("h"), nm("s"), "eta", "c")
  names(se)    <- names(theta)
  
  structure(list(
    spec      = spec,
    regs      = regs,
    z_var     = z_var,
    v         = v,
    gfun      = if (spec == "minimal_estar") "ESTAR" else "LSTAR",
    coef      = theta,
    se        = se,
    fitted_h  = pr$yhat_h,
    resid_h   = res_h,
    fitted_s  = pr$yhat_s,
    resid_s   = res_s,
    logLik    = ll, AIC = AIC, BIC = BIC,
    nobs      = n,
    z_center  = z_mu,      # for reference
    z_scale   = z_sd       # for reference
  ), class = "vstar_fit")
}




#### Linear & Polynomial Bivariate Models

# Always create OLS versions for diagnostics
form_lin_h <- h ~ x_lag1 + x_lag2 + X_lag1 + y_lag1 + y_lag2 + Y_lag1
form_lin_s <- s ~ y_lag1 + y_lag2 + Y_lag1 + x_lag1 + x_lag2 + X_lag1
m_lin_h <- stats::lm(form_lin_h, data = df_samp_bi)
m_lin_s <- stats::lm(form_lin_s, data = df_samp_bi)

# If SUR available, fit it additionally (no spaces/underscores in labels!)
if (requireNamespace("systemfit", quietly = TRUE)) {
  eqs <- list(hEq = form_lin_h, sEq = form_lin_s)  # <- labels OK
  sys_lin <- systemfit::systemfit(eqs, data = df_samp_bi, method = "SUR")
  print(summary(sys_lin))
} else {
  print(summary(m_lin_h)); print(summary(m_lin_s))
}


# AR(2) per block (univariate-by-equation)
m_ar2_h <- lm(h ~ x_lag1 + x_lag2, data = df_samp_bi)
m_ar2_s <- lm(s ~ y_lag1 + y_lag2, data = df_samp_bi)

# Minimal (one cubic, + linear cross-effects) 
# keep one cubic on own-cycle; admit *linear* cross-block terms in both eqs
m_min_h <- lm(h ~ x_lag1 + x_lag2 + X_lag1 +
                y_lag1 + y_lag2 + Y_lag1 +
                x_lag1_cub,
              data = df_samp_bi)

m_min_s <- lm(s ~ y_lag1 + y_lag2 + Y_lag1 +
                x_lag1 + x_lag2 + X_lag1 +
                y_lag1_cub,
              data = df_samp_bi)

# Intermediate (own third-order + selected cross interactions)
# own-block: pure cubics + key (x,X) third-order interactions
# cross-block: bilinear x1*y1 plus four stock–cycle 3rd-order cross terms per equation
m_int_h <- lm(h ~ x_lag1 + x_lag2 + X_lag1 +
                y_lag1 + y_lag2 + Y_lag1 +
                x_lag1_cub + x_lag2_cub + X_lag1_cub +
                x_lag1_sq_X_lag1 + x_lag1_X_lag1_sq +
                x_lag1_y_lag1 +
                x_lag1_sq_Y_lag1 + x_lag1_Y_lag1_sq +
                y_lag1_sq_X_lag1 + y_lag1_X_lag1_sq,
              data = df_samp_bi)

m_int_s <- lm(s ~ y_lag1 + y_lag2 + Y_lag1 +
                x_lag1 + x_lag2 + X_lag1 +
                y_lag1_cub + y_lag2_cub + Y_lag1_cub +
                y_lag1_sq_Y_lag1 + y_lag1_Y_lag1_sq +
                y_lag1_x_lag1 +
                y_lag1_sq_X_lag1 + y_lag1_X_lag1_sq +
                x_lag1_sq_Y_lag1 + x_lag1_Y_lag1_sq,
              data = df_samp_bi)

# Full (third-order with cross-block interactions)
# all within-block third-order terms + a compact but rich set of cross-block terms
m_full_h <- lm(h ~
                 # linear (own + cross)
                 x_lag1 + x_lag2 + X_lag1 + y_lag1 + y_lag2 + Y_lag1 +
                 # within-block quadratics & cubics
                 x_lag1_sq + x_lag2_sq + X_lag1_sq +
                 x_lag1_cub + x_lag2_cub + X_lag1_cub +
                 x_lag1_x_lag2 + x_lag1_X_lag1 + x_lag2_X_lag1 +
                 x_lag1_sq_x_lag2 + x_lag1_x_lag2_sq +
                 x_lag1_sq_X_lag1 + x_lag2_sq_X_lag1 +
                 x_lag1_X_lag1_sq + x_lag2_X_lag1_sq +
                 x_lag1_x_lag2_X_lag1 +
                 # cross-block bilinear
                 x_lag1_y_lag1 + x_lag1_Y_lag1 + X_lag1_y_lag1 +
                 # cross-block third-order (cycle–stock & cycle–cycle)
                 x_lag1_sq_y_lag1 + x_lag1_y_lag1_sq +
                 x_lag1_sq_Y_lag1 + x_lag1_Y_lag1_sq +
                 X_lag1_sq_y_lag1 + X_lag1_y_lag1_sq,
               data = df_samp_bi)

m_full_s <- lm(s ~
                 # linear (own + cross)
                 y_lag1 + y_lag2 + Y_lag1 + x_lag1 + x_lag2 + X_lag1 +
                 # within-block quadratics & cubics
                 y_lag1_sq + y_lag2_sq + Y_lag1_sq +
                 y_lag1_cub + y_lag2_cub + Y_lag1_cub +
                 y_lag1_y_lag2 + y_lag1_Y_lag1 + y_lag2_Y_lag1 +
                 y_lag1_sq_y_lag2 + y_lag1_y_lag2_sq +
                 y_lag1_sq_Y_lag1 + y_lag2_sq_Y_lag1 +
                 y_lag1_Y_lag1_sq + y_lag2_Y_lag1_sq +
                 y_lag1_y_lag2_Y_lag1 +
                 # cross-block bilinear
                 y_lag1_x_lag1 + y_lag1_X_lag1 + Y_lag1_x_lag1 +
                 # cross-block third-order (cycle–stock & cycle–cycle)
                 y_lag1_sq_x_lag1 + y_lag1_x_lag1_sq +
                 y_lag1_sq_X_lag1 + y_lag1_X_lag1_sq +
                 Y_lag1_sq_x_lag1 + Y_lag1_x_lag1_sq,
               data = df_samp_bi)

print(summary(m_min_h)); print(summary(m_min_s))
print(summary(m_int_h)); print(summary(m_int_s))
print(summary(m_full_h)); print(summary(m_full_s))



#### VSTAR fits (shared transition)

vstar_std  <- estimate_vstar(df_samp_bi, spec = "standard_lstar",
                             regs = .vstar_default_regs, z_var = "x_lag1",
                             trace = TRUE)
vstar_min  <- estimate_vstar(df_samp_bi, spec = "minimal_lstar",
                             regs = .vstar_default_regs, z_var = "x_lag1",
                             trace = TRUE)
vestar_min <- estimate_vstar(df_samp_bi, spec = "minimal_estar",
                             regs = .vstar_default_regs, z_var = "x_lag1",
                             trace = TRUE)

print(vstar_std); print(vstar_min); print(vestar_min)







### AUTOCORRELATION TESTS (linear-in-parameters — run on OLS per eq)
if (!requireNamespace("lmtest", quietly = TRUE)) stop("Please install lmtest")
library(lmtest)

dw_h <- lmtest::dwtest(m_lin_h); dw_s <- lmtest::dwtest(m_lin_s)
cat("\nDurbin–Watson (h):\n"); print(dw_h)
cat("\nDurbin–Watson (s):\n"); print(dw_s)

bg_h <- lmtest::bgtest(m_lin_h, order = 4); bg_s <- lmtest::bgtest(m_lin_s, order = 4)
cat("\nBreusch–Godfrey (h):\n"); print(bg_h)
cat("\nBreusch–Godfrey (s):\n"); print(bg_s)

lb_h <- stats::Box.test(stats::residuals(m_lin_h), lag = 20, type = "Ljung-Box")
lb_s <- stats::Box.test(stats::residuals(m_lin_s), lag = 20, type = "Ljung-Box")
cat("\nLjung–Box (h):\n"); print(lb_h)
cat("\nLjung–Box (s):\n"); print(lb_s)
stats::acf(stats::residuals(m_lin_h), main = "ACF of Residuals (h)")
stats::acf(stats::residuals(m_lin_s), main = "ACF of Residuals (s)")


# NORMALITY 
resid_h <- stats::residuals(m_full_h); graphics::hist(resid_h)
stats::qqnorm(resid_h); stats::qqline(resid_h); print(stats::shapiro.test(resid_h))
resid_s <- stats::residuals(m_full_s); graphics::hist(resid_s)
stats::qqnorm(resid_s); stats::qqline(resid_s); print(stats::shapiro.test(resid_s))


# VSTAR restriction tests (indicative Wald - eta,c not ID under H0)
R_vstar_standard <- function(theta_names, regs = .vstar_default_regs) {
  tar <- c("a1_h","a1_s", paste0("b1_", regs, "_h"), paste0("b1_", regs, "_s"))
  idx <- match(tar, theta_names)
  R <- matrix(0, nrow = length(idx), ncol = length(theta_names))
  R[cbind(seq_along(idx), idx)] <- 1
  R
}
R_vstar_minimal <- function(theta_names, regs = .vstar_default_regs) {
  tar <- c(paste0("d_", regs, "_h"), paste0("d_", regs, "_s"))
  idx <- match(tar, theta_names)
  R <- matrix(0, nrow = length(idx), ncol = length(theta_names))
  R[cbind(seq_along(idx), idx)] <- 1
  R
}
vcov_vstar <- function(vstar_obj, df, regs = .vstar_default_regs) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Please install numDeriv")
  Z <- as.matrix(df[, regs, drop=FALSE])
  z <- make_z_bi(df, z_var = vstar_obj$z_var, v = vstar_obj$v)
  ok<- is.finite(df$h) & is.finite(df$s) & is.finite(z) & apply(Z,1,function(r) all(is.finite(r)))
  Z <- Z[ok,,drop=FALSE]; z <- z[ok]
  f <- function(th) {
    pr <- vstar_eval(th, as.data.frame(Z), regs=colnames(Z),
                     spec=vstar_obj$spec, z=z)
    c(pr$yhat_h, pr$yhat_s)
  }
  G <- numDeriv::jacobian(f, vstar_obj$coef)
  s2 <- (mean(vstar_obj$resid_h^2) + mean(vstar_obj$resid_s^2))/2
  V  <- chol2inv(chol(crossprod(G))) * s2
  dimnames(V) <- list(names(vstar_obj$coef), names(vstar_obj$coef))
  V
}
wald_vstar <- function(vstar_obj, df, which = c("standard","minimal")) {
  which <- match.arg(which)
  V  <- vcov_vstar(vstar_obj, df, regs = vstar_obj$regs)
  th <- vstar_obj$coef
  R  <- if (which=="standard") R_vstar_standard(names(th), regs = vstar_obj$regs)
  else                   R_vstar_minimal(names(th),  regs = vstar_obj$regs)
  Rth <- drop(R %*% th)
  RVRT<- R %*% V %*% t(R)
  stat<- drop(crossprod(Rth, solve(RVRT, Rth)))
  dfq <- nrow(R)
  p   <- stats::pchisq(stat, df=dfq, lower.tail=FALSE)
  list(stat=stat, df=dfq, p.value=p,
       note="Indicative Wald; Davies issue applies (eta,c not ID under H0). Prefer bootstrap LR.")
}
wstd <- wald_vstar(vstar_std,  df_samp_bi, which="standard")
wmin <- wald_vstar(vstar_min,  df_samp_bi, which="minimal")
west <- wald_vstar(vestar_min, df_samp_bi, which="minimal")
cat("\n[Indicative] Wald VSTAR standard (all G parts = 0):  chi2(", wstd$df, ")=",
    round(wstd$stat,3), "  p=", signif(wstd$p.value,3), "\n", sep="")
cat("[Indicative] Wald VSTAR minimal (all d's = 0):       chi2(", wmin$df, ")=",
    round(wmin$stat,3), "  p=", signif(wmin$p.value,3), "\n", sep="")
cat("[Indicative] Wald VESTAR minimal (all d's = 0):      chi2(", west$df, ")=",
    round(west$stat,3), "  p=", signif(west$p.value,3), "\n", sep="")
cat("\nNOTE: Wald is only indicative for VSTAR due to Davies (c,gamma not ID under H0).\n",
    "Use a bootstrap LR adapted to the joint system for reliable p-values.\n", sep="")












### LASSO fits for h and s (glmnet) → coef_lasso_h / coef_lasso_s

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Please install glmnet: install.packages('glmnet')")
}

# Full regressor set the simulator understands (must match features_vec names)
regs_all <- c(need_x, need_y, need_cross)
regs_all <- intersect(regs_all, names(df_samp_bi))
nzv <- vapply(df_samp_bi[regs_all], function(v) sd(v, na.rm = TRUE) > 0, logical(1))
regs_all <- regs_all[nzv]

# Common mask of complete cases for all predictors + both targets (safe & simple)
ok <- stats::complete.cases(df_samp_bi[, c("h","s", regs_all), drop = FALSE])
X  <- as.matrix(df_samp_bi[ok, regs_all])
y_h <- df_samp_bi$h[ok]
y_s <- df_samp_bi$s[ok]

# Cross-validated LASSO; coefficients come back on the original scale
set.seed(123)
cv_h <- glmnet::cv.glmnet(X, y_h, alpha = 1, family = "gaussian",
                          intercept = TRUE, standardize = TRUE, nfolds = 10)
cv_s <- glmnet::cv.glmnet(X, y_s, alpha = 1, family = "gaussian",
                          intercept = TRUE, standardize = TRUE, nfolds = 10)

# Pull the chosen coefficients (lambda.min or lambda.1se)
co_h <- coef(cv_h, s = "lambda.min")
co_s <- coef(cv_s, s = "lambda.min")

# Turn into named numeric vectors 
coef_lasso_h <- as.numeric(co_h); names(coef_lasso_h) <- rownames(co_h)
coef_lasso_s <- as.numeric(co_s); names(coef_lasso_s) <- rownames(co_s)
















################################################################################
# LOCAL STABILITY & STEADY STATES - Bivariate (h,s) with stocks (X,Y)
################################################################################

`%||%` <- function(a,b) if (length(a)==0 || is.na(a)) b else a
rho     <- function(M)  max(Mod(eigen(M, only.values = TRUE)$values))
eigvals <- function(M)  eigen(M, only.values = TRUE)$values

# 0) Pull coefficients from  fitted models 
get_coefs <- function(obj) {
  # systemfit eq object OR lm: coef() works for both
  out <- stats::coef(obj)
  out[!is.na(out)]
}

specs_bi <- list(
  AR2          = list(h = get_coefs(m_ar2_h),  s = get_coefs(m_ar2_s)),
  Linear       = list(h = get_coefs(m_lin_h),  s = get_coefs(m_lin_s)),
  Minimal      = list(h = get_coefs(m_min_h),  s = get_coefs(m_min_s)),
  Intermediate = list(h = get_coefs(m_int_h),  s = get_coefs(m_int_s)),
  Full         = list(h = get_coefs(m_full_h), s = get_coefs(m_full_s))
)

# 1) Regressor factory at a *state* (h1,s1,h2,s2,X1,Y1)
#    Includes every regressor 
# 1) Regressor factory at a *state* (h1,s1,h2,s2,X1,Y1)
#    Now includes the *same cross-block features* as the bivariate polynomial specs.
row_from_state <- function(h1, s1, h2, s2, X1, Y1) {
  # x-block (real)
  xl1 <- h1; xl2 <- h2; XL <- X1
  xl1_sq <- xl1^2; xl2_sq <- xl2^2; XL_sq <- XL^2
  xl1_cu <- xl1^3; xl2_cu <- xl2^3; XL_cu <- XL^3
  
  # y-block (financial)
  yl1 <- s1; yl2 <- s2; YL <- Y1
  yl1_sq <- yl1^2; yl2_sq <- yl2^2; YL_sq <- YL^2
  yl1_cu <- yl1^3; yl2_cu <- yl2^3; YL_cu <- YL^3
  
  c(
    `(Intercept)` = 1,
    
    # ----- linear (own + cross, to match your model formulas) -----
    x_lag1 = xl1, x_lag2 = xl2, X_lag1 = XL,
    y_lag1 = yl1, y_lag2 = yl2, Y_lag1 = YL,
    
    # ----- x-block nonlinear (own-block) -----
    x_lag1_sq = xl1_sq, x_lag2_sq = xl2_sq, X_lag1_sq = XL_sq,
    x_lag1_cub = xl1_cu, x_lag2_cub = xl2_cu, X_lag1_cub = XL_cu,
    x_lag1_x_lag2 = xl1*xl2, x_lag1_X_lag1 = xl1*XL, x_lag2_X_lag1 = xl2*XL,
    x_lag1_sq_x_lag2 = xl1_sq*xl2, x_lag1_x_lag2_sq = xl1*xl2_sq,
    x_lag1_sq_X_lag1 = xl1_sq*XL, x_lag2_sq_X_lag1 = xl2_sq*XL,
    x_lag1_X_lag1_sq = xl1*XL_sq, x_lag2_X_lag1_sq = xl2*XL_sq,
    x_lag1_x_lag2_X_lag1 = xl1*xl2*XL,
    
    # ----- y-block nonlinear (own-block) -----
    y_lag1_sq = yl1_sq, y_lag2_sq = yl2_sq, Y_lag1_sq = YL_sq,
    y_lag1_cub = yl1_cu, y_lag2_cub = yl2_cu, Y_lag1_cub = YL_cu,
    y_lag1_y_lag2 = yl1*yl2, y_lag1_Y_lag1 = yl1*YL, y_lag2_Y_lag1 = yl2*YL,
    y_lag1_sq_y_lag2 = yl1_sq*yl2, y_lag1_y_lag2_sq = yl1*yl2_sq,
    y_lag1_sq_Y_lag1 = yl1_sq*YL, y_lag2_sq_Y_lag1 = yl2_sq*YL,
    y_lag1_Y_lag1_sq = yl1*YL_sq, y_lag2_Y_lag1_sq = yl2*YL_sq,
    y_lag1_y_lag2_Y_lag1 = yl1*yl2*YL,
    
    # ===== cross-block interactions (new; keep names aligned with your fit) =====
    # bilinear
    x_lag1_y_lag1 = xl1*yl1,
    x_lag1_Y_lag1 = xl1*YL,
    X_lag1_y_lag1 = XL*yl1,
    
    # third-order: cycle–cycle and cycle–stock (h on s-side)
    x_lag1_sq_y_lag1 = (xl1^2)*yl1,
    x_lag1_y_lag1_sq = xl1*(yl1^2),
    x_lag1_sq_Y_lag1 = (xl1^2)*YL,
    x_lag1_Y_lag1_sq = xl1*(YL^2),
    X_lag1_sq_y_lag1 = (XL^2)*yl1,
    X_lag1_y_lag1_sq = XL*(yl1^2),
    
    # symmetric aliases for s-equation (same numerics, clearer labels)
    y_lag1_x_lag1 = yl1*xl1,
    y_lag1_X_lag1 = yl1*XL,
    Y_lag1_x_lag1 = YL*xl1,
    
    y_lag1_sq_x_lag1 = (yl1^2)*xl1,
    y_lag1_x_lag1_sq = yl1*(xl1^2),
    y_lag1_sq_X_lag1 = (yl1^2)*XL,
    y_lag1_X_lag1_sq = yl1*(XL^2),
    Y_lag1_sq_x_lag1 = (YL^2)*xl1,
    Y_lag1_x_lag1_sq = YL*(xl1^2)
  )
}

# 2) Next-period maps h_t, s_t given coefficients (works for any subset of names)
predict_h <- function(beta_h, h1,s1,h2,s2,X1,Y1) {
  r <- row_from_state(h1,s1,h2,s2,X1,Y1)
  nm <- intersect(names(beta_h), names(r))
  sum(beta_h[nm] * r[nm])
}
predict_s <- function(beta_s, h1,s1,h2,s2,X1,Y1) {
  r <- row_from_state(h1,s1,h2,s2,X1,Y1)
  nm <- intersect(names(beta_s), names(r))
  sum(beta_s[nm] * r[nm])
}

# 3) State transition H(z) and a simple finite-difference Jacobian
H_builder <- function(beta_h, beta_s, delta) {
  function(z) {
    h1 <- z[1]; s1 <- z[2]; h2 <- z[3]; s2 <- z[4]; X1 <- z[5]; Y1 <- z[6]
    h_next <- predict_h(beta_h, h1,s1,h2,s2,X1,Y1)
    s_next <- predict_s(beta_s, h1,s1,h2,s2,X1,Y1)
    X_next <- (1 - delta) * X1 + delta * h_next
    Y_next <- (1 - delta) * Y1 + delta * s_next
    c(h_next, s_next, h1, s1, X_next, Y_next)
  }
}

jacobian_fd <- function(F, z0, eps = 1e-6) {
  f0 <- F(z0)
  p  <- length(z0); m <- length(f0)
  J  <- matrix(0, m, p)
  for (j in seq_len(p)) {
    ej <- rep(0, p); ej[j] <- eps
    fp <- F(z0 + ej); fm <- F(z0 - ej)
    J[, j] <- (fp - fm) / (2 * eps)
  }
  J
}

# 4) Steady states: solve F(h*,s*) - (h*,s*) = 0  (reduced system with X*=h*, Y*=s*, h2=h*, s2=s*)
#    We minimize the squared norm; multi-start; no extra packages required
F_reduced <- function(beta_h, beta_s, hs) {
  h <- hs[1]; s <- hs[2]
  h_next <- predict_h(beta_h, h, s, h, s, h, s)
  s_next <- predict_s(beta_s, h, s, h, s, h, s)
  c(h_next - h, s_next - s)
}
obj_reduced <- function(beta_h, beta_s, hs) sum(F_reduced(beta_h, beta_s, hs)^2)

find_steady_states <- function(beta_h, beta_s, starts, tol = 1e-8) {
  sols <- list()
  vals <- c()
  for (k in seq_len(nrow(starts))) {
    st <- as.numeric(starts[k, ])
    opt <- try(stats::optim(st, fn = function(u) obj_reduced(beta_h, beta_s, u),
                            method = "BFGS", control = list(reltol = 1e-12, maxit = 5000)),
               silent = TRUE)
    if (inherits(opt, "try-error")) next
    if (opt$convergence != 0) next
    if (opt$value > 1e-6) next
    cand <- opt$par
    # de-duplicate
    if (!length(sols)) {
      sols[[1]] <- cand; vals[1] <- opt$value
    } else {
      d <- vapply(sols, function(v) sqrt(sum((v - cand)^2)), 0.0)
      if (min(d) > 1e-4) { sols[[length(sols)+1]] <- cand; vals[length(vals)+1] <- opt$value }
    }
  }
  sols
}

# Build multi-starts around sample means ± 5σ
mu_h <- mean(df_samp_bi$h, na.rm = TRUE);  sd_h <- sd(df_samp_bi$h, na.rm = TRUE) %||% 1
mu_s <- mean(df_samp_bi$s, na.rm = TRUE);  sd_s <- sd(df_samp_bi$s, na.rm = TRUE) %||% 1
grid_h <- seq(mu_h - 5*sd_h, mu_h + 5*sd_h, length.out = 7L)
grid_s <- seq(mu_s - 5*sd_s, mu_s + 5*sd_s, length.out = 7L)
starts  <- as.matrix(expand.grid(h = grid_h[c(1,4,7)], s = grid_s[c(1,4,7)]))  # 3×3 = 9 starts

# 5) MAIN LOOP: spectra at zero and at each steady state
results_zero <- list()
results_ss   <- list()

for (nm in names(specs_bi)) {
  cat("\n================  MODEL:", nm, "  ================\n")
  bh <- specs_bi[[nm]]$h
  bs <- specs_bi[[nm]]$s
  H  <- H_builder(bh, bs, delta)
  
  # Jacobian at the zero state (all components = 0)
  z0 <- rep(0, 6L)
  J0 <- jacobian_fd(H, z0)
  ev0 <- eigvals(J0)
  cat("Eigenvalues at (0,0,0,0,0,0):", round(ev0, 6), "\n")
  cat("  → max |λ|:", round(max(Mod(ev0)), 6),
      " | max Im:", round(max(abs(Im(ev0))), 6), "\n")
  results_zero[[nm]] <- list(J = J0, ev = ev0, rho = rho(J0), imag_max = max(abs(Im(ev0))))
  
  # Non-zero steady states (reduced 2D system with X*=h*, Y*=s*, and lags=h*,s*)
  ss_list <- find_steady_states(bh, bs, starts)
  if (!length(ss_list)) {
    cat("No non-zero steady states found in the search region.\n")
  } else {
    ss_mat <- do.call(rbind, ss_list)
    colnames(ss_mat) <- c("h_star", "s_star")
    print(round(ss_mat, 6))
    for (i in seq_len(nrow(ss_mat))) {
      hss <- ss_mat[i, 1]; sss <- ss_mat[i, 2]
      zss <- c(hss, sss, hss, sss, hss, sss)  # (h1,s1,h2,s2,X1,Y1) at SS
      Jss <- jacobian_fd(H, zss)
      evs <- eigvals(Jss)
      cat(sprintf("SS #%d at (h*,s*) = (%.6f, %.6f)\n", i, hss, sss))
      cat("  Eigenvalues:", round(evs, 6), "\n")
      cat("  → max |λ|:", round(max(Mod(evs)), 6),
          " | max Im:", round(max(abs(Im(evs))), 6), "\n")
      results_ss[[paste(nm, i, sep = "#")]] <- list(
        h_star = hss, s_star = sss, J = Jss, ev = evs,
        rho = rho(Jss), imag_max = max(abs(Im(evs)))
      )
    }
  }
}

# 6) Panel B (λmax & Imag part at zero) and (optional) merge with Panel A
panelB <- data.frame(
  Model              = names(results_zero),
  Lambda_max_at_0    = vapply(results_zero, function(x) x$rho, 0.0),
  Imag_part_max_at_0 = vapply(results_zero, function(x) x$imag_max, 0.0),
  row.names = NULL
)

cat("\n**Panel B – Local Stability at (h*,s*) = (0,0)**\n")
print(panelB)

if (exists("results_df")) {
  extended3 <- merge(results_df, panelB, by = "Model", all = TRUE)
  cat("\n**Extended Table 3 – Panel A + Panel B (bivariate)**\n")
  print(extended3)
}

# 7) STAR/VSTAR local stability via numerical Jacobian
#     If want entries for VSTAR, it can evaluate H using vstar_eval().
add_vstar_to_panelB <- function(vstar_obj, label, z_var = vstar_obj$z_var, regs = vstar_obj$regs) {
  # Build a one-row reg DF from a state vector z
  make_df_from_z <- function(z) {
    h1 <- z[1]; s1 <- z[2]; h2 <- z[3]; s2 <- z[4]; X1 <- z[5]; Y1 <- z[6]
    r  <- as.list(row_from_state(h1,s1,h2,s2,X1,Y1))
    as.data.frame(as.list(r[regs]), check.names = FALSE)
  }
  H_vstar <- function(z) {
    df1 <- make_df_from_z(z)
    # transition variable value
    zval <- switch(z_var,
                   x_lag1 = z[1],
                   y_lag1 = z[2],
                   combo  = stop("combo z_var needs 'v' weights; extend as needed."),
                   z[1])
    pr  <- vstar_eval(vstar_obj$coef, df1, regs = regs, spec = vstar_obj$spec, z = zval)
    h_next <- pr$yhat_h
    s_next <- pr$yhat_s
    X_next <- (1 - delta) * z[5] + delta * h_next
    Y_next <- (1 - delta) * z[6] + delta * s_next
    c(h_next, s_next, z[1], z[2], X_next, Y_next)
  }
  z0 <- rep(0, 6L)
  J0 <- jacobian_fd(H_vstar, z0)
  ev <- eigvals(J0)
  data.frame(
    Model = label,
    Lambda_max_at_0    = rho(J0),
    Imag_part_max_at_0 = max(abs(Im(ev))),
    row.names = NULL
  )
}

 panelB_star <- NULL
 if (exists("vstar_std"))  panelB_star <- rbind(panelB_star, add_vstar_to_panelB(vstar_std,  "VSTAR standard"))
 if (exists("vstar_min"))  panelB_star <- rbind(panelB_star, add_vstar_to_panelB(vstar_min,  "VSTAR minimal"))
 if (exists("vestar_min")) panelB_star <- rbind(panelB_star, add_vstar_to_panelB(vestar_min, "VESTAR minimal"))
 if (!is.null(panelB_star)) {
   cat("\n**Panel B – Local Stability at 0 (including VSTAR)**\n"); print(panelB_star)
   panelB_all <- rbind(panelB, panelB_star)
   if (exists("results_df")) {
     extended3_star <- merge(results_df, panelB_all, by = "Model", all = TRUE)
     cat("\n**Extended Table 3 – with VSTAR**\n"); print(extended3_star)
   }
 }
 
 # helpers 
 dw_stat <- function(e) {
   e <- as.numeric(e); e <- e[is.finite(e)]
   if (length(e) < 3) return(NA_real_)
   sum(diff(e)^2, na.rm = TRUE) / sum(e^2, na.rm = TRUE)
 }
 
 panelA_from_vstar <- function(vobj, df, label = vobj$spec) {
   # rebuild the mask used inside estimate_vstar so y aligns with resid_*
   Z <- as.matrix(df[, vobj$regs, drop = FALSE])
   z <- make_z_bi(df, z_var = vobj$z_var, v = vobj$v)
   ok <- is.finite(df$h) & is.finite(df$s) & is.finite(z) &
     apply(Z, 1, function(r) all(is.finite(r)))
   y_h <- df$h[ok]; y_s <- df$s[ok]
   eh  <- vobj$resid_h; es <- vobj$resid_s
   n   <- length(y_h); k <- length(vobj$coef)
   
   SSE <- sum(eh^2) + sum(es^2)
   SST <- sum((y_h - mean(y_h))^2) + sum((y_s - mean(y_s))^2)
   R2  <- if (SST > 0) 1 - SSE / SST else NA_real_
   
   # system-style adjusted R2
   df_num <- max(2*n - k, 1)
   df_den <- max(2*n - 1, 1)
   Adj_R2 <- 1 - ( (SSE/df_num) / (SST/df_den) )
   
   DW <- mean(c(dw_stat(eh), dw_stat(es)), na.rm = TRUE)
   
   data.frame(Model = label, R2 = R2, Adj_R2 = Adj_R2, DW = DW, row.names = NULL)
 }
 
 # build Panel A rows for VSTAR
 panelA_vstar <- do.call(rbind, Filter(Negate(is.null), list(
   if (exists("vstar_std"))  panelA_from_vstar(vstar_std,  df_samp_bi, "VSTAR standard"),
   if (exists("vstar_min"))  panelA_from_vstar(vstar_min,  df_samp_bi, "VSTAR minimal"),
   if (exists("vestar_min")) panelA_from_vstar(vestar_min, df_samp_bi, "VESTAR minimal")
 )))
 
 # append to existing Panel A
 if (exists("results_df")) {
   results_df <- dplyr::bind_rows(results_df, panelA_vstar)
 } else {
   results_df <- panelA_vstar
 }
 
 # harmonize labels to avoid duplicate AR rows 
 fix_labels <- function(d) { d$Model <- sub("^AR2$", "AR(2)", d$Model); d }
 results_df   <- fix_labels(results_df)
 panelB       <- fix_labels(panelB)
 if (exists("panelB_star") && !is.null(panelB_star)) panelB_star <- fix_labels(panelB_star)
 
 # final merge & print 
 panelB_all <- if (exists("panelB_star") && !is.null(panelB_star)) rbind(panelB, panelB_star) else panelB
 extended3_star <- merge(results_df, panelB_all, by = "Model", all = TRUE)
 cat("\n**Extended Table 3 – with VSTAR**\n"); print(extended3_star)
 
 
 
 
 

 

 
 
 
 
 
 
 ################################################################################
 # DETERMINISTIC FORECAST - BIVARIATE (h,s) with stocks (H,Y) 
 ################################################################################
 
 # 0) GLOBALS & helpers 
 horizon <- 1000 * obs_per_year_real / 4  # simulate ≈ 1,000 quarters
 
 # Small helpers
 `%||%` <- function(a, b) {
   if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
 }
 
 # Step size for date sequences
 ## when need a single‑step date increment
 increment <- if (freq_real=="q") months(3) else
   if (freq_real=="m") months(1) else weeks(1)
 
 
 cat(obs_per_year_real, obs_unit_real, delta, N, horizon, "\n")
 
 
 # 2) Robust initializer: works even if s/y_lag*/Y_lag1 are absent in df
 build_initial_objects_bi <- function(df, origin_date) {
   stopifnot("date" %in% names(df))
   if (!(origin_date %in% df$date)) {
     # pick the last date <= origin_date if exact match is missing
     cand <- which(df$date <= origin_date)
     if (!length(cand)) stop("origin_date precedes all rows in df")
     origin_date <- df$date[max(cand)]
   }
   i0 <- which(df$date == origin_date)[1]
   
   get1 <- function(nm, default = 0) {
     if (nm %in% names(df)) {
       v <- as.numeric(df[[nm]][i0])
       if (!is.finite(v)) default else v
     } else default
   }
   
   list(
     # x-block (must exist in setup)
     h0  = get1("h",      0),
     h_1 = get1("x_lag1", 0),
     h_2 = get1("x_lag2", 0),
     H_1 = get1("X_lag1", 0),
     # y-block (may be missing in data_merged_real → safely default to 0)
     s0  = get1("s",      0),
     s_1 = get1("y_lag1", 0),
     s_2 = get1("y_lag2", 0),
     Y_1 = get1("Y_lag1", 0)
   )
 }
 

  # 3) Feature factory that ALWAYS returns numeric scalars (no NULLs)
 features_vec <- function(h1,s1,h2,s2,H1,Y1) {
   h1 <- as.numeric(h1 %||% 0); s1 <- as.numeric(s1 %||% 0)
   h2 <- as.numeric(h2 %||% 0); s2 <- as.numeric(s2 %||% 0)
   H1 <- as.numeric(H1 %||% 0); Y1 <- as.numeric(Y1 %||% 0)
   
   # x-block
   xl1 <- h1; xl2 <- h2; XL <- H1
   xl1_sq <- xl1^2; xl2_sq <- xl2^2; XL_sq <- XL^2
   xl1_cu <- xl1^3; xl2_cu <- xl2^3; XL_cu <- XL^3
   # y-block
   yl1 <- s1; yl2 <- s2; YL <- Y1
   yl1_sq <- yl1^2; yl2_sq <- yl2^2; YL_sq <- YL^2
   yl1_cu <- yl1^3; yl2_cu <- yl2^3; YL_cu <- YL^3
   
   c(
     `(Intercept)` = 1,
     # linear
     x_lag1 = xl1, x_lag2 = xl2, X_lag1 = XL,
     y_lag1 = yl1, y_lag2 = yl2, Y_lag1 = YL,
     # x-block poly & interactions
     x_lag1_sq = xl1_sq, x_lag2_sq = xl2_sq, X_lag1_sq = XL_sq,
     x_lag1_cub = xl1_cu, x_lag2_cub = xl2_cu, X_lag1_cub = XL_cu,
     x_lag1_x_lag2 = xl1*xl2, x_lag1_X_lag1 = xl1*XL, x_lag2_X_lag1 = xl2*XL,
     x_lag1_sq_x_lag2 = xl1_sq*xl2, x_lag1_x_lag2_sq = xl1*xl2_sq,
     x_lag1_sq_X_lag1 = xl1_sq*XL, x_lag2_sq_X_lag1 = xl2_sq*XL,
     x_lag1_X_lag1_sq = xl1*XL_sq, x_lag2_X_lag1_sq = xl2*XL_sq,
     x_lag1_x_lag2_X_lag1 = xl1*xl2*XL,
     # y-block poly & interactions
     y_lag1_sq = yl1_sq, y_lag2_sq = yl2_sq, Y_lag1_sq = YL_sq,
     y_lag1_cub = yl1_cu, y_lag2_cub = yl2_cu, Y_lag1_cub = YL_cu,
     y_lag1_y_lag2 = yl1*yl2, y_lag1_Y_lag1 = yl1*YL, y_lag2_Y_lag1 = yl2*YL,
     y_lag1_sq_y_lag2 = yl1_sq*yl2, y_lag1_y_lag2_sq = yl1*yl2_sq,
     y_lag1_sq_Y_lag1 = yl1_sq*YL, y_lag2_sq_Y_lag1 = yl2_sq*YL,
     y_lag1_Y_lag1_sq = yl1*YL_sq, y_lag2_Y_lag1_sq = yl2*YL_sq,
     y_lag1_y_lag2_Y_lag1 = yl1*yl2*YL
   )
 }
 
 predict_from_coefs <- function(beta, h1,s1,h2,s2,H1,Y1) {
   f  <- features_vec(h1,s1,h2,s2,H1,Y1)
   nm <- intersect(names(beta), names(f))
   if (!length(nm)) return(0)  # all coefs missing here → predict 0
   sum(beta[nm] * f[nm])
 }
 
 
 # 4) TRUE bivariate recursion — coerce scalars once up-front
 simulate_core_recursive_bi <- function(h1,h2,H1, s1,s2,Y1, horizon, delta,
                                        law_h, law_s) {
   h1 <- as.numeric(h1 %||% 0); h2 <- as.numeric(h2 %||% 0); H1 <- as.numeric(H1 %||% 0)
   s1 <- as.numeric(s1 %||% 0); s2 <- as.numeric(s2 %||% 0); Y1 <- as.numeric(Y1 %||% 0)
   
   out_h <- out_s <- out_H <- out_Y <- numeric(horizon)
   for (t in 1:horizon) {
     h_t <- as.numeric(law_h(h1,h2,H1, s1,s2,Y1)); if (!is.finite(h_t)) h_t <- 0
     s_t <- as.numeric(law_s(h1,h2,H1, s1,s2,Y1)); if (!is.finite(s_t)) s_t <- 0
     H_t <- (1 - delta) * H1 + delta * h_t
     Y_t <- (1 - delta) * Y1 + delta * s_t
     out_h[t] <- h_t; out_s[t] <- s_t; out_H[t] <- H_t; out_Y[t] <- Y_t
     # roll state
     h2 <- h1; h1 <- h_t; H1 <- H_t
     s2 <- s1; s1 <- s_t; Y1 <- Y_t
   }
   data.frame(Time = 1:horizon, h = out_h, s = out_s, H = out_H, Y = out_Y)
 }
 

 # 5) Wrapper unchanged — just uses predict_from_coefs() above
 simulate_forecast_from_coefs_bi <- function(beta_h, beta_s,
                                             h0,h_1,H_1, s0,s_1,Y_1,
                                             horizon, delta) {
   law_h <- function(h1,h2,H1, s1,s2,Y1) predict_from_coefs(beta_h, h1,s1,h2,s2,H1,Y1)
   law_s <- function(h1,h2,H1, s1,s2,Y1) predict_from_coefs(beta_s, h1,s1,h2,s2,H1,Y1)
   simulate_core_recursive_bi(h0,h_1,H_1, s0,s_1,Y_1, horizon, delta, law_h, law_s)
 }
 
 simulate_forecast_from_models_bi <- function(model_h, model_s,
                                              h0,h_1,H_1, s0,s_1,Y_1,
                                              horizon, delta) {
   bh <- stats::coef(model_h); bh <- bh[!is.na(bh)]
   bs <- stats::coef(model_s); bs <- bs[!is.na(bs)]
   simulate_forecast_from_coefs_bi(bh, bs, h0,h_1,H_1, s0,s_1,Y_1, horizon, delta)
 }
 
 
 # 4.3 VSTAR / VESTAR 
 simulate_forecast_vstar_bi <- function(vstar_obj,
                                        h0,h_1,H_1, s0,s_1,Y_1,
                                        horizon, delta) {
   regs <- vstar_obj$regs
   v    <- vstar_obj$v %||% c(1,0,0,0,0,0)
   zpick<- vstar_obj$z_var %||% "x_lag1"
   law_h <- function(h1,h2,H1, s1,s2,Y1) {
     # build one-row reg DF expected by vstar_eval
     df1 <- as.data.frame(setNames(list(h1,h2,H1,s1,s2,Y1),
                                   c("x_lag1","x_lag2","X_lag1","y_lag1","y_lag2","Y_lag1")))[, regs, drop = FALSE]
     zval <- if (zpick=="x_lag1") h1 else if (zpick=="y_lag1") s1 else
       sum(v * c(h1,h2,H1,s1,s2,Y1))
     pr   <- vstar_eval(vstar_obj$coef, df1, regs=regs, spec=vstar_obj$spec, z=zval)
     as.numeric(pr$yhat_h)
   }
   law_s <- function(h1,h2,H1, s1,s2,Y1) {
     df1 <- as.data.frame(setNames(list(h1,h2,H1,s1,s2,Y1),
                                   c("x_lag1","x_lag2","X_lag1","y_lag1","y_lag2","Y_lag1")))[, regs, drop = FALSE]
     zval <- if (zpick=="x_lag1") h1 else if (zpick=="y_lag1") s1 else
       sum(v * c(h1,h2,H1,s1,s2,Y1))
     pr   <- vstar_eval(vstar_obj$coef, df1, regs=regs, spec=vstar_obj$spec, z=zval)
     as.numeric(pr$yhat_s)
   }
   simulate_core_recursive_bi(h0,h_1,H_1, s0,s_1,Y_1, horizon, delta, law_h, law_s)
 }
 
 #  LASSO with separate coefficient vectors for each equation
 simulate_forecast_lasso_bi <- function(coef_h, coef_s,
                                        h0,h_1,H_1, s0,s_1,Y_1,
                                        horizon, delta) {
   coef_h <- coef_h[!is.na(coef_h)]; coef_s <- coef_s[!is.na(coef_s)]
   simulate_forecast_from_coefs_bi(coef_h, coef_s, h0,h_1,H_1, s0,s_1,Y_1, horizon, delta)
 }
 
 
 # 5) Run every simulation
 # (re)compute forecast origin & init
 if (!exists("df_samp_bi")) stop("df_samp_bi not in memory.")
 forecast_starter <- 3 * obs_per_year_real
 
 if (!exists("forecast_origin") || is.null(forecast_origin)) {
   if (1 + forecast_starter > nrow(df_samp_bi))
     stop("Not enough data after harmonization for the requested forecast_starter.")
   forecast_origin <- df_samp_bi$date[1 + forecast_starter]
 }
 
 init <- build_initial_objects_bi(df_samp_bi, origin_date = forecast_origin)
 
 # (optional sanity)
 print(forecast_origin)
 str(init)
 
  fc_list <- list()
 
 if (exists("m_ar2_h") && exists("m_ar2_s")) {
   fc_list$`AR(2)` <- simulate_forecast_from_models_bi(
     m_ar2_h, m_ar2_s,
     init$h0, init$h_1, init$H_1,
     init$s0, init$s_1, init$Y_1,
     horizon, delta
   )
 }
 if (exists("m_lin_h") && exists("m_lin_s")) {
   fc_list$Linear <- simulate_forecast_from_models_bi(
     m_lin_h, m_lin_s,
     init$h0, init$h_1, init$H_1,
     init$s0, init$s_1, init$Y_1,
     horizon, delta
   )
 }
 if (exists("m_min_h") && exists("m_min_s")) {
   fc_list$Minimal <- simulate_forecast_from_models_bi(
     m_min_h, m_min_s,
     init$h0, init$h_1, init$H_1,
     init$s0, init$s_1, init$Y_1,
     horizon, delta
   )
 }
 if (exists("m_int_h") && exists("m_int_s")) {
   fc_list$Intermediate <- simulate_forecast_from_models_bi(
     m_int_h, m_int_s,
     init$h0, init$h_1, init$H_1,
     init$s0, init$s_1, init$Y_1,
     horizon, delta
   )
 }
 if (exists("m_full_h") && exists("m_full_s")) {
   fc_list$Full <- simulate_forecast_from_models_bi(
     m_full_h, m_full_s,
     init$h0, init$h_1, init$H_1,
     init$s0, init$s_1, init$Y_1,
     horizon, delta
   )
 }
 # LASSO 
 if (exists("coef_lasso_h") && exists("coef_lasso_s")) {
   fc_list$LASSO <- simulate_forecast_lasso_bi(
     coef_lasso_h, coef_lasso_s,
     init$h0, init$h_1, init$H_1,
     init$s0, init$s_1, init$Y_1,
     horizon, delta
   )
 }
 
 # VSTAR family
 if (exists("vstar_std"))  fc_list$STAR_std   <- simulate_forecast_vstar_bi(vstar_std,
                                                                            init$h0, init$h_1, init$H_1, init$s0, init$s_1, init$Y_1, horizon, delta)
 if (exists("vstar_min"))  fc_list$STAR_min   <- simulate_forecast_vstar_bi(vstar_min,
                                                                            init$h0, init$h_1, init$H_1, init$s0, init$s_1, init$Y_1, horizon, delta)
 if (exists("vestar_min")) fc_list$ESTAR_min  <- simulate_forecast_vstar_bi(vestar_min,
                                                                            init$h0, init$h_1, init$H_1, init$s0, init$s_1, init$Y_1, horizon, delta)
 
 # 6) Stamp dates & trim to horizon_data
 origin_date <- as.Date(forecast_origin)
 
 from <- if (freq_real %in% c("q","m")) origin_date %m+% increment else origin_date + increment
 forecast_dates <- seq(from = from, by = obs_unit_real, length.out = horizon_data)
 
 trim_and_stamp_bi <- function(df, horizon_data, dates) {
   out <- df[1:horizon_data, , drop = FALSE]
   out$date <- dates
   out
 }
 fc_list <- lapply(fc_list, trim_and_stamp_bi, horizon_data = horizon_data, dates = forecast_dates)
 
 # Real data slice (for h and s)
 df_plot_h <- data_merged_real %>%
   dplyr::arrange(date) %>%
   dplyr::filter(date >= forecast_origin)
 if (nrow(df_plot_h) > horizon_data) df_plot_h <- df_plot_h[1:horizon_data, ]
 real_h <- df_plot_h %>%
   dplyr::select(date, h) %>%
   dplyr::rename(value = h) %>%
   dplyr::mutate(model = "Data")
 
 
 # Real data slice (for s)
 df_plot_s <- data_merged_fin %>%
   dplyr::arrange(date) %>%
   dplyr::filter(date >= forecast_origin)
 if (nrow(df_plot_s) > horizon_data) df_plot_s <- df_plot_s[1:horizon_data, ]
 real_s <- df_plot_s %>%
   dplyr::select(date, .data$s) %>%     
   dplyr::rename(value = .data$s) %>%
   dplyr::mutate(model = "Data")

 

 # 7) STATE-SPACE VISUALS (3D h-block & s-block; 2D h–H and s–Y)
 #    + h–s phase-portrait 
 for (model_name in names(fc_list)) {
   df <- fc_list[[model_name]]
   h <- df$h; s <- df$s; H <- df$H; Y <- df$Y
   n <- length(h)
   if (n < 3) next
   
   # 3D state space (h-block): (h_{t-2}, H_{t-1}, h_{t-1})
   state3d_h <- data.frame(h_tm2 = h[1:(n-2)], H_tm1 = H[2:(n-1)], h_tm1 = h[2:(n-1)])
   scatterplot3d::scatterplot3d(
     x = state3d_h$h_tm2, y = state3d_h$H_tm1, z = state3d_h$h_tm1,
     pch = 8, color = "black", angle = 20,
     main = paste0("State Space (h-block) — ", model_name),
     xlab = expression(h[t-2]), ylab = expression(H[t-1]), zlab = expression(h[t-1])
   )
   readline(prompt = "Press <Enter> for s-block 3D…")
   
   # 3D state space (s-block): (s_{t-2}, Y_{t-1}, s_{t-1})
   state3d_s <- data.frame(s_tm2 = s[1:(n-2)], Y_tm1 = Y[2:(n-1)], s_tm1 = s[2:(n-1)])
   scatterplot3d::scatterplot3d(
     x = state3d_s$s_tm2, y = state3d_s$Y_tm1, z = state3d_s$s_tm1,
     pch = 8, color = "black", angle = 20,
     main = paste0("State Space (s-block) — ", model_name),
     xlab = expression(s[t-2]), ylab = expression(Y[t-1]), zlab = expression(s[t-1])
   )
   readline(prompt = "Press <Enter> for (h_t, H_t) path…")
   
   # 2D (h_t, H_t)
   p_hH <- ggplot2::ggplot(data.frame(h_t = h, H_t = H), ggplot2::aes(h_t, H_t)) +
     ggplot2::geom_path(linewidth = 0.8) +
     ggplot2::geom_point(shape = 8, size = 1) +
     ggplot2::labs(title = paste0("Projection (h_t, H_t) — ", model_name),
                   x = expression(h[t]), y = expression(H[t])) +
     ggplot2::theme_minimal()
   print(p_hH)
   readline(prompt = "Press <Enter> for (s_t, Y_t) path…")
   
   # 2D (s_t, Y_t)
   p_sY <- ggplot2::ggplot(data.frame(s_t = s, Y_t = Y), ggplot2::aes(s_t, Y_t)) +
     ggplot2::geom_path(linewidth = 0.8) +
     ggplot2::geom_point(shape = 8, size = 1) +
     ggplot2::labs(title = paste0("Projection (s_t, Y_t) — ", model_name),
                   x = expression(s[t]), y = expression(Y[t])) +
     ggplot2::theme_minimal()
   print(p_sY)
   readline(prompt = "Press <Enter> for (h_t, s_t) phase portrait…")
   
   # ★ Requested: h–s phase portrait
   p_hs <- ggplot2::ggplot(data.frame(h_t = h, s_t = s), ggplot2::aes(h_t, s_t)) +
     ggplot2::geom_path(linewidth = 0.9) +
     ggplot2::geom_point(shape = 8, size = 1) +
     ggplot2::labs(title = paste0("Phase Portrait (h_t, s_t) — ", model_name),
                   x = expression(h[t]), y = expression(s[t])) +
     ggplot2::theme_minimal()
   print(p_hs)
   readline(prompt = "Press <Enter> for next model…")
 }
 

 
 # 8) SUPERPOSITION: 2D (h_t, H_t) for baseline vs STAR (optional buckets)
 all_fc <- fc_list
 base_models <- c("Linear","Minimal","Intermediate","Full","LASSO")
 star_models <- c("STAR_std","STAR_min","ESTAR_min")
 
 combined_proj_base <- dplyr::bind_rows(all_fc[intersect(names(all_fc), base_models)], .id = "model") |>
   dplyr::transmute(model, h_t = h, H_t = H)
 combined_proj_star <- dplyr::bind_rows(all_fc[intersect(names(all_fc), star_models)], .id = "model") |>
   dplyr::transmute(model, h_t = h, H_t = H)
 
 p_base <- ggplot2::ggplot(combined_proj_base, ggplot2::aes(h_t, H_t, color = model)) +
   ggplot2::geom_path(linewidth = 1) +
   ggplot2::geom_point(shape = 8, size = 1, alpha = 0.6) +
   ggplot2::scale_color_brewer("Model", palette = "Dark2") +
   ggplot2::labs(title = "Superposed (h_t, H_t) — Baseline",
                 x = expression(h[t]), y = expression(H[t])) +
   ggplot2::theme_minimal()
 
 p_star <- ggplot2::ggplot(combined_proj_star, ggplot2::aes(h_t, H_t, color = model)) +
   ggplot2::geom_path(linewidth = 1) +
   ggplot2::geom_point(shape = 8, size = 1, alpha = 0.6) +
   ggplot2::scale_color_brewer("Model", palette = "Set1") +
   ggplot2::labs(title = "Superposed (h_t, H_t) — STAR family",
                 x = expression(h[t]), y = expression(H[t])) +
   ggplot2::theme_minimal()
 
 gridExtra::grid.arrange(p_base, p_star, ncol = 2)
 
 
 
 

 # SUPERPOSITION: 2D (s_t, Y_t) for baseline vs STAR (optional buckets)
 all_fc <- fc_list
 
 # choose which models to show
 base_models_s <- c("Linear","Minimal","Intermediate","LASSO")
 star_models_s <- c("STAR_std","STAR_min","ESTAR_min")
 
 keep_base_s <- intersect(names(fc_list), base_models_s)
 keep_star_s <- intersect(names(fc_list), star_models_s)
 
 combined_proj_base_s <- if (length(keep_base_s)) {
   dplyr::bind_rows(fc_list[keep_base_s], .id = "model") |>
     dplyr::transmute(model, s_t = s, Y_t = Y)
 } else NULL
 
 combined_proj_star_s <- if (length(keep_star_s)) {
   dplyr::bind_rows(fc_list[keep_star_s], .id = "model") |>
     dplyr::transmute(model, s_t = s, Y_t = Y)
 } else NULL
 
 if (!is.null(combined_proj_base_s)) {
   p_base_s <- ggplot2::ggplot(combined_proj_base_s, ggplot2::aes(s_t, Y_t, color = model)) +
     ggplot2::geom_path(linewidth = 1) +
     ggplot2::geom_point(shape = 8, size = 1, alpha = 0.6) +
     ggplot2::scale_color_brewer("Model", palette = "Dark2") +
     ggplot2::labs(title = "Superposed (s_t, Y_t) — Baseline",
                   x = expression(s[t]), y = expression(Y[t])) +
     ggplot2::theme_minimal()
   print(p_base_s)
 }
 
 if (!is.null(combined_proj_star_s)) {
   p_star_s <- ggplot2::ggplot(combined_proj_star_s, ggplot2::aes(s_t, Y_t, color = model)) +
     ggplot2::geom_path(linewidth = 1) +
     ggplot2::geom_point(shape = 8, size = 1, alpha = 0.6) +
     ggplot2::scale_color_brewer("Model", palette = "Set1") +
     ggplot2::labs(title = "Superposed (s_t, Y_t) — STAR family",
                   x = expression(s[t]), y = expression(Y[t])) +
     ggplot2::theme_minimal()
   print(p_star_s)
 }
 
 
 gridExtra::grid.arrange(p_base_s, p_star_s, ncol = 2)
 
 
 
 

 #
 # SUPERPOSITION: CROSSED PHASE PORTRAITS 2D (h_t, s_t) 
 # 
 
 all_fc <- fc_list
 base_models_hs <- c("Linear","Minimal","Intermediate","Full","LASSO")
 star_models_hs <- c("STAR_std","STAR_min","ESTAR_min")
 
 keep_base_hs <- intersect(names(all_fc), base_models_hs)
 keep_star_hs <- intersect(names(all_fc), star_models_hs)
 
 combined_base_hs <- if (length(keep_base_hs)) {
   dplyr::bind_rows(all_fc[keep_base_hs], .id = "model") |>
     dplyr::transmute(model, h_t = h, s_t = s) |>
     dplyr::filter(is.finite(h_t), is.finite(s_t))
 } else NULL
 
 combined_star_hs <- if (length(keep_star_hs)) {
   dplyr::bind_rows(all_fc[keep_star_hs], .id = "model") |>
     dplyr::transmute(model, h_t = h, s_t = s) |>
     dplyr::filter(is.finite(h_t), is.finite(s_t))
 } else NULL
 
 p_base_hs <- NULL; p_star_hs <- NULL
 if (!is.null(combined_base_hs)) {
   p_base_hs <- ggplot2::ggplot(combined_base_hs, ggplot2::aes(h_t, s_t, color = model)) +
     ggplot2::geom_path(linewidth = 1) +
     ggplot2::geom_point(shape = 8, size = 1, alpha = 0.6) +
     ggplot2::scale_color_brewer("Model", palette = "Dark2") +
     ggplot2::labs(title = "Superposed (h_t, s_t) — Baseline",
                   x = expression(h[t]), y = expression(s[t])) +
     ggplot2::theme_minimal()
 }
 if (!is.null(combined_star_hs)) {
   p_star_hs <- ggplot2::ggplot(combined_star_hs, ggplot2::aes(h_t, s_t, color = model)) +
     ggplot2::geom_path(linewidth = 1) +
     ggplot2::geom_point(shape = 8, size = 1, alpha = 0.6) +
     ggplot2::scale_color_brewer("Model", palette = "Set1") +
     ggplot2::labs(title = "Superposed (h_t, s_t) — STAR family",
                   x = expression(h[t]), y = expression(s[t])) +
     ggplot2::theme_minimal()
 }
 
 if (!is.null(p_base_hs) && !is.null(p_star_hs)) {
   gridExtra::grid.arrange(p_base_hs, p_star_hs, ncol = 2)
 } else {
   if (!is.null(p_base_hs)) print(p_base_hs)
   if (!is.null(p_star_hs)) print(p_star_hs)
 }
 
 # ─
 # SUPERPOSITION: CROSSED PHASE PORTRAITS OF SLOW STOCKS 2D (H_t, Y_t) 
 # 
 
 base_models_HY <- c("Linear","Minimal","Intermediate","Full","LASSO")
 star_models_HY <- c("STAR_std","STAR_min","ESTAR_min")
 
 keep_base_HY <- intersect(names(all_fc), base_models_HY)
 keep_star_HY <- intersect(names(all_fc), star_models_HY)
 
 combined_base_HY <- if (length(keep_base_HY)) {
   dplyr::bind_rows(all_fc[keep_base_HY], .id = "model") |>
     dplyr::transmute(model, H_t = H, Y_t = Y) |>
     dplyr::filter(is.finite(H_t), is.finite(Y_t))
 } else NULL
 
 combined_star_HY <- if (length(keep_star_HY)) {
   dplyr::bind_rows(all_fc[keep_star_HY], .id = "model") |>
     dplyr::transmute(model, H_t = H, Y_t = Y) |>
     dplyr::filter(is.finite(H_t), is.finite(Y_t))
 } else NULL
 
 p_base_HY <- NULL; p_star_HY <- NULL
 if (!is.null(combined_base_HY)) {
   p_base_HY <- ggplot2::ggplot(combined_base_HY, ggplot2::aes(H_t, Y_t, color = model)) +
     ggplot2::geom_path(linewidth = 1) +
     ggplot2::geom_point(shape = 8, size = 1, alpha = 0.6) +
     ggplot2::scale_color_brewer("Model", palette = "Dark2") +
     ggplot2::labs(title = "Superposed (H_t, Y_t) — Baseline",
                   x = expression(H[t]), y = expression(Y[t])) +
     ggplot2::theme_minimal()
 }
 if (!is.null(combined_star_HY)) {
   p_star_HY <- ggplot2::ggplot(combined_star_HY, ggplot2::aes(H_t, Y_t, color = model)) +
     ggplot2::geom_path(linewidth = 1) +
     ggplot2::geom_point(shape = 8, size = 1, alpha = 0.6) +
     ggplot2::scale_color_brewer("Model", palette = "Set1") +
     ggplot2::labs(title = "Superposed (H_t, Y_t) — STAR family",
                   x = expression(H[t]), y = expression(Y[t])) +
     ggplot2::theme_minimal()
 }
 
 if (!is.null(p_base_HY) && !is.null(p_star_HY)) {
   gridExtra::grid.arrange(p_base_HY, p_star_HY, ncol = 2)
 } else {
   if (!is.null(p_base_HY)) print(p_base_HY)
   if (!is.null(p_star_HY)) print(p_star_HY)
 }
 
 
 
 
 
 # SELECTORS 
 # Time-series superpositions (h)
 models_h_panel_a <- c("AR(2)","Minimal","Linear")
 models_h_panel_b <- c("Intermediate","Full","LASSO")
 models_h_star    <- c("STAR_std","STAR_min","ESTAR_min")
 
 # Time-series superpositions (s)
 models_s_panel_a <- c("AR(2)","Minimal","Linear")
 models_s_panel_b <- c("Intermediate","LASSO", "Full")
 models_s_star    <- c("STAR_std","STAR_min","ESTAR_min")
 
 pick_available <- function(wanted, pool) intersect(wanted, names(pool))
 
 filter_long_by_models <- function(long_df, wanted, pool) {
   keep <- pick_available(wanted, pool)
   dplyr::filter(long_df, model %in% keep)
 }
 
 bind_phase <- function(pool, wanted, xy = c("h","H")) {
   keep <- pick_available(wanted, pool)
   stopifnot(length(xy) == 2)
   nm_x <- xy[1]; nm_y <- xy[2]
   dplyr::bind_rows(pool[keep], .id = "model") |>
     dplyr::transmute(model, x = .data[[nm_x]], y = .data[[nm_y]])
 }

 
 # 9) SERIES TIME PLOTS - h
 # Long frames for h
 fc_long_h <- lapply(names(fc_list), function(nm) {
   fc_list[[nm]] |>
     dplyr::select(date, h) |>
     dplyr::rename(value = h) |>
     dplyr::mutate(model = nm)
 }) |> dplyr::bind_rows()
 
 model_cols <- c(
   Data = "black",
   "AR(2)" = "#1f78b4",
   Minimal = "#33a02c",
   Linear  = "#fb9a99",
   Intermediate = "#e31a1c",
   Full    = "#ff7f00",
   LASSO   = "#6a3d9a",
   STAR_std  = "#008b8b",
   STAR_min  = "#b15928",
   ESTAR_min = "#a6cee3"
 )
 
 # Example split panels (adjust subsets if some models are missing)
 panel_a <- ggplot2::ggplot(
   dplyr::bind_rows(real_h, filter_long_by_models(fc_long_h, models_h_panel_a, fc_list)),
   ggplot2::aes(date, value, color = model)
 ) + ggplot2::geom_line(linewidth = 1) +
   ggplot2::scale_color_manual(values = model_cols) +
   ggplot2::labs(title = "Panel (a): Data + selected (h)", x = "Date", y = "h") +
   ggplot2::theme_minimal()
 
 panel_b <- ggplot2::ggplot(
   dplyr::bind_rows(real_h, filter_long_by_models(fc_long_h, models_h_panel_b, fc_list)),
   ggplot2::aes(date, value, color = model)
 ) + ggplot2::geom_line(linewidth = 1) +
   ggplot2::scale_color_manual(values = model_cols) +
   ggplot2::labs(title = "Panel (b): Data + selected (h)", x = "Date", y = "h") +
   ggplot2::theme_minimal()
 
 gridExtra::grid.arrange(panel_a, panel_b, ncol = 2)
 
 # STAR-only panel (h)
 star_only_h <- dplyr::bind_rows(
   real_h, filter_long_by_models(fc_long_h, models_h_star, fc_list)
 )
 p_star_h <- ggplot2::ggplot(star_only_h, ggplot2::aes(date, value, color = model)) +
   ggplot2::geom_line(linewidth = 1) +
   ggplot2::scale_color_manual(values = model_cols[names(model_cols) %in% c("Data", models_h_star)]) +
   ggplot2::labs(title = "STAR (h)", x = "Date", y = "h", color = NULL) +
   ggplot2::theme_minimal()
 
 p_star_h
 
 # (Optional) overlay: Data + AR(2) + Minimal (h)
 if (all(c("AR(2)","Minimal") %in% names(fc_list))) {
   combo_h <- dplyr::bind_rows(
     real_h,
     dplyr::filter(fc_long_h, model %in% c("AR(2)","Minimal"))
   )
   p_combo <- ggplot2::ggplot(combo_h, ggplot2::aes(date, value, color = model)) +
     ggplot2::geom_line(linewidth = 1) +
     ggplot2::labs(title = "Real Data with AR(2) and Minimal Forecasts (h)",
                   x = "Date", y = "h (cycle)") +
     ggplot2::theme_minimal() +
     ggplot2::scale_color_manual(values = c("Data"="black","AR(2)"="blue","Minimal"="red"))
   print(p_combo)
 }
 
 

 
 # SERIES TIME PLOTS — s 
 
 # 0) Build the real series frame for s (Data)
 real_s <- if ("s" %in% names(data_merged_real)) {
   data_merged_real %>%
     dplyr::arrange(date) %>%
     dplyr::filter(date >= forecast_origin) %>%
     { if (exists("horizon_data")) dplyr::slice_head(., n = horizon_data) else . } %>%
     dplyr::select(date, s) %>%
     dplyr::rename(value = s) %>%
     dplyr::mutate(model = "Data")
 } else {
   # fallback if the merged panel has no 's'
   df_samp_bi %>%
     dplyr::arrange(date) %>%
     dplyr::filter(date >= forecast_origin) %>%
     { if (exists("horizon_data")) dplyr::slice_head(., n = horizon_data) else . } %>%
     dplyr::select(date, s) %>%
     dplyr::rename(value = s) %>%
     dplyr::mutate(model = "Data")
 }
 
 # 1) Long forecast frame for s
 fc_long_s <- lapply(names(fc_list), function(nm) {
   this <- fc_list[[nm]]
   if (!all(c("date","s") %in% names(this))) return(NULL)
   this |>
     dplyr::select(date, s) |>
     dplyr::rename(value = s) |>
     dplyr::mutate(model = nm)
 }) |>
   purrr::compact() |>
   dplyr::bind_rows()
 
 # 2) Colors 
 model_cols <- c(
   Data = "black",
   "AR(2)" = "#1f78b4",
   Minimal = "#33a02c",
   Linear  = "#fb9a99",
   Intermediate = "#e31a1c",
   Full    = "#ff7f00",
   LASSO   = "#6a3d9a",
   STAR_std  = "#008b8b",
   STAR_min  = "#b15928",
   ESTAR_min = "#a6cee3"
 )
 
 # 3) Split panels for s
 panel_a_s <- ggplot2::ggplot(
   dplyr::bind_rows(real_s, filter_long_by_models(fc_long_s, models_s_panel_a, fc_list)),
   ggplot2::aes(date, value, color = model)
 ) + ggplot2::geom_line(linewidth = 1) +
   ggplot2::scale_color_manual(values = model_cols) +
   ggplot2::labs(title = "Panel (a): Data + selected (s)", x = "Date", y = "s") +
   ggplot2::theme_minimal()
 
 panel_b_s <- ggplot2::ggplot(
   dplyr::bind_rows(real_s, filter_long_by_models(fc_long_s, models_s_panel_b, fc_list)),
   ggplot2::aes(date, value, color = model)
 ) + ggplot2::geom_line(linewidth = 1) +
   ggplot2::scale_color_manual(values = model_cols) +
   ggplot2::labs(title = "Panel (b): Data + selected (s)", x = "Date", y = "s") +
   ggplot2::theme_minimal()
 
 gridExtra::grid.arrange(panel_a_s, panel_b_s, ncol = 2)
 
 
 
 # 4) STAR-only panel (s)
 star_only_s <- dplyr::bind_rows(
   real_s, filter_long_by_models(fc_long_s, models_s_star, fc_list)
 )
 p_star_s <- ggplot2::ggplot(star_only_s, ggplot2::aes(date, value, color = model)) +
   ggplot2::geom_line(linewidth = 1) +
   ggplot2::scale_color_manual(values = model_cols[names(model_cols) %in% c("Data", models_s_star)]) +
   ggplot2::labs(title = "STAR (s)", x = "Date", y = "s", color = NULL) +
   ggplot2::theme_minimal()
 
 p_star_s
 
 # 5) Optional overlay: Data + AR(2) + Minimal (s)
 if (all(c("AR(2)","Minimal") %in% names(fc_list))) {
   combo_s <- dplyr::bind_rows(
     real_s,
     dplyr::filter(fc_long_s, model %in% c("AR(2)","Minimal"))
   )
   p_combo_s <- ggplot2::ggplot(combo_s, ggplot2::aes(date, value, color = model)) +
     ggplot2::geom_line(linewidth = 1) +
     ggplot2::labs(title = "Real Data with AR(2) and Minimal Forecasts (s)",
                   x = "Date", y = "s (cycle)") +
     ggplot2::theme_minimal() +
     ggplot2::scale_color_manual(values = c("Data"="black","AR(2)"="blue","Minimal"="red"))
   print(p_combo_s)
 }
 