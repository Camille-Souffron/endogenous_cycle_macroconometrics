################################################################################
# CHAPTER 3 - ENDOUGENOUS-DGP TEST
################################################################################

# Camille Souffron - MASTER THESIS APE (PSE & ENS)

#NB: all packages in the DATA_LOADING file


####### 1. HELPER FUNCTIONS
serial_pass <- function(VARobj, lags_pt = 1, alpha = 0.10) {
  pv <- suppressWarnings(tryCatch({
    st <- serial.test(VARobj, lags.pt = lags_pt, type = "BG")
    candidates <- c(
      tryCatch(as.numeric(st$serial$p.value), error=function(e) NA_real_),
      tryCatch(as.numeric(st$serial$BG$p.value), error=function(e) NA_real_),
      tryCatch(as.numeric(st$serial[["BG test"]]$p.value),error=function(e) NA_real_)
    )
    cand <- candidates[is.finite(candidates)]
    if (length(cand) >= 1L) cand[1L] else NA_real_
  }, error=function(e) NA_real_))
  if (!is.finite(pv)) return(TRUE)
  pv > alpha
}


cycle_length <- function(lambda, tol_im = 1e-10, tol_ang = 1e-12){
  re <- Re(lambda); im <- Im(lambda)
  if (abs(im) < tol_im) return(Inf)
  theta <- atan2(abs(im), re)
  if (theta < tol_ang) return(Inf)
  2*pi/theta
}


cycle_length_acos <- function(lambda, tol_im = 1e-10, tol_ang = 1e-12){
  if (abs(Im(lambda)) < tol_im) return(Inf)
  x <- Re(lambda)/Mod(lambda); x <- max(-1, min(1, x))
  theta <- acos(x)
  if (theta < tol_ang) return(Inf)
  2*pi/theta
}


## a2*b1 test for VAR (delta-method using joint vcov from vars::VAR)
(test_cycle_condition <- function(VARobj, alpha=0.10){
  ey <- VARobj$varresult$y; ef <- VARobj$varresult$f
  a2 <- unname(coef(ey)["f.l1"]); b1 <- unname(coef(ef)["y.l1"]) ; est <- a2*b1
  if (!all(is.finite(c(a2,b1)))) return(list(stat=est,se=NA_real_,t=NA_real_,pval=NA_real_))
  V <- vcov(VARobj); rn <- rownames(V)
  pick <- function(pats){ i <- unique(unlist(lapply(pats, function(p) grep(p,rn,perl=TRUE)))) ; if (length(i)==1L) i else NA_integer_ }
  i_a2 <- pick(c("^y[.:]f\\.l1$")); i_b1 <- pick(c("^f[.:]y\\.l1$"))
  if (anyNA(c(i_a2,i_b1))) return(list(stat=est,se=NA_real_,t=NA_real_,pval=NA_real_))
  var_a2 <- as.numeric(V[i_a2,i_a2]); var_b1 <- as.numeric(V[i_b1,i_b1]); cov_ab <- as.numeric(V[i_a2,i_b1])
  se2 <- (b1^2)*var_a2 + (a2^2)*var_b1 + 2*a2*b1*cov_ab
  if (!is.finite(se2) || se2<0) return(list(stat=est,se=NA_real_,t=NA_real_,pval=NA_real_))
  se <- sqrt(se2); if (!is.finite(se) || se==0) return(list(stat=est,se=se,t=NA_real_,pval=1.0))
  tval <- est/se; pval <- pt(tval, df = ey$df.residual, lower.tail = TRUE)
  list(stat=est, se=se, t=tval, pval=pval)
})

## 1.4 TRUE a2*b1 for each DGP (off-diagonal product of the first-lag Jacobian at the steady state)
true_a2b1 <- function(DGP){
  switch(DGP,
         "SM39" = {
           # Matches simulate_SM39() defaults
           c <- 0.4; beta <- 2
           c * beta * (c - 1)
         },
         "Metzler" = {
           # From linearized mapping: Y_t = (b+kb)Y_{t-1} - S_{t-1}; S_t = kb Y_{t-1} - S_{t-1}
           b <- 0.2; k <- 0.99
           -k * b
         },
         "BGP20" = {
           alpha1 <- 0.4; delta <- 0.2
           -alpha1 * (1 - delta)
         },
         "Goodwin" = { 
           # Goodwin form; y = prey x, f = predator y
           r <- 3.1; c <- 2.1
           xstar <- 1/c
           ystar <- r*(c-1)/c^2
           a2 <- -c * xstar      
           b1 <-  c * ystar     
           a2 * b1             
         },
         "Puu" = {
           a <- 1.1; b <- 0.1
           # X_{t+1} = X_t + Z_t, Z_{t+1} = a Z_t - (a+1)Z_t^3 - b X_t ; at (0,0): a2=1, b1=-b
           -b
         },
         "Henon" = {
           a <- 1.4; b <- 0.3
           # x_t = 1 - a x_{t-1}^2 + y_{t-1}; y_t = b x_{t-1}  ⇒ a2=1, b1=b
           1 * b                  # > 0 ⇒ ‘false positive’ check
         },
         stop("Unknown DGP: ", DGP)
  )
}




####### 2.  ENDOGENOUS-CYCLE DATA‑GENERATING PROCESSES

# helper: uniform jitter around steady‑state
jitter_ss <- function(ss, width = 0.05) ss * runif(1, 1-width, 1+width)

# 2·1  Samuelson (1939) SM39 
simulate_SM39 <- function(T, c = 0.4, beta = 2, G = 0, # 1.99 rather than 1.4, for persistence
                          shock_rho = 0.8, sd_eps = 1, eps0 = 0.1){
  # Steady state: C* = I* = 0
  C <- I <- numeric(T)
  uy <- uf <- numeric(T)
  # random init in [-eps0, +eps0]
  C[1] <- runif(1, -eps0, eps0)
  I[1] <- runif(1, -eps0, eps0)
  for(t in 2:T){
    uy[t] <- shock_rho*uy[t-1] + rnorm(1, 0, sd_eps)
    uf[t] <- shock_rho*uf[t-1] + rnorm(1, 0, sd_eps)
    C[t] <- c*(C[t-1] + I[t-1] + G) + uy[t]
    I[t] <- beta*( c*(C[t-1] + I[t-1] + G) - C[t-1] ) + uf[t]
  }
  data.frame(C = C, I = I)
}


# 2.2  Beaudry–Galizia–Portier (2020) BGP20
simulate_BGP20 <- function(
    T,
    alpha0 = -1,  # -1 normal
    alpha1 = 0.4,  # 0.4 normal ; 3 (high acc dislike) if want chaos
    alpha2 = 0.9, # 0.9 normal
    alpha3 = 2, # 2 normal ; small if want focus (<=1.4)
    delta  = 0.2, # 0.2 normal 
    shock_rho = 0.8, 
    sd_eps = 1, 
    eps0 = 1
){
  # Steady state: I* = X* = 0 (for these parameters)
  I <- X <- numeric(T)
  uI <- uX <- numeric(T)
  # random init in [-eps0, +eps0]
  I[1] <- runif(1, -eps0, eps0)
  X[1] <- runif(1, -eps0, eps0)
  for(t in 2:T){
    # AR(1) shocks if desired (else sd_eps=0 for deterministic)
    uI[t] <- shock_rho*uI[t-1] + rnorm(1, 0, sd_eps)
    uX[t] <- shock_rho*uX[t-1] + rnorm(1, 0, sd_eps)
    # Main model equations
    I[t] <- alpha0 - alpha1*(1-delta)*X[t-1] + (alpha2-alpha1)*I[t-1] +
      alpha3/(1+exp(-I[t-1])) + uI[t]
    X[t] <- (1-delta)*X[t-1] + I[t-1] + uX[t]
  }
  data.frame(I = I, X = X)
}


# Simulate 200 periods, default parameters, no noise
sim <- simulate_BGP20(150, sd_eps =0, shock_rho = 0)

# Plot the time series
matplot(sim, type = "l", lty = 1, col = c("blue", "red"), ylab = "Value", xlab = "Time", main = "BGP(2020) simulation")
legend("topright", legend = c("I", "X"), col = c("blue", "red"), lty = 1)



# 2.3 Adapted Goodwin (1967) - modified Lotka–Volterra

simulate_Goodwin <- function(T,
                          r = 3.1, c = 2.1,
                          shock_rho = 0.8, sd_eps = 0.01) {
  
  x <- y <- numeric(T)
  # coexistence steady‑state 
  ssx <- 1/c;  ssy <- r*(c-1)/c^2
  x[1] <- 0.6;   y[1] <- 0.6
  
  ux <- uy <- numeric(T)
  for (t in 2:T) {
    ux[t] <- shock_rho*ux[t-1] + rnorm(1, 0, sd_eps)
    uy[t] <- shock_rho*uy[t-1] + rnorm(1, 0, sd_eps)
    
    y[t] <- c*x[t-1]*y[t-1]                      + uy[t]     # predator 
    x[t] <- (r+1)*x[t-1] - r*x[t-1]^2 -
      c*x[t-1]*y[t-1]                     + ux[t]     # prey 
  }
  data.frame(y = x, f = y)           # y ≡ prey , f ≡ predator 
}




# 2.4  Puu (2002) nonlinear accelerator‑multiplier
simulate_Puu <- function(
    T,
    a         = 1.1,     # net accelerator parameter (v - s) # 2 = Chaotic, 1.1 = limit cycle, <1 = focus
    b         = 0.1,   # saving parameter b = (1 - ε) s
    shock_rho = 0.8,     # AR(1) persistence of shocks
    sd_eps    = 0.1,   # shock volatility
    eps0      = 1      # initial draw range
) {
  # allocate
  X <- numeric(T)
  Z <- numeric(T)
  u <- numeric(T)
  v <- numeric(T)
  
  # initial states
  X[1] <- runif(1, -eps0, eps0)
  Z[1] <- runif(1, -eps0, eps0)
  u[1] <- 0
  v[1] <- 0
  
  for (t in 1:(T-1)) {
    # AR(1) shocks
    u[t+1] <- shock_rho * u[t] + rnorm(1, 0, sd_eps)
    v[t+1] <- shock_rho * v[t] + rnorm(1, 0, sd_eps)
    
    # Puu-Sushko cubic multiplier–accelerator map with shocks
    X[t+1] <- X[t] + Z[t] + u[t+1]
    Z[t+1] <- a * Z[t] - (a + 1) * Z[t]^3 - b * X[t] + v[t+1]
  }
  
  data.frame(X = X, Z = Z)
}


# run
sim_Puu_test  <- simulate_Puu(T = 100, a=1.1, b=0.1, shock_rho = 0.8, sd_eps = 0.1) #max rho 0.4 and sd 0.1

# Plot the time series
matplot(sim_Puu_test, type = "l", lty = 1, col = c("blue", "red"), ylab = "Value", xlab = "Time", main = "Puu (1990")
legend("topright", legend = c("X", "Z"), col = c("blue", "red"), lty = 1)




# 2.5 Hénon chaotic map  
simulate_Henon <- function(T,
                           a = 1.4, b = 0.3,
                           shock_rho = 0.8, sd_eps = 0.01) {
  # pre‐allocate
  x <- numeric(T)
  y <- numeric(T)
  
  # AR(1) shocks for x‐ and y‐equations
  ux <- uy <- numeric(T)
  
  # initial condition
  x[1] <- 0.1
  y[1] <- 0.1
  
  for (t in 2:T) {
    # shock updates
    ux[t] <- shock_rho * ux[t-1] + rnorm(1, 0, sd_eps)
    uy[t] <- shock_rho * uy[t-1] + rnorm(1, 0, sd_eps)
    
    # Hénon map + additive shocks
    x[t] <- 1 - a * x[t-1]^2 + y[t-1] + ux[t]
    y[t] <- b * x[t-1] + uy[t]
  }
  
  data.frame(x = x, y = y)
}

# run
sim_Henon_test  <- simulate_Henon(T = 100,
                                      a = 1.4, b = 0.3,
                                      shock_rho = 0.3, sd_eps = 0.01)
# Plot the time series
matplot(sim_Henon_test, type = "l", lty = 1, col = c("blue", "red"), ylab = "Value", xlab = "Time", main = "Henon")
legend("topright", legend = c("X", "Z"), col = c("blue", "red"), lty = 1)












####### 3. FULL RUN (ALL MODELS) 

simulate_data <- function(DGP,
                          T,
                          shock_rho = NULL,
                          sd_eps    = NULL,
                          eps0      = NULL) {
  fn <- switch(DGP,
               "SM39"    = simulate_SM39,
               "Metzler" = simulate_Metzler,
               "BGP20"   = simulate_BGP20,
               "Goodwin" = simulate_Goodwin,
               "Puu"     = simulate_Puu,
               "Henon"   = simulate_Henon,
               stop("Unknown DGP: ", DGP)
  )
  ## always pass T
  args <- list(T = T)
  ## only override if non‐NULL
  if (!is.null(shock_rho)) args$shock_rho <- shock_rho
  if (!is.null(sd_eps))    args$sd_eps    <- sd_eps
  if (!is.null(eps0))      args$eps0      <- eps0
  do.call(fn, args)
}







##### 4. STAR (bivariate, common transition) 
## Parameterizations:
##  - standard_lstar: y_t = a0y + a1y*G + (b0y + b1y*G)' X;  f_t analogous
##  - minimal_lstar : y_t = a0y + (b0y + d_y*G)' X;           f_t analogous
##  - minimal_estar : same as minimal_lstar but G=1-exp(-gamma*(z-c)^2)
## Here X = (y_{t-1}, f_{t-1})'. Transition z = q'X (optionally standardized).

.star_G_logistic <- function(z, eta, c){ gamma <- exp(eta); plogis(gamma*(z-c)) }
.star_dGdz_logistic <- function(z, eta, c){ gamma <- exp(eta); g <- plogis(gamma*(z-c)); gamma*g*(1-g) }
.star_G_estar <- function(z, eta, c){ gamma <- exp(eta); 1 - exp(-gamma*(z-c)^2) }
.star_dGdz_estar <- function(z, eta, c){ gamma <- exp(eta); 2*gamma*(z-c)*exp(-gamma*(z-c)^2) }

.star_build_X <- function(df){
  stopifnot(all(c("y","f","y_l1","f_l1") %in% names(df)))
  list(y = as.numeric(df$y), f = as.numeric(df$f), X = cbind(y_l1=df$y_l1, f_l1=df$f_l1))
}

## STAR evaluated means given theta and data matrices
.star_eval <- function(theta, X, spec, z_std, q, names_out=FALSE){
  n <- nrow(X); R <- ncol(X)
  idx <- 0L
  if (spec=="standard_lstar"){
    a0y <- theta[idx<-idx+1]; a1y <- theta[idx<-idx+1]
    b0y <- theta[(idx+1):(idx+R)]; idx<-idx+R
    b1y <- theta[(idx+1):(idx+R)]; idx<-idx+R
    a0f <- theta[idx<-idx+1]; a1f <- theta[idx<-idx+1]
    b0f <- theta[(idx+1):(idx+R)]; idx<-idx+R
    b1f <- theta[(idx+1):(idx+R)]; idx<-idx+R
    eta <- theta[idx<-idx+1]; c0 <- theta[idx<-idx+1]
    G   <- .star_G_logistic(z_std, eta, c0)
    yhat <- a0y + a1y*G + rowSums( matrix(rep(b0y,each=n),n)*X + matrix(rep(b1y,each=n),n)*X*G )
    fhat <- a0f + a1f*G + rowSums( matrix(rep(b0f,each=n),n)*X + matrix(rep(b1f,each=n),n)*X*G )
  } else {
    a0y <- theta[idx<-idx+1]
    b0y <- theta[(idx+1):(idx+R)]; idx<-idx+R
    dy  <- theta[(idx+1):(idx+R)]; idx<-idx+R
    a0f <- theta[idx<-idx+1]
    b0f <- theta[(idx+1):(idx+R)]; idx<-idx+R
    df  <- theta[(idx+1):(idx+R)]; idx<-idx+R
    eta <- theta[idx<-idx+1]; c0 <- theta[idx<-idx+1]
    G   <- if (spec=="minimal_lstar") .star_G_logistic(z_std, eta, c0) else .star_G_estar(z_std, eta, c0)
    yhat <- a0y + rowSums( matrix(rep(b0y,each=n),n)*X + matrix(rep(dy,each=n),n)*X*G )
    fhat <- a0f + rowSums( matrix(rep(b0f,each=n),n)*X + matrix(rep(df,each=n),n)*X*G )
  }
  if (isTRUE(names_out)) return(list(yhat=yhat,fhat=fhat, names=names(theta)))
  list(yhat=yhat,fhat=fhat)
}


## Feasible GLS for STAR with AR(P) errors (Cochrane–Orcutt)
.fit_eq_gls_AR <- function(y, X, ar_order = 1L, maxit = 5L, tol = 1e-6){
  # OLS start
  ols  <- stats::lm.fit(X, y)
  beta <- as.numeric(ols$coefficients)
  res  <- y - drop(X %*% beta)
  p    <- max(0L, as.integer(ar_order))
  phi  <- if (p > 0L) rep(0, p) else numeric(0)
  
  # fast AR(p) on residuals: Yule–Walker (AR(1) closed form)
  yw_phi <- function(e){
    if (p == 0L) return(numeric(0))
    if (p == 1L) {
      num <- sum(e[-1L] * e[-length(e)])
      den <- sum(e[-length(e)]^2)
      return(if (den == 0) 0 else as.numeric(num/den))
    } else {
      fit <- try(stats::ar.yw(e, order.max = p, aic = FALSE, demean = FALSE),
                 silent = TRUE)
      if (inherits(fit, "try-error") || is.null(fit$ar)) rep(0, p) else as.numeric(fit$ar)
    }
  }
  
  if (p > 0L){
    for (it in 1:maxit){
      phi_new <- yw_phi(res)
      tr      <- .gls_transform_AR(y, X, phi_new)
      gls     <- stats::lm.fit(tr$X, tr$y)
      beta_new<- as.numeric(gls$coefficients)
      res_new <- tr$y - drop(tr$X %*% beta_new)
      if (max(abs(beta_new - beta), abs(phi_new - phi)) < tol){
        beta <- beta_new; phi <- phi_new; res <- res_new
        break
      }
      beta <- beta_new; phi <- phi_new; res <- res_new
    }
    tr  <- .gls_transform_AR(y, X, phi)       # final transform
    X_  <- tr$X; y_ <- tr$y
    gls <- stats::lm.fit(X_, y_)
    drop_p <- tr$drop
  } else {
    X_ <- X; y_ <- y; gls <- ols; drop_p <- 0L
  }
  
  df <- nrow(X_) - ncol(X_)
  s2 <- sum((y_ - drop(X_ %*% gls$coefficients))^2) / max(df, 1L)
  V  <- tryCatch(s2 * chol2inv(chol(crossprod(X_))),
                 error = function(e) matrix(NA_real_, ncol = ncol(X_), nrow = ncol(X_)))
  list(beta = as.numeric(gls$coefficients),
       phi = phi, s2 = s2, V = V, drop = drop_p)
}



## Fit by joint SSE minimization with robust starts; returns theta, vcov, and data needed for delta-method
estimate_star_bivar <- function(df, spec=c("standard_lstar","minimal_lstar","minimal_estar"),
                                z_choice=c("y","f","combo"), q=c(1,0),
                                starts_gamma=c(0.5,1,2,5,10), starts_c=c("median","q25","q75"),
                                maxit=2000, trace=TRUE, compute_vcov=TRUE){
  spec <- match.arg(spec); z_choice <- match.arg(z_choice)
  stopifnot(all(c("y","f","y_l1","f_l1") %in% names(df)))
  d <- .star_build_X(df); y <- d$y; f <- d$f; X <- d$X; n <- nrow(X); R <- ncol(X)
  ## Build z = q'X and standardize for numerical stability
  if (z_choice=="y") q <- c(1,0) else if (z_choice=="f") q <- c(0,1) else { stopifnot(length(q)==R) }
  z_raw <- drop(X %*% matrix(q, ncol=1))
  z_mu <- stats::median(z_raw, na.rm=TRUE); z_sd <- stats::sd(z_raw, na.rm=TRUE)
  if (!is.finite(z_sd) || z_sd==0) z_sd <- stats::mad(z_raw, na.rm=TRUE)
  z_std <- (z_raw - z_mu)/z_sd
  
  ## Warm-start via linear SUR/OLS
  lin_y <- stats::lm(y ~ X); lin_f <- stats::lm(f ~ X)
  a0y <- unname(coef(lin_y)[1]); b0y <- unname(coef(lin_y)[-1])
  a0f <- unname(coef(lin_f)[1]); b0f <- unname(coef(lin_f)[-1])
  a0y[!is.finite(a0y)]<-0; a0f[!is.finite(a0f)]<-0; b0y[!is.finite(b0y)]<-0; b0f[!is.finite(b0f)]<-0
  
  par_linear <- function(){
    if (spec=="standard_lstar"){
      c(a0y, 0, b0y, rep(0,R),  a0f, 0, b0f, rep(0,R),  0, 0)
    } else {
      c(a0y, b0y, rep(0,R),      a0f, b0f, rep(0,R),    0, 0)
    }
  }
  make_start <- function(g0,c0){
    if (spec=="standard_lstar"){
      c(a0y, 0, b0y, rep(0,R), a0f, 0, b0f, rep(0,R), log(g0), c0)
    } else {
      c(a0y, b0y, rep(0,R),     a0f, b0f, rep(0,R),    log(g0), c0)
    }
  }
  
  ## Objective: joint SSE across equations
  obj <- function(th){ pr <- .star_eval(th, X, spec, z_std, q); sum((y-pr$yhat)^2) + sum((f-pr$fhat)^2) }
  
  st0 <- par_linear(); best <- list(val=obj(st0), par=st0, conv=0L)
  if (trace) message(sprintf("STAR[%s] baseline SSE=%.6f", spec, best$val))
  
  qfun <- function(z,key) switch(key,
                                 median=stats::median(z), q25=stats::quantile(z,0.25), q75=stats::quantile(z,0.75))
  ccands <- unique(vapply(starts_c, function(k) qfun(z_std,k), numeric(1)))
  
  for (g0 in starts_gamma) for (c0 in ccands){
    st <- make_start(g0,c0)
    opt1 <- try(stats::optim(st, obj, method="BFGS", control=list(maxit=maxit,reltol=1e-10)), silent=TRUE)
    cand <- NULL
    if (!inherits(opt1,"try-error") && is.finite(opt1$value)) cand <- opt1
    if (is.null(cand)){
      opt2 <- try(stats::optim(st, obj, method="Nelder-Mead", control=list(maxit=maxit,reltol=1e-8)), silent=TRUE)
      if (!inherits(opt2,"try-error") && is.finite(opt2$value)) cand <- opt2
    }
    if (!is.null(cand) && cand$value < best$val){ best <- list(val=cand$value, par=cand$par, conv=cand$convergence)
    if (trace) message(sprintf("  improved SSE=%.6f (gamma0=%.3g, c0=%.3f)", best$val, g0, c0))
    }
  }
  
  theta <- best$par
  pr <- .star_eval(theta, X, spec, z_std, q)
  res_y <- y - pr$yhat; res_f <- f - pr$fhat
  s2 <- (mean(res_y^2) + mean(res_f^2))/2
  
  V <- matrix(NA_real_, length(theta), length(theta))
  if (isTRUE(compute_vcov)){
    if (!requireNamespace("numDeriv", quietly=TRUE)) stop("Please install numDeriv")
    fstack <- function(th){ p <- .star_eval(th, X, spec, z_std, q); c(p$yhat, p$fhat) }
    G <- try(numDeriv::jacobian(fstack, theta), silent=TRUE)  ## (2n x k)
    if (!inherits(G,"try-error")){
      XtX <- crossprod(G)
      ok  <- try(chol2inv(chol(XtX)), silent=TRUE)
      if (!inherits(ok,"try-error")) V <- ok * s2
    }
  }
  
  structure(list(spec=spec, z_choice=z_choice, q=q, z_mu=z_mu, z_sd=z_sd,
                 coef=theta, vcov=V, s2=s2, n=n, X=X, y=y, f=f, z_std=z_std),
            class="star_bivar_fit")
}


estimate_star_bivar_sep <- function(
    df, spec = c("standard_lstar","minimal_lstar","minimal_estar"),
    z_choice = c("y","f","combo"), q = c(1,0),
    grid_gamma = c(0.5, 1, 2, 4),       # slightly leaner default grid
    grid_c     = c("q25","median","q75"),
    refine     = TRUE, maxit_refine = 150, trace = FALSE,
    errors     = c("iid","AR"),
    ar_order   = 1L,
    vcov_mode  = c("plugin_linear","none")
){
  spec <- match.arg(spec); z_choice <- match.arg(z_choice)
  errors <- match.arg(errors); vcov_mode <- match.arg(vcov_mode)
  stopifnot(all(c("y","f","y_l1","f_l1") %in% names(df)))
  
  n <- nrow(df)
  y <- as.numeric(df$y); f <- as.numeric(df$f)
  X <- as.matrix(df[,c("y_l1","f_l1")]); R <- ncol(X)
  
  if (z_choice=="y") q <- c(1,0) else if (z_choice=="f") q <- c(0,1) else stopifnot(length(q)==R)
  z_raw <- drop(X %*% q)
  z_mu  <- stats::median(z_raw); z_sd <- stats::sd(z_raw); if (!is.finite(z_sd) || z_sd==0) z_sd <- stats::mad(z_raw)
  z_std <- (z_raw - z_mu)/z_sd
  
  G_logistic <- function(z, eta, c) { plogis(exp(eta)*(z - c)) }
  G_estar    <- function(z, eta, c) { 1 - exp(-exp(eta)*(z - c)^2) }
  
  build_design <- function(G){
    if (spec=="standard_lstar") {
      Z <- cbind(1, G, X, X*G); W <- Z
      ny <- 2 + 2*R; nf <- ny
      idx <- list(
        y = list(a0=1, a1=2, b0=3:(2+R), b1=(3+R):(2+2*R)),
        f = list(a0=1, a1=2, b0=3:(2+R), b1=(3+R):(2+2*R))
      )
    } else {
      Z <- cbind(1, X, X*G); W <- Z
      ny <- 1 + 2*R; nf <- ny
      idx <- list(
        y = list(a0=1, b0=2:(1+R), d=(2+R):(1+2*R)),
        f = list(a0=1, b0=2:(1+R), d=(2+R):(1+2*R))
      )
    }
    list(Z=Z, W=W, idx=idx, ny=ny, nf=nf)
  }
  
  # i.i.d. fit (cheap)
  fit_two_eq_iid <- function(G){
    des <- build_design(G)
    ly <- stats::lm.fit(des$Z, y)
    lf <- stats::lm.fit(des$W, f)
    resy <- y - drop(des$Z %*% ly$coefficients)
    resf <- f - drop(des$W %*% lf$coefficients)
    SSE  <- sum(resy^2) + sum(resf^2)
    Vy <- Vf <- NULL
    if (vcov_mode=="plugin_linear"){
      s2y <- sum(resy^2)/max(n - ncol(des$Z),1L)
      s2f <- sum(resf^2)/max(n - ncol(des$W),1L)
      Vy  <- tryCatch(s2y * chol2inv(chol(crossprod(des$Z))), error=function(e) NULL)
      Vf  <- tryCatch(s2f * chol2inv(chol(crossprod(des$W))), error=function(e) NULL)
    }
    list(SSE=SSE, beta_y=as.numeric(ly$coefficients), beta_f=as.numeric(lf$coefficients),
         Z=des$Z, W=des$W, idx=des$idx, Vy=Vy, Vf=Vf)
  }
  
  # AR(P) GLS (used only once at final (eta,c))
  fit_two_eq_AR <- function(G){
    des <- build_design(G)
    gy <- .fit_eq_gls_AR(y, des$Z, ar_order = ar_order)
    gf <- .fit_eq_gls_AR(f, des$W, ar_order = ar_order)
    SSE <- gy$s2 * (nrow(des$Z) - gy$drop - ncol(des$Z)) +
      gf$s2 * (nrow(des$W) - gf$drop - ncol(des$W))
    list(SSE=SSE,
         beta_y=gy$beta, beta_f=gf$beta,
         Z=des$Z, W=des$W, idx=des$idx,
         Vy = if (vcov_mode=="plugin_linear") gy$V else NULL,
         Vf = if (vcov_mode=="plugin_linear") gf$V else NULL,
         phi_y = gy$phi, phi_f = gf$phi,
         drop = max(gy$drop, gf$drop))
  }
  
  qfun <- function(z,key) switch(key, median=stats::median(z), q25=stats::quantile(z,0.25), q75=stats::quantile(z,0.75))
  c_grid <- unique(vapply(grid_c, function(k) qfun(z_std,k), numeric(1)))
  
  # Coarse grid (ALWAYS i.i.d. for speed when errors=="AR") 
  best <- list(SSE = Inf, eta = NA_real_, c0 = NA_real_)
  for (g0 in grid_gamma) for (c0 in c_grid){
    G <- if (spec=="minimal_estar") G_estar(z_std, log(g0), c0) else G_logistic(z_std, log(g0), c0)
    cand <- fit_two_eq_iid(G)
    if (cand$SSE < best$SSE) best <- c(cand, list(eta = log(g0), c0 = c0))
  }
  if (isTRUE(trace)) message(sprintf("STAR[%s,%s] coarse best SSE(iid)=%.4f",
                                     spec, errors, best$SSE))
  
  # Optional refinement in (eta,c) (use i.i.d. SSE; cheap)
  if (isTRUE(refine)){
    sse_of_iid <- function(par){
      G <- if (spec=="minimal_estar") G_estar(z_std, par[1], par[2]) else G_logistic(z_std, par[1], par[2])
      fit_two_eq_iid(G)$SSE
    }
    lower_c <- min(z_std); upper_c <- max(z_std)
    opt <- try(stats::optim(c(best$eta, best$c0), sse_of_iid, method="L-BFGS-B",
                            control=list(maxit=maxit_refine, factr=1e7),
                            lower=c(-5, lower_c), upper=c(5, upper_c)), silent=TRUE)
    if (!inherits(opt,"try-error") && is.finite(opt$value)){
      best$eta <- opt$par[1]; best$c0 <- opt$par[2]
      if (isTRUE(trace)) message(sprintf(" refined SSE(iid)=%.4f", opt$value))
    }
  }
  
  # Final fit at (eta,c): i.i.d. or single GLS-AR, depending on 'errors' 
  G_final <- if (spec=="minimal_estar") G_estar(z_std, best$eta, best$c0) else G_logistic(z_std, best$eta, best$c0)
  final <- if (errors=="iid") fit_two_eq_iid(G_final) else fit_two_eq_AR(G_final)
  
  # pack theta in original order + append (eta,c)
  if (spec=="standard_lstar"){
    th <- c(
      final$beta_y[1], final$beta_y[2], final$beta_y[3:(2+R)], final$beta_y[(3+R):(2+2*R)],
      final$beta_f[1], final$beta_f[2], final$beta_f[3:(2+R)], final$beta_f[(3+R):(2+2*R)],
      best$eta, best$c0
    )
    k_lin_y <- (2+2*R); k_lin_f <- (2+2*R)
  } else {
    th <- c(
      final$beta_y[1], final$beta_y[2:(1+R)], final$beta_y[(2+R):(1+2*R)],
      final$beta_f[1], final$beta_f[2:(1+R)], final$beta_f[(2+R):(1+2*R)],
      best$eta, best$c0
    )
    k_lin_y <- (1+2*R); k_lin_f <- (1+2*R)
  }
  
  # assemble block-diagonal vcov for linear parts (eta,c fixed)
  V <- matrix(0, length(th), length(th))
  if (!is.null(final$Vy)){
    V[1:k_lin_y, 1:k_lin_y] <- final$Vy
    i0 <- k_lin_y
    V[(i0+1):(i0+k_lin_f), (i0+1):(i0+k_lin_f)] <- final$Vf
  } else {
    V[,] <- NA_real_
  }
  
  structure(list(spec=spec, z_choice=z_choice, q=q, z_mu=z_mu, z_sd=z_sd,
                 coef=th, vcov=V,
                 s2 = NA_real_, n = n,
                 X = X, y = y, f = f, z_std = z_std,
                 eta = best$eta, c0 = best$c0,
                 errors = errors, ar_order = ar_order,
                 phi_y = final$phi_y %||% numeric(0),
                 phi_f = final$phi_f %||% numeric(0)),
            class="star_bivar_fit")
}



## Local A1 at evaluation point x_eval = (y_{t-1}, f_{t-1}). Uses standardized z with stored (mu, sd).
star_local_A1 <- function(fit, x_eval){
  stopifnot(inherits(fit, "star_bivar_fit"))
  Xrow <- as.numeric(x_eval)                 # ensure plain numeric vector
  R    <- length(Xrow)
  q    <- as.numeric(fit$q)
  stopifnot(length(q) == R)
  
  # standardized transition at x_eval
  z    <- (sum(q * Xrow) - fit$z_mu) / fit$z_sd
  dzdx <- q / fit$z_sd                      
  
  th   <- fit$coef
  spec <- fit$spec
  idx  <- 0L
  
  if (spec == "standard_lstar") {
    a0y <- th[idx <- idx+1]; a1y <- th[idx <- idx+1]
    b0y <- th[(idx+1):(idx+R)]; idx <- idx+R
    b1y <- th[(idx+1):(idx+R)]; idx <- idx+R
    
    a0f <- th[idx <- idx+1]; a1f <- th[idx <- idx+1]
    b0f <- th[(idx+1):(idx+R)]; idx <- idx+R
    b1f <- th[(idx+1):(idx+R)]; idx <- idx+R
    
    eta <- th[idx <- idx+1]; c0 <- th[idx <- idx+1]
    
    G    <- .star_G_logistic(z, eta, c0)
    dGdz <- .star_dGdz_logistic(z, eta, c0)
    
    # term_y = a1y + b1y' x ; term_f = a1f + b1f' x
    term_y <- a1y + sum(b1y * Xrow)
    term_f <- a1f + sum(b1f * Xrow)
    
    #  = b0y + b1y*G + (a1y + b1y'x) * dGdz * dzdx
    dy_dx <- b0y + b1y * G + term_y * dGdz * dzdx
    #  = b0f + b1f*G + (a1f + b1f'x) * dGdz * dzdx
    df_dx <- b0f + b1f * G + term_f * dGdz * dzdx
    
  } else {
    a0y <- th[idx <- idx+1]
    b0y <- th[(idx+1):(idx+R)]; idx <- idx+R
    dy  <- th[(idx+1):(idx+R)]; idx <- idx+R
    
    a0f <- th[idx <- idx+1]
    b0f <- th[(idx+1):(idx+R)]; idx <- idx+R
    df  <- th[(idx+1):(idx+R)]; idx <- idx+R
    
    eta <- th[idx <- idx+1]; c0 <- th[idx <- idx+1]
    
    if (spec == "minimal_lstar") {
      G    <- .star_G_logistic(z, eta, c0)
      dGdz <- .star_dGdz_logistic(z, eta, c0)
    } else {
      G    <- .star_G_estar(z, eta, c0)
      dGdz <- .star_dGdz_estar(z, eta, c0)
    }
    
    # term_y = dy' x ; term_f = df' x
    term_y <- sum(dy * Xrow)
    term_f <- sum(df * Xrow)
    
    # = b0y + dy*G + (dy'x) * dGdz * dzdx
    dy_dx <- b0y + dy * G + term_y * dGdz * dzdx
    # = b0f + df*G + (df'x) * dGdz * dzdx
    df_dx <- b0f + df * G + term_f * dGdz * dzdx
  }
  
  # Assemble local Jacobian
  A1 <- rbind(dy_dx, df_dx)
  
  # Safe, informative names
  rownames(A1) <- c("y", "f")
  if (R == 2) {
    colnames(A1) <- c("y.l1", "f.l1")
  } else {
    colnames(A1) <- paste0("x", seq_len(R), ".l1")
  }
  A1
}

## Delta-method test for a2*b1<0 using numerical gradient wrt theta 
star_test_cycle_condition <- function(fit, x_eval, alpha=0.10){ 
if (!requireNamespace("numDeriv", quietly=TRUE)) stop("Please install numDeriv")
A1 <- star_local_A1(fit, x_eval) 
est <- as.numeric(A1[1,2] * A1[2,1]) 
if (!all(is.finite(c(est))) || !all(is.finite(fit$vcov))) 
return(list(stat=est,se=NA_real_,t=NA_real_,pval=NA_real_)) 
f_theta <- function(th){ f2 <- fit; f2$coef <- th; A <- 
star_local_A1(f2, x_eval); as.numeric(A[1,2]*A[2,1]) } 
g <- try(numDeriv::grad(f_theta, fit$coef), silent=TRUE) 
if (inherits(g,"try-error") || any(!is.finite(g))) 
return(list(stat=est,se=NA_real_,t=NA_real_,pval=NA_real_)) 
V <- fit$vcov 
se2 <- drop(t(g) %*% V %*% g) 
if (!is.finite(se2) || se2<0) 
return(list(stat=est,se=NA_real_,t=NA_real_,pval=NA_real_)) 
se <- sqrt(se2); if (se==0) 
return(list(stat=est,se=se,t=NA_real_,pval=1.0)) 
tval <- est/se; pval <- pt(tval, df = 2*fit$n - length(fit$coef),
lower.tail=TRUE) 
list(stat=est,se=se,t=tval,pval=pval) }


star_test_cycle_condition_fast <- function(fit, x_eval, alpha=0.10){
  if (!requireNamespace("numDeriv", quietly=TRUE))
    stop("Please install numDeriv")
  
  # indices: all except the last two (eta, c0) are linear
  k <- length(fit$coef)
  idx_lin <- seq_len(k-2L)
  Vlin <- fit$vcov[idx_lin, idx_lin, drop=FALSE]
  
  A1 <- star_local_A1(fit, x_eval)
  est <- as.numeric(A1[1,2]*A1[2,1])
  
  f_lin <- function(th_lin){
    th <- fit$coef; th[idx_lin] <- th_lin
    f2 <- fit; f2$coef <- th
    A <- star_local_A1(f2, x_eval)
    as.numeric(A[1,2]*A[2,1])
  }
  g <- try(numDeriv::grad(f_lin, fit$coef[idx_lin]), silent=TRUE)
  if (inherits(g,"try-error") || any(!is.finite(g)) || any(!is.finite(Vlin))) {
    return(list(stat=est,se=NA_real_,t=NA_real_,pval=NA_real_))
  }
  se2 <- drop(t(g) %*% Vlin %*% g); if (!is.finite(se2) || se2<0) se2 <- NA_real_
  if (!is.finite(se2)) return(list(stat=est,se=NA_real_,t=NA_real_,pval=NA_real_))
  se <- sqrt(se2); tval <- if (se==0) 0 else est/se
  pval <- pt(tval, df = 2*fit$n - length(idx_lin), lower.tail=TRUE)
  list(stat=est,se=se,t=tval,pval=pval)
}






##### 5. Single replication: VAR or STAR
## Build lagged bivariate sample from raw simulated y,f (T x 2 matrix or data.frame)
make_lag_df <- function(yf){
  stopifnot(ncol(yf)==2)
  colnames(yf) <- c("y","f")
  data.frame(y = yf[-1,1], f = yf[-1,2], y_l1 = yf[-nrow(yf),1], f_l1 = yf[-nrow(yf),2])
}

one_replication_VAR <- function(yf, shock_rho, max_lags=4, alpha_ser=0.10, alpha_norm=0.10){
  colnames(yf) <- c("y","f")
  if (isTRUE(all.equal(shock_rho,0))){ p<-1; VARobj <- vars::VAR(yf, p=p, type="const"); pass<-TRUE } else {
    p <- 2
    repeat{
      VARobj <- vars::VAR(yf, p=p, type="const"); pass <- serial_pass(VARobj, lags_pt=p, alpha=alpha_ser)
      if (pass || p==max_lags) break; p<-p+1
    }
    if (!pass && p==max_lags) return(NULL)
  }
  jb_y <- try(tseries::jarque.bera.test(stats::residuals(VARobj$varresult$y)), silent=TRUE)
  jb_f <- try(tseries::jarque.bera.test(stats::residuals(VARobj$varresult$f)), silent=TRUE)
  p_jb_y <- if(!inherits(jb_y,"try-error")) jb_y$p.value else NA_real_
  p_jb_f <- if(!inherits(jb_f,"try-error")) jb_f$p.value else NA_real_
  nonnorm_flag <- as.integer(any(c(p_jb_y,p_jb_f) < alpha_norm, na.rm=TRUE))
  coefs_y <- coef(VARobj$varresult$y); coefs_f <- coef(VARobj$varresult$f)
  A1 <- matrix(c(coefs_y["y.l1"], coefs_y["f.l1"], coefs_f["y.l1"], coefs_f["f.l1"]), 2,2, byrow=TRUE)
  interaction <- A1[1,2]*A1[2,1]
  test_out <- test_cycle_condition(VARobj)
  k <- ncol(yf); p <- VARobj$p
  As <- array(0, c(k,k,p), dimnames=list(colnames(yf),colnames(yf),NULL))
  for (j in 1:p){ As[,,j] <- sapply(colnames(yf), function(eqn) sapply(colnames(yf), function(lv) coef(VARobj$varresult[[eqn]])[paste0(lv, ".l", j)])) }
  top <- do.call(cbind, lapply(1:p, function(j) As[,,j])); bot <- cbind(diag(k*(p-1)), matrix(0,k*(p-1),k)); Cp <- rbind(top,bot)
  eig_vals <- eigen(Cp, only.values=TRUE)$values
  # All eigenvalues (real + complex)
  eig_vals <- eigen(Cp, only.values = TRUE)$values
  # Dominant modulus over ALL eigenvalues
  dom_all      <- eig_vals[ which.max(Mod(eig_vals)) ]
  dom_mod_all  <- Mod(dom_all)
  unstable_all <- as.integer(dom_mod_all > 1 + 1e-10)
  # Keep existing complex-only block for cycle lengths
  is_cplx    <- abs(Im(eig_vals)) > 1e-12
  cplx_vals  <- eig_vals[is_cplx]
  dom_lambda <- if (length(cplx_vals)) cplx_vals[ which.max(Mod(cplx_vals)) ] else NA_complex_
  L_all <- if(length(cplx_vals)) sapply(cplx_vals, cycle_length) else numeric(0)
  L_avg <- if(length(L_all)) mean(L_all[L_all>3 & L_all<20]) else NA_real_
  L_dom <- if(!is.na(dom_lambda)) cycle_length(dom_lambda) else NA_real_
  list(
    estimator = "VAR", lags = p, interaction = interaction, test = test_out,
    dom_lambda = dom_lambda,                # (complex, for L)
    dom_mod_all = dom_mod_all,              # NEW: |lambda|_max (all roots)
    unstable = unstable_all,                # NEW: indicator(|lambda|_max > 1)
    cycle_avg = L_avg, cycle_dom = L_dom,
    jb_p_y = p_jb_y, jb_p_f = p_jb_f, nonnorm = nonnorm_flag
  )
}

one_replication_STAR <- function(yf, 
                                 spec,
                                 z_choice, q,
                                 eval_point, alpha_star,
                                 errors = c("iid","AR"), ar_order = 1L){
  errors <- match.arg(errors)
  df <- make_lag_df(yf)
  
  fit <- estimate_star_bivar_sep(
    df,
    spec       = spec,
    z_choice   = z_choice,
    q          = q,
    errors     = errors,
    ar_order   = ar_order,
    refine     = (errors == "iid"),     # <- only refine for iid (fast)
    vcov_mode  = "plugin_linear",
    trace      = FALSE
  )
  
  x_eval <- switch(eval_point,
                   mean = colMeans(df[,c("y_l1","f_l1")]),
                   zero = c(0,0))
  
  A1  <- star_local_A1(fit, x_eval)
  a2b1 <- A1[1,2] * A1[2,1]
  tst <- star_test_cycle_condition_fast(fit, x_eval, alpha = alpha_star)
  ev  <- eigen(A1, only.values = TRUE)$values
  # Dominant modulus over ALL eigenvalues of A1
  dom_all      <- ev[ which.max(Mod(ev)) ]
  dom_mod_all  <- Mod(dom_all)
  unstable_all <- as.integer(dom_mod_all > 1 + 1e-10)
  # Complex-only for cycles
  evc <- ev[ abs(Im(ev)) > 1e-12 ]
  dom <- if (length(evc)) evc[ which.max(Mod(evc)) ] else NA_complex_
  evc <- ev[abs(Im(ev)) > 1e-12]
  dom <- if (length(evc)) evc[which.max(Mod(evc))] else NA_complex_
  Ls  <- if (length(evc)) sapply(evc, cycle_length) else numeric(0)
  Lavg <- if (length(Ls)) mean(Ls[Ls > 3 & Ls < 20]) else NA_real_
  Ldom <- if (!is.na(dom)) cycle_length(dom) else NA_real_
  
  list(
    estimator = paste0("STAR-", spec), lags = 1L, interaction = a2b1, test = tst,
    dom_lambda = dom,                    # (complex, for L)
    dom_mod_all = dom_mod_all,           # NEW
    unstable = unstable_all,             # NEW
    cycle_avg = Lavg, cycle_dom = Ldom,
    jb_p_y = NA_real_, jb_p_f = NA_real_, nonnorm = NA_integer_
  )
}


one_replication <- function(DGP, T=50, shock_rho=NULL, sd_eps=NULL, eps0=NULL,
                            max_lags=4, alpha_ser=0.10, alpha_star=0.10, alpha_norm=0.10,
                            estimator=c("VAR","STAR-standard_lstar","STAR-minimal_lstar","STAR-minimal_estar"),
                            z_choice=c("y","f","combo"), q=c(1,0), eval_point=c("mean","zero"),
                            ## NEW: STAR error structure
                            star_errors = c("iid","AR"), star_ar_order = 1L){
  estimator  <- match.arg(estimator)
  star_errors <- match.arg(star_errors)
  yf <- simulate_data(DGP, T, shock_rho = shock_rho, sd_eps = sd_eps, eps0 = eps0)
  colnames(yf) <- c("y","f")
  if (startsWith(estimator, "VAR")){
    return(one_replication_VAR(yf, shock_rho = shock_rho, max_lags = max_lags,
                               alpha_ser = alpha_ser, alpha_norm = alpha_norm))
  } else {
    spec <- sub("STAR-","", estimator)
    return(one_replication_STAR(
      yf, spec = spec, z_choice = match.arg(z_choice), q = q,
      eval_point = match.arg(eval_point), alpha_star = alpha_star,
      errors = star_errors, ar_order = star_ar_order
    ))
  }
}



## Helper to print a compact line
print_rep <- function(tag, res) {
  cat(sprintf(
    "%-18s | lags=%s | a2b1=% .4f | L_avg=%6.2f | L_dom=%6.2f | p(a2*b1<0)=% .4f\n",
    tag,
    if (!is.null(res$lags)) res$lags else NA_integer_,
    res$interaction,
    ifelse(is.finite(res$cycle_avg), res$cycle_avg, NA_real_),
    ifelse(is.finite(res$cycle_dom), res$cycle_dom, NA_real_),
    if (!is.null(res$test$pval)) res$test$pval else NA_real_
  ))
}

## Single DGP, choose estimator 
run_example <- function(DGP, estimator = c("VAR","STAR-standard_lstar","STAR-minimal_lstar","STAR-minimal_estar"),
                        T = 50, shock_rho = 0, max_lags = 4,
                        alpha_ser = 0.10, alpha_star = 0.10, alpha_norm = 0.10,
                        z_choice = c("y","f","combo"), q = c(1,0), eval_point = c("mean","zero")) {
  estimator <- match.arg(estimator)
  z_choice  <- match.arg(z_choice)
  eval_point<- match.arg(eval_point)
  res <- one_replication(
    DGP = DGP, T = T, shock_rho = shock_rho, max_lags = max_lags,
    alpha_ser = alpha_ser, alpha_star = alpha_star, alpha_norm = alpha_norm,
    estimator = estimator, z_choice = z_choice, q = q, eval_point = eval_point
  )
  print_rep(paste0(DGP, " [", estimator, "]"), res)
  invisible(res)
}
## Multiple estimators for one DGP 
run_example_all_specs <- function(DGP, estimators = c("VAR","STAR-standard_lstar","STAR-minimal_lstar","STAR-minimal_estar"),
                                  T = 50, shock_rho = 0, ...) {
  out <- lapply(estimators, function(est)
    run_example(DGP, estimator = est, T = T, shock_rho = shock_rho, ...)
  )
  names(out) <- estimators
  invisible(out)
}



############################################################
## EXAMPLE SINGLE RUN.
## EX1: SM39 — pick VAR or any STAR
ex1_VAR   <- run_example("SM39", estimator = "VAR",                  T = 50, shock_rho = 0)
ex1_STD   <- run_example("SM39", estimator = "STAR-standard_lstar",  T = 50, shock_rho = 0, z_choice="y")
ex1_MIN   <- run_example("SM39", estimator = "STAR-minimal_lstar",   T = 50, shock_rho = 0, z_choice="y")
ex1_ESTAR <- run_example("SM39", estimator = "STAR-minimal_estar",   T = 50, shock_rho = 0, z_choice="y")

## EX2: BGP20 — same pattern
ex2_VAR   <- run_example("BGP20", estimator = "VAR",                 T = 50, shock_rho = 0)
ex2_STD   <- run_example("BGP20", estimator = "STAR-standard_lstar", T = 50, shock_rho = 0, z_choice="y")
ex2_MIN   <- run_example("BGP20", estimator = "STAR-minimal_lstar",  T = 50, shock_rho = 0, z_choice="y")
ex2_ESTAR <- run_example("BGP20", estimator = "STAR-minimal_estar",  T = 50, shock_rho = 0, z_choice="y")

## EX3: Puu — same pattern
ex3_VAR   <- run_example("Puu", estimator = "VAR",                   T = 50, shock_rho = 0)
ex3_STD   <- run_example("Puu", estimator = "STAR-standard_lstar",   T = 50, shock_rho = 0, z_choice="y")
ex3_MIN   <- run_example("Puu", estimator = "STAR-minimal_lstar",    T = 50, shock_rho = 0, z_choice="y")
ex3_ESTAR <- run_example("Puu", estimator = "STAR-minimal_estar",    T = 50, shock_rho = 0, z_choice="y")

# HENON IMPORTANT AS WON'T DO IT IN MONTE-CARLO!!! 
## EXFALSE: Henon - same pattern
exFalse_VAR_WN   <- run_example("Henon", estimator = "VAR",                   T = 50, shock_rho = 0)
exFalse_STD_WN    <- run_example("Henon", estimator = "STAR-standard_lstar",   T = 50, shock_rho = 0, z_choice="y")
exFalse_MIN_WN    <- run_example("Henon", estimator = "STAR-minimal_lstar",    T = 50, shock_rho = 0, z_choice="y")
exFalse_ESTAR_WN  <- run_example("Henon", estimator = "STAR-minimal_estar",    T = 50, shock_rho = 0, z_choice="y")

exFalse_VAR_AR   <- run_example("Henon", estimator = "VAR",                   T = 50, shock_rho = 0.8)
exFalse_STD_AR   <- run_example("Henon", estimator = "STAR-standard_lstar",   T = 50, shock_rho = 0.8, z_choice="y")
exFalse_MIN_AR   <- run_example("Henon", estimator = "STAR-minimal_lstar",    T = 50, shock_rho = 0.8, z_choice="y")
exFalse_ESTAR_AR <- run_example("Henon", estimator = "STAR-minimal_estar",    T = 50, shock_rho = 0.8, z_choice="y")







##### 6.  MONTE-CARLO WRAPPER + (STAR AR-errors support)

run_mcarlo <- function(
    DGP, N = 1000, T = 50, shock_rho = NULL, max_lags = 4,
    alpha_ser = 0.10, alpha_norm = 0.10, alpha_star = 0.10,
    estimator = c("VAR","STAR-standard_lstar","STAR-minimal_lstar","STAR-minimal_estar"),
    z_choice = c("y","f","combo"), q = c(1,0), eval_point = c("mean","zero"),
    parallel = TRUE, workers = NULL, seed = 123, progress = TRUE,
    ## NEW:
    star_errors = c("iid","AR"), star_ar_order = 1L
){
  estimator   <- match.arg(estimator)
  z_choice    <- match.arg(z_choice)
  eval_point  <- match.arg(eval_point)
  star_errors <- match.arg(star_errors)
  
  if (!requireNamespace("pbapply", quietly = TRUE)) stop("Please install pbapply")
  if (isTRUE(parallel) && !requireNamespace("parallel", quietly = TRUE)) stop("parallel package is required")
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # Per-replication runner (pure function in terms of args)
  runner <- function(i){
    ok <- try(
      one_replication(
        DGP         = DGP,
        T           = T,
        shock_rho   = shock_rho,
        max_lags    = max_lags,
        alpha_ser   = alpha_ser,
        alpha_norm  = alpha_norm,
        alpha_star  = alpha_star,
        estimator   = estimator,
        z_choice    = z_choice,
        q           = q,
        eval_point  = eval_point,
        ## NEW
        star_errors = star_errors,
        star_ar_order = star_ar_order
      ),
      silent = TRUE
    )
    if (inherits(ok, "try-error") || is.null(ok)) return(NULL)
    
    data.frame(
      estimator   = ok$estimator,
      lags        = ok$lags,
      a2b1        = ok$interaction,
      detect_sign = as.integer(ok$interaction < 0),
      detect_star = as.integer((ok$interaction < 0) && is.finite(ok$test$pval) && (ok$test$pval < 0.10)),
      lambda_re   = if (!is.na(ok$dom_lambda)) Re(ok$dom_lambda) else NA_real_,
      lambda_im   = if (!is.na(ok$dom_lambda)) Im(ok$dom_lambda) else NA_real_,
      lambda_mod  = if (!is.na(ok$dom_lambda)) Mod(ok$dom_lambda) else NA_real_,
      lambda_dom_mod_all = ok$dom_mod_all,     # NEW: |lambda|_max per rep
      unstable    = ok$unstable,               # NEW: 1{|lambda|_max>1}
      cycle_avg   = ok$cycle_avg,
      cycle_dom   = ok$cycle_dom,
      nonnorm     = ok$nonnorm
    )
  }
  
  # Progress bar style
  pbapply::pboptions(type = if (isTRUE(progress)) "timer" else "none")
  
  if (!isTRUE(parallel)) {
    # SEQUENTIAL
    lst <- pbapply::pblapply(seq_len(N), runner)
  } else {
    # PARALLEL (portable PSOCK cluster)
    workers <- workers %||% max(1L, parallel::detectCores() - 1L)
    cl <- parallel::makeCluster(workers, type = "PSOCK")
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    
    # Reproducible parallel RNG
    if (!is.null(seed)) {
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)
      parallel::clusterSetRNGStream(cl, iseed = seed)
    }
    
    # Ensure required packages are attached on workers
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(stats)
        library(pbapply)
        library(vars)
        library(tseries)
      })
      NULL
    })
    
    # Export functions defined in the global environment
    funs <- c(
      # core MC API
      "one_replication","one_replication_VAR","one_replication_STAR","make_lag_df",
      # STAR machinery (used within one_replication_STAR / estimate_star_bivar_sep)
      "estimate_star_bivar","estimate_star_bivar_sep",
      "star_local_A1","star_test_cycle_condition","star_test_cycle_condition_fast",
      ".star_G_logistic",".star_dGdz_logistic",".star_G_estar",".star_dGdz_estar",
      ".star_build_X",".star_eval",
      # NEW AR-GLS helpers
      ".gls_transform_AR",".fit_eq_gls_AR",
      # helpers & tests
      "serial_pass","cycle_length","test_cycle_condition","true_a2b1",
      # simulators
      "simulate_data","simulate_SM39","simulate_Metzler","simulate_BGP20",
      "simulate_Goodwin","simulate_Puu","simulate_Henon"
    )
    funs <- funs[funs %in% ls(envir = .GlobalEnv)]
    if (length(funs)) parallel::clusterExport(cl, varlist = funs, envir = .GlobalEnv)
    
    # Export scalar arguments from *this* function's environment
    scalars <- c("DGP","T","shock_rho","max_lags","alpha_ser","alpha_norm","alpha_star",
                 "estimator","z_choice","q","eval_point","star_errors","star_ar_order")
    parallel::clusterExport(cl, varlist = scalars, envir = environment())
    
    # Parallel apply with progress bar
    lst <- pbapply::pblapply(seq_len(N), runner, cl = cl)
  }
  
  valid <- Filter(Negate(is.null), lst)
  if (!length(valid)) stop("All replications failed.")
  results <- do.call(rbind, valid)
  if (nrow(results) == 0) stop("All replications failed.")
  
  a_true <- true_a2b1(DGP)
  rmse   <- sqrt(mean((results$a2b1 - a_true)^2))
  sdhat  <- stats::sd(results$a2b1)
  nrmse  <- if (is.finite(sdhat) && sdhat > 0) rmse/sdhat else NA_real_
  n_acc  <- nrow(results)
  nonn   <- mean(results$nonnorm, na.rm = TRUE)
  
  dom_mod_mean <- mean(results$lambda_dom_mod_all, na.rm = TRUE)
  dom_mod_sd   <-  sd(results$lambda_dom_mod_all,  na.rm = TRUE)
  instab_rate  <- mean(results$unstable, na.rm = TRUE)   # share(|lambda|_max>1)
  
  
  data.frame(
    DGP, estimator, n = n_acc,
    mean_a2b1 = mean(results$a2b1), true_a2b1 = a_true,
    RMSE_a2b1 = rmse, NRMSE_a2b1 = nrmse,
    detect_sign = mean(results$detect_sign),
    detect_star = mean(results$detect_star, na.rm = TRUE),
    mean_cycle_avg = mean(results$cycle_avg, na.rm = TRUE),
    sd_cycle_avg   = sd(results$cycle_avg,   na.rm = TRUE),
    mean_cycle_dom = mean(results$cycle_dom, na.rm = TRUE),
    sd_cycle_dom   = sd(results$cycle_dom,   na.rm = TRUE),
    nonnorm_pct    = 100 * nonn,
    ## NEW:
    mean_dom_mod   = dom_mod_mean,
    sd_dom_mod     = dom_mod_sd,
    instab_pct     = 100 * instab_rate
  )
  
}





######### 7. RESULTS: SYSTEMATIC MONTE CARLO (VAR + STAR) across all DGPs and estimators
 library(pbapply)
  library(dplyr)
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  pboptions(type = "timer")  # nice progress bars
  
  models     <- c("SM39", "BGP20", "Goodwin", "Puu")  # (keep Hénon out of MC to avoid bugs)
  estimators <- c("VAR", "STAR-standard_lstar", "STAR-minimal_lstar", "STAR-minimal_estar")
  
  ## Tune N per estimator (STAR now much faster with separable fit)
  N_per <- c(
    "VAR"                 = 1000,
    "STAR-standard_lstar" = 1000,
    "STAR-minimal_lstar"  = 1000,
    "STAR-minimal_estar"  = 1000
  )
  
  ## Recommended: parallelize INSIDE run_mcarlo, keep outer loop sequential.
  ## (Avoid nested clusters.)
  workers_mc <- max(1L, parallel::detectCores() - 1L)
  seed_base  <- 123
  
  .run_grid <- function(models, estimators, shock_rho_value, scenario_tag,
                        T = 50, max_lags = 4,
                        alpha_ser = 0.10, alpha_norm = 0.10, alpha_star = 0.10,
                        z_choice = "y", q = c(1,0), eval_point = "mean") {
    
    grid <- expand.grid(DGP = models, estimator = estimators, stringsAsFactors = FALSE)
    
    # Outer loop: one (DGP × estimator) at a time.
    res_list <- pblapply(seq_len(nrow(grid)), function(i) {
      g   <- grid[i, ]
      N_i <- N_per[[g$estimator]] %||% 1000
      
      # Important: give a different base seed per grid job for reproducibility
      seed_i <- seed_base + i
      
      df <- run_mcarlo(
        DGP        = g$DGP,
        N          = N_i,
        T          = T,
        shock_rho  = shock_rho_value,   # 0 for WN; NULL => use DGP default AR(1)
        max_lags   = max_lags,
        alpha_ser  = alpha_ser,
        alpha_norm = alpha_norm,
        alpha_star = alpha_star,
        estimator  = g$estimator,
        z_choice   = z_choice,
        q          = q,
        eval_point = eval_point,
        parallel   = FALSE,              # <—— turn on inner parallel MC
        workers    = workers_mc,        # all but one core
        seed       = seed_i,            # reproducible parallel RNG
        progress   = TRUE
      )
      
      df$scenario <- scenario_tag
      df
    })
    
    do.call(rbind, res_list)
  }
  
  ## 1) White-noise shocks (ideal VAR(1) case)
  res_wn <- .run_grid(
    models, estimators,
    shock_rho_value = 0,
    scenario_tag    = "WN",
    z_choice = "y", q = c(1,0), eval_point = "mean"
  )
  
  ## 2) AR(1) shocks (use each DGP's default persistence)
  res_ar1 <- .run_grid(
    models, estimators,
    shock_rho_value = NULL,
    scenario_tag    = "AR1",
    z_choice = "y", q = c(1,0), eval_point = "mean"
  )
  
  ## 3) Combine & tidy
  res_all <- rbind(res_wn, res_ar1)
  rownames(res_all) <- NULL
  
  summary_table <- res_all %>%
    dplyr::select(
      DGP, estimator, scenario,
      n,
      mean_a2b1, true_a2b1,
      RMSE_a2b1, NRMSE_a2b1,
      detect_sign, detect_star,
      mean_cycle_avg, sd_cycle_avg,
      mean_cycle_dom, sd_cycle_dom,
      nonnorm_pct,
      ## NEW: dominant-modulus & instability summary
      mean_dom_mod, sd_dom_mod, instab_pct
    ) %>%
    dplyr::arrange(DGP, scenario, estimator)
  
  print(summary_table)
  
  # (Optional) Save results
  saveRDS(list(summary = summary_table, raw = res_all), file = "MC_fast_parallel.rds")
  
 
















##########################################################################
## FIGURES
##########################################################################

## FIGURE 1) Deterministic vs WN vs AR(1), all DGPs #######
model_names <- c("SM39", "BGP20", "Goodwin", "Puu", "Henon")

# Per-model AR(1) persistence
rho_for <- function(m) 0.8

# Per-model shock standard deviations (match earlier choices)
sd_eps_for <- function(m) switch(m,
                                 "Henon"   = 0.01,
                                 "Puu"     = 0.1,
                                 "Goodwin" = 0.01,
                                 1  # default for SM39, Metzler, BGP20
)

# Small helper to reshape a simulator output into long panel form
make_panel <- function(df, model, scenario) {
  stopifnot(ncol(df) >= 2)
  names(df)[1:2] <- c("y","f")
  df |>
    dplyr::mutate(t = dplyr::row_number(),
                  model = model,
                  scenario = factor(scenario,
                                    levels = c("Deterministic","Deterministic + WN shocks","Deterministic + AR(1) shocks"))) |>
    tidyr::pivot_longer(c("y","f"), names_to = "series", values_to = "value")
}

# Build the three scenarios
det_panels <- lapply(model_names, function(m) {
  df <- do.call(paste0("simulate_", m), list(T = 100, shock_rho = 0,   sd_eps = 0))
  make_panel(df, m, "Deterministic")
})

wn_panels <- lapply(model_names, function(m) {
  df <- do.call(paste0("simulate_", m), list(T = 100, shock_rho = 0,   sd_eps = sd_eps_for(m)))
  make_panel(df, m, "Deterministic + WN shocks")
})

ar1_panels <- lapply(model_names, function(m) {
  df <- do.call(paste0("simulate_", m), list(T = 100, shock_rho = rho_for(m), sd_eps = sd_eps_for(m)))
  make_panel(df, m, "Deterministic + AR(1) shocks")
})

df3 <- dplyr::bind_rows(det_panels, wn_panels, ar1_panels)


use_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)

p3 <- ggplot2::ggplot(df3, ggplot2::aes(t, value, color = series)) +
  ggplot2::geom_line() +
  (
    if (use_ggh4x) {
      ggh4x::facet_grid2(model ~ scenario, scales = "free", independent = "all")
    } else {
      ggplot2::facet_grid(model ~ scenario, scales = "free")
    }
  ) +
  ggplot2::labs(
    title = "Figure 3: Deterministic vs White Noise vs AR(1) (all DGPs)",
    x = "t", y = "", color = "Series"
  ) +
  ggplot2::theme_minimal()

print(p3)






## FIGURE 2) Histograms, BGP20 & SM39, VAR(1) ans STAR with AR(0) shocks #######
# N=1000, T=50; vertical lines = analytical truth.

## Conjugate with non-negative imag part (ties by modulus)
pick_pos_conj <- function(evals) {
  cplx <- evals[abs(Im(evals)) > 1e-12]
  if (!length(cplx)) return(NA_complex_)
  pos <- cplx[Im(cplx) > 0]
  if (!length(pos)) pos <- cplx
  pos[which.max(Mod(pos))]
}

## BGP20 analytical truth in (y=H(I), X) coordinates 
J_BGP20_phi_true <- function(alpha0=-1, alpha1=0.4, alpha2=0.9, alpha3=2, delta=0.2){
  Istar <- 0
  num <- (1 + exp(-Istar))^2
  den <- num - alpha3 * exp(-Istar)
  phi <- num / den  # = 2 at baseline
  matrix(c((alpha2 - alpha1)*phi, -alpha1*(1 - delta)*phi,
           1,                     1 - delta), 2, 2, byrow = TRUE)
}

## SM39 analytical truth 
J_SM39_true <- function(c = 0.4, beta = 2){
  matrix(c(c, c, beta*(c-1), beta*c), 2, 2, byrow = TRUE)
}

## Compact 4-stat extractor from a 2x2 A1
extract_stats_from_A1 <- function(A1){
  ev  <- eigen(A1, only.values = TRUE)$values
  lam <- pick_pos_conj(ev)
  c(
    a2b1 = A1[1,2] * A1[2,1],
    re   = Re(lam),
    im   = Im(lam),
    mod  = Mod(lam)
  )
}

## Generic 4-panel histogram plotter with vertical truth lines
plot_hist_four <- function(df_long, truth_df, title_text, subtitle_text){
  ggplot2::ggplot(df_long, ggplot2::aes(x = val)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey70", color = "black") +
    ggplot2::geom_vline(data = truth_df, ggplot2::aes(xintercept = true),
                        linetype = 2, linewidth = 0.9) +
    ggplot2::facet_wrap(~ var, scales = "free", ncol = 2) +
    ggplot2::labs(title = title_text, subtitle = subtitle_text,
                  x = NULL, y = "Frequency") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}



############# SM39 - VAR(1) histograms ############# 
compute_hist_SM39_VAR <- function(N = 1000, Tn = 50, sd_eps = 1){
  ## Truth
  A_true <- J_SM39_true(c = 0.4, beta = 2)
  lamT   <- pick_pos_conj(eigen(A_true, only.values = TRUE)$values)
  truth_df <- data.frame(
    var  = c("a2b1","re","im","mod"),
    true = c(A_true[1,2]*A_true[2,1], Re(lamT), Im(lamT), Mod(lamT))
  )
  
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  sims <- matrix(NA_real_, nrow = N, ncol = 4,
                 dimnames = list(NULL, c("a2b1","re","im","mod")))
  for (i in seq_len(N)){
    dat <- simulate_SM39(T = Tn, shock_rho = 0, sd_eps = sd_eps)
    names(dat)[1:2] <- c("y","f")
    v1 <- vars::VAR(dat, p = 1, type = "const")
    A1 <- matrix(c(
      coef(v1$varresult$y)["y.l1"], coef(v1$varresult$y)["f.l1"],
      coef(v1$varresult$f)["y.l1"], coef(v1$varresult$f)["f.l1"]
    ), 2, 2, byrow = TRUE)
    sims[i, ] <- extract_stats_from_A1(A1)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  df_long <- tidyr::pivot_longer(as.data.frame(sims), everything(),
                                 names_to = "var", values_to = "val")
  p <- plot_hist_four(df_long, truth_df,
                      "SM39 — VAR(1) with AR(0) shocks",
                      "N = 1,000, T = 50. Dashed = analytical truth.")
  list(data = df_long, truth = truth_df, plot = p)
}

## RUN 
res_SM39_VAR <- compute_hist_SM39_VAR(N = 1000, Tn = 50, sd_eps = 1)
res_SM39_VAR$plot



############# SM39 - STAR histograms ############# 
compute_hist_SM39_STAR <- function(spec = c("standard_lstar","minimal_lstar","minimal_estar"),
                                   z_choice = c("y","f","combo"),
                                   q = c(1,0), eval_point = c("mean","zero"),
                                   N = 1000, Tn = 50, sd_eps = 1){
  spec       <- match.arg(spec)
  z_choice   <- match.arg(z_choice)
  eval_point <- match.arg(eval_point)
  
  ## Truth
  A_true <- J_SM39_true(c = 0.4, beta = 2)
  lamT   <- pick_pos_conj(eigen(A_true, only.values = TRUE)$values)
  truth_df <- data.frame(
    var  = c("a2b1","re","im","mod"),
    true = c(A_true[1,2]*A_true[2,1], Re(lamT), Im(lamT), Mod(lamT))
  )
  
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  sims <- matrix(NA_real_, nrow = N, ncol = 4,
                 dimnames = list(NULL, c("a2b1","re","im","mod")))
  keep <- 0L
  
  for (i in seq_len(N)){
    dat <- simulate_SM39(T = Tn, shock_rho = 0, sd_eps = sd_eps)
    yf  <- setNames(dat[,1:2], c("y","f"))
    df  <- make_lag_df(yf)
    
    fit <- estimate_star_bivar_sep(df, spec = spec, z_choice = z_choice, q = q,
                                    grid_gamma = c(0.5, 1, 2, 4), grid_c = c("q25","median","q75"),
                                    refine = FALSE, trace = FALSE, vcov_mode = "none")
    x_eval <- if (eval_point == "mean") colMeans(df[, c("y_l1","f_l1")]) else c(0,0)
    
    
    A1 <- star_local_A1(fit, x_eval)
    keep <- keep + 1L
    sims[keep, ] <- extract_stats_from_A1(A1)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  sims <- sims[seq_len(keep), , drop = FALSE]
  df_long <- tidyr::pivot_longer(as.data.frame(sims), everything(),
                                 names_to = "var", values_to = "val")
  
  subtitle <- sprintf("STAR-%s (z=%s, q=[%s], eval=%s); N=%d, T=%d. Dashed = analytical truth.",
                      spec, z_choice, paste(q, collapse=","), eval_point, keep, Tn)
  p <- plot_hist_four(df_long, truth_df,
                      "SM39 — STAR with AR(0) shocks",
                      subtitle)
  list(data = df_long, truth = truth_df, plot = p)
}

## RUN ONE STAR SPEC AT A TIME
res_SM39_STD  <- compute_hist_SM39_STAR(spec="standard_lstar", z_choice="y", q=c(1,0), eval_point="mean",
                                        N=1000, Tn=50, sd_eps=1)
res_SM39_STD$plot

res_SM39_MINIMAL  <- compute_hist_SM39_STAR(spec="minimal_lstar", z_choice="y", q=c(1,0), eval_point="mean",
                                        N=1000, Tn=50, sd_eps=1)
res_SM39_MINIMAL$plot

res_SM39_ESTAR  <- compute_hist_SM39_STAR(spec="minimal_estar", z_choice="y", q=c(1,0), eval_point="mean",
                                        N=1000, Tn=50, sd_eps=1)
res_SM39_ESTAR$plot




############# BGP20 - VAR(1) histograms ############# 
compute_hist_BGP20_VAR <- function(N = 1000, Tn = 50, sd_eps = 1){
  ## Truth
  A_true  <- J_BGP20_phi_true(alpha0=-1, alpha1=0.4, alpha2=0.9, alpha3=2, delta=0.2)
  lamT    <- pick_pos_conj(eigen(A_true, only.values = TRUE)$values)
  truth_df <- data.frame(
    var  = c("a2b1","re","im","mod"),
    true = c(A_true[1,2]*A_true[2,1], Re(lamT), Im(lamT), Mod(lamT))
  )
  
  ## MC with progress bar (AR(0) errors as in Appendix E)
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  sims <- matrix(NA_real_, nrow = N, ncol = 4,
                 dimnames = list(NULL, c("a2b1","re","im","mod")))
  for (i in seq_len(N)){
    dat <- simulate_BGP20(T = Tn,
                          alpha0=-1, alpha1=0.4, alpha2=0.9, alpha3=2,
                          delta=0.2, shock_rho = 0, sd_eps = sd_eps)
    ## Transform to y = H(I) so estimates are in the same coordinates as truth
    yH <- dat$I - 2/(1 + exp(-dat$I)) - 1  # alpha3=2, alpha0=-1
    df <- data.frame(y = yH, f = dat$X)
    
    v1 <- vars::VAR(df, p = 1, type = "const")
    A1 <- matrix(c(
      coef(v1$varresult$y)["y.l1"], coef(v1$varresult$y)["f.l1"],
      coef(v1$varresult$f)["y.l1"], coef(v1$varresult$f)["f.l1"]
    ), 2, 2, byrow = TRUE)
    
    sims[i, ] <- extract_stats_from_A1(A1)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  df_long <- tidyr::pivot_longer(as.data.frame(sims), everything(),
                                 names_to = "var", values_to = "val")
  p <- plot_hist_four(df_long, truth_df,
                      "Figure A2: Histograms, BGP20 — VAR(1) with AR(0) shocks",
                      "N = 1,000, T = 50. Dashed = analytical truth (structural A1).")
  list(data = df_long, truth = truth_df, plot = p)
}

## RUN (plots immediately) 
res_BGP20_VAR <- compute_hist_BGP20_VAR(N = 1000, Tn = 50, sd_eps = 1)
res_BGP20_VAR$plot




############# BGP20 - STAR histograms ############# 
compute_hist_BGP20_STAR <- function(spec = c("standard_lstar","minimal_lstar","minimal_estar"),
                                    z_choice = c("y","f","combo"),
                                    q = c(1,0), eval_point = c("mean","zero"),
                                    N = 1000, Tn = 50, sd_eps = 1){
  spec       <- match.arg(spec)
  z_choice   <- match.arg(z_choice)
  eval_point <- match.arg(eval_point)

  ## Truth (same as VAR)
  A_true  <- J_BGP20_phi_true(alpha0=-1, alpha1=0.4, alpha2=0.9, alpha3=2, delta=0.2)
  lamT    <- pick_pos_conj(eigen(A_true, only.values = TRUE)$values)
  truth_df <- data.frame(
    var  = c("a2b1","re","im","mod"),
    true = c(A_true[1,2]*A_true[2,1], Re(lamT), Im(lamT), Mod(lamT))
  )

  pb <- txtProgressBar(min = 0, max = N, style = 3)
  sims <- matrix(NA_real_, nrow = N, ncol = 4,
                 dimnames = list(NULL, c("a2b1","re","im","mod")))
  keep <- 0L

  for (i in seq_len(N)){
    ## Simulate with AR(0) errors (appendix setting)
    dat <- simulate_BGP20(T = Tn,
                          alpha0=-1, alpha1=0.4, alpha2=0.9, alpha3=2,
                          delta=0.2, shock_rho = 0, sd_eps = sd_eps)
    ## Transform to y = H(I)
    yH <- dat$I - 2/(1 + exp(-dat$I)) - 1
    yf <- data.frame(y = yH, f = dat$X)

    ## Estimate STAR on lagged df (local A1 at chosen eval point)
    df  <- make_lag_df(yf)
    fit <- estimate_star_bivar_sep(df, spec = spec, z_choice = z_choice, q = q,
                                   grid_gamma = c(0.5, 1, 2, 4), grid_c = c("q25","median","q75"),
                                   refine = FALSE, trace = FALSE, vcov_mode = "none")
    x_eval <- if (eval_point == "mean") colMeans(df[, c("y_l1","f_l1")]) else c(0,0)
    
    A1 <- star_local_A1(fit, x_eval)
    keep <- keep + 1L
    sims[keep, ] <- extract_stats_from_A1(A1)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  sims <- sims[seq_len(keep), , drop = FALSE]
  df_long <- tidyr::pivot_longer(as.data.frame(sims), everything(),
                                 names_to = "var", values_to = "val")

  subtitle <- sprintf("STAR-%s (z=%s, q=[%s], eval=%s); N=%d, T=%d. Dashed = analytical truth.",
                      spec, z_choice, paste(q, collapse=","), eval_point, keep, Tn)
  p <- plot_hist_four(df_long, truth_df,
                      "Figure A2*: Histograms, BGP20 — STAR with AR(0) shocks",
                      subtitle)
  list(data = df_long, truth = truth_df, plot = p)
}

## RUN ONE STAR SPEC AT A TIME (plots immediately)
res_BGP20_STD   <- compute_hist_BGP20_STAR(spec="standard_lstar", z_choice="y", q=c(1,0), eval_point="mean",
                                           N=1000, Tn=50, sd_eps=1)
res_BGP20_STD$plot

res_BGP20_MINL  <- compute_hist_BGP20_STAR(spec="minimal_lstar",  z_choice="y", q=c(1,0), eval_point="mean",
                                           N=1000, Tn=50, sd_eps=1)
res_BGP20_MINL$plot

res_BGP20_ESTAR  <- compute_hist_BGP20_STAR(spec="minimal_estar",  z_choice="y", q=c(1,0), eval_point="mean",
                                           N=1000, Tn=50, sd_eps=1)
res_BGP20_ESTAR$plot





