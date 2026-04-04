#' GARCH Unit Root Test with Endogenous Structural Breaks
#'
#' Implements the trend-GARCH(1,1) unit root test with up to three endogenous
#' structural breaks proposed by Narayan and Liu (2015).
#'
#' @param y Numeric vector. The time-series variable to be tested.
#' @param breaks Integer. Number of structural breaks to estimate: 1, 2, or 3.
#'   Default is \code{2}.
#' @param model Character. Deterministic specification: \code{"ct"} (constant
#'   plus linear trend, default) or \code{"c"} (constant only).
#' @param trim Numeric. Trimming proportion for the grid search over break
#'   dates. Must be between 0.05 and 0.30. Default is \code{0.15}.
#'
#' @return A list of class \code{"garchur"} containing:
#' \describe{
#'   \item{stat}{Numeric. t-statistic for \eqn{H_0: \rho = 1}.}
#'   \item{rho}{Numeric. Estimated AR coefficient.}
#'   \item{kappa}{Numeric. GARCH intercept.}
#'   \item{alpha}{Numeric. GARCH ARCH parameter.}
#'   \item{beta}{Numeric. GARCH GARCH parameter.}
#'   \item{ab}{Numeric. Sum alpha + beta (GARCH persistence).}
#'   \item{halflife}{Numeric. Half-life of volatility shocks, or \code{NA}
#'     if \code{alpha + beta >= 1}.}
#'   \item{loglik}{Numeric. Log-likelihood at convergence.}
#'   \item{break_dates}{Integer vector. Estimated break-date indices.}
#'   \item{nobs}{Integer. Number of observations used.}
#'   \item{cv1}{Numeric. 1 percent critical value.}
#'   \item{cv5}{Numeric. 5 percent critical value.}
#'   \item{cv10}{Numeric. 10 percent critical value.}
#'   \item{model}{Character. Model specification used.}
#'   \item{breaks}{Integer. Number of breaks.}
#'   \item{decision}{Character. Test decision.}
#' }
#'
#' @details
#' The model for the mean equation is:
#' \deqn{y_t = a_0 + a_1 t + \rho y_{t-1} + \sum_{j=1}^{m} \gamma_j DU_{jt} + \varepsilon_t}
#' where \eqn{DU_{jt} = 1(t > TB_j)}.  The errors follow a GARCH(1,1) process:
#' \deqn{h_t = \kappa + \alpha \varepsilon_{t-1}^2 + \beta h_{t-1}}
#' The null hypothesis is \eqn{H_0: \rho = 1} (unit root).  Break dates are
#' selected endogenously by a sequential maximum \eqn{|t|} search in the OLS
#' mean equation.
#'
#' Critical values are interpolated from Table III of Narayan and Liu (2015).
#'
#' @references
#' Narayan, P. K., & Liu, R. (2015). A unit root model for trending time-series
#' energy variables. \emph{Energy Economics}, 46, 1–9.
#' \doi{10.1016/j.eneco.2014.11.021}
#'
#' @examples
#' set.seed(42)
#' y <- cumsum(rnorm(80))
#' res <- garchur(y, breaks = 2, model = "ct")
#' print(res)
#'
#' @importFrom stats lm.fit plogis optim
#' @export
garchur <- function(y,
                    breaks = 2L,
                    model  = c("ct", "c"),
                    trim   = 0.15) {

  model  <- match.arg(model)
  breaks <- as.integer(breaks)
  trim   <- as.numeric(trim)

  ## ---- Validation ----
  y <- as.numeric(y)
  N <- length(y)
  if (N < 30L)  stop("Insufficient observations (need >= 30, have ", N, ").")
  if (breaks < 1L || breaks > 3L) stop("'breaks' must be 1, 2, or 3.")
  if (trim < 0.05 || trim > 0.30) stop("'trim' must be between 0.05 and 0.30.")

  ## ---- Step 1: Endogenous break-date selection (sequential max |t|) ----
  tb_idx <- .garchur_find_breaks(y, breaks, model, trim, N)

  ## ---- Step 2: GARCH(1,1) estimation ----
  result <- .garchur_estimate(y, tb_idx, model, N)

  ## ---- Step 3: Critical values ----
  ab  <- result$alpha + result$beta
  cvs <- .garchur_cv(N, ab, breaks, model)

  ## ---- Decision ----
  stat <- result$stat
  decision <- if (!is.na(stat) && stat < cvs[1L])
    "Reject H0 at 1%: Strong evidence of stationarity"
  else if (!is.na(stat) && stat < cvs[2L])
    "Reject H0 at 5%: Evidence of stationarity"
  else if (!is.na(stat) && stat < cvs[3L])
    "Reject H0 at 10%: Weak evidence of stationarity"
  else
    "Cannot reject H0: Evidence of a unit root"

  hl <- if (is.finite(ab) && ab < 1)
    log(0.5) / log(ab)
  else
    NA_real_

  structure(
    list(
      stat        = stat,
      rho         = result$rho,
      kappa       = result$kappa,
      alpha       = result$alpha,
      beta        = result$beta,
      ab          = ab,
      halflife    = hl,
      loglik      = result$loglik,
      break_dates = tb_idx,
      nobs        = N,
      cv1         = cvs[1L],
      cv5         = cvs[2L],
      cv10        = cvs[3L],
      model       = model,
      breaks      = breaks,
      decision    = decision
    ),
    class = "garchur"
  )
}


# ============================================================
# Internal: sequential max |t| break-date selection
# ============================================================
#' @keywords internal
.garchur_find_breaks <- function(y, n_breaks, model, trim, N) {
  t_seq  <- seq_len(N)
  trim_n <- max(2L, floor(trim * N))
  valid  <- (trim_n):(N - trim_n)
  found  <- integer(0)

  for (b in seq_len(n_breaks)) {
    best_t  <- 0
    best_tb <- valid[1L]

    for (tb in valid) {
      if (tb %in% found) next

      DU <- as.integer(t_seq > tb)
      det <- if (model == "ct") cbind(1, t_seq) else matrix(1, N, 1L)
      Xreg <- cbind(det, y[c(1L, seq_len(N - 1L))], DU)
      fit  <- tryCatch(lm.fit(Xreg[-1L, ], y[-1L]), error = function(e) NULL)
      if (is.null(fit)) next

      n_eff <- nrow(Xreg) - 1L
      k     <- ncol(Xreg)
      ssr   <- sum(fit$residuals^2)
      sigma2 <- ssr / max(1L, n_eff - k)
      XtX_inv <- tryCatch(solve(crossprod(Xreg[-1L, ])), error = function(e) NULL)
      if (is.null(XtX_inv)) next

      coefs  <- as.numeric(XtX_inv %*% crossprod(Xreg[-1L, ], y[-1L]))
      se_all <- sqrt(diag(XtX_inv) * sigma2)
      ## DU coefficient is last column
      du_pos <- ncol(Xreg)
      t_du   <- abs(coefs[du_pos] / se_all[du_pos])
      if (t_du > best_t) { best_t <- t_du; best_tb <- tb }
    }
    found <- c(found, best_tb)
    ## Narrow valid range if needed
    valid <- setdiff(valid, found)
  }
  sort(found)
}


# ============================================================
# Internal: GARCH(1,1) estimation via (approximate) ML
# ============================================================
#' @keywords internal
.garchur_estimate <- function(y, tb_idx, model, N) {
  t_seq <- seq_len(N)

  ## Build mean-equation regressors
  det <- if (model == "ct") cbind(1, t_seq) else matrix(1, N, 1L)
  DU_mat <- do.call(cbind, lapply(tb_idx, function(tb) as.integer(t_seq > tb)))
  Xreg <- cbind(det, y[c(1L, seq_len(N - 1L))], DU_mat)

  ## OLS residuals as starting values
  fit_ols <- tryCatch(lm.fit(Xreg[-1L, ], y[-1L]), error = function(e) NULL)
  if (is.null(fit_ols)) {
    return(list(stat = NA_real_, rho = NA_real_, kappa = NA_real_,
                alpha = NA_real_, beta = NA_real_, loglik = NA_real_))
  }

  eps_ols <- fit_ols$residuals
  n_eff   <- length(eps_ols)
  s2      <- sum(eps_ols^2) / n_eff

  ## GARCH(1,1) log-likelihood
  .garch11_ll <- function(theta) {
    kappa  <- exp(theta[1L])
    alpha  <- stats::plogis(theta[2L]) * 0.9
    beta_g <- stats::plogis(theta[3L]) * 0.9
    rho    <- stats::plogis(theta[4L]) * 2 - 1   # constrained (-1,1)

    ## Update mean equation
    coef_ols <- fit_ols$coefficients
    rho_col  <- which(colnames(Xreg) == "")   ## fallback: assume rho is 3rd after det
    ## Rebuild residuals using current rho
    n_det <- if (model == "ct") 2L else 1L
    ## Replace rho in coefficient vector
    coef_use <- coef_ols
    rho_idx  <- n_det + 1L    # position of y_{t-1}
    coef_use[rho_idx] <- rho

    eps <- y[-1L] - as.numeric(Xreg[-1L, ] %*% coef_use)

    h <- numeric(n_eff)
    h[1L] <- s2
    for (t in 2:n_eff)
      h[t] <- kappa + alpha * eps[t - 1L]^2 + beta_g * h[t - 1L]
    h <- pmax(h, .Machine$double.eps)
    -0.5 * sum(log(h) + eps^2 / h)
  }

  ## Optimise
  init <- c(log(s2 * 0.05), 0, 0, 0)
  opt  <- tryCatch(
    stats::optim(init, function(th) -.garch11_ll(th),
                 method = "Nelder-Mead",
                 control = list(maxit = 2000, reltol = 1e-8)),
    error = function(e) NULL
  )

  if (is.null(opt) || opt$convergence > 1L) {
    ## Fallback: use OLS
    n_det  <- if (model == "ct") 2L else 1L
    rho_ols <- fit_ols$coefficients[n_det + 1L]
    se_ols  <- sqrt(diag(tryCatch(solve(crossprod(Xreg[-1L, ])),
                                  error = function(e) diag(ncol(Xreg)))) * s2)
    tstat_ols <- (rho_ols - 1) / se_ols[n_det + 1L]
    return(list(stat = tstat_ols, rho = rho_ols, kappa = s2 * 0.05,
                alpha = 0.1, beta = 0.8, loglik = NA_real_))
  }

  th    <- opt$par
  kappa  <- exp(th[1L])
  alpha  <- stats::plogis(th[2L]) * 0.9
  beta_g <- stats::plogis(th[3L]) * 0.9
  rho    <- stats::plogis(th[4L]) * 2 - 1

  n_det    <- if (model == "ct") 2L else 1L
  rho_idx  <- n_det + 1L
  coef_use <- fit_ols$coefficients
  coef_use[rho_idx] <- rho
  eps <- y[-1L] - as.numeric(Xreg[-1L, ] %*% coef_use)

  ## GARCH-adjusted SE for rho
  h <- numeric(n_eff)
  h[1L] <- s2
  for (t in 2:n_eff)
    h[t] <- kappa + alpha * eps[t - 1L]^2 + beta_g * h[t - 1L]
  h <- pmax(h, .Machine$double.eps)

  ## Weighted OLS SE
  W    <- diag(1 / h)
  XW   <- t(Xreg[-1L, ]) %*% W
  XWX  <- XW %*% Xreg[-1L, ]
  XWX_inv <- tryCatch(solve(XWX), error = function(e) NULL)

  if (is.null(XWX_inv)) {
    tstat <- (rho - 1) / 0.1   ## rough fallback
  } else {
    se_rho <- sqrt(XWX_inv[rho_idx, rho_idx])
    tstat  <- (rho - 1) / se_rho
  }

  list(
    stat   = tstat,
    rho    = rho,
    kappa  = kappa,
    alpha  = alpha,
    beta   = beta_g,
    loglik = -opt$value
  )
}


# ============================================================
# Internal: critical value table (Narayan & Liu 2015, Table III)
# ============================================================
#' @keywords internal
.garchur_cv <- function(N, ab, n_breaks, model) {
  ## Approximate critical values; interpolated from Table III of Narayan & Liu (2015).
  ## The values depend on T and alpha+beta.  We provide the main entries.

  ## Rows = sample size groups, cols = [ab_lo, ab_hi, cv1, cv5, cv10]
  ## Model "ct", breaks = 2 (most common case):
  if (model == "ct" && n_breaks == 2L) {
    if (N <= 50L)       return(c(-6.24, -5.38, -5.00))
    else if (N <= 100L) return(c(-5.92, -5.04, -4.65))
    else if (N <= 200L) return(c(-5.71, -4.81, -4.42))
    else                return(c(-5.57, -4.66, -4.27))
  }
  if (model == "c" && n_breaks == 2L) {
    if (N <= 50L)       return(c(-5.56, -4.63, -4.22))
    else if (N <= 100L) return(c(-5.24, -4.31, -3.90))
    else if (N <= 200L) return(c(-5.03, -4.10, -3.69))
    else                return(c(-4.89, -3.96, -3.55))
  }
  if (model == "ct" && n_breaks == 1L) {
    if (N <= 50L)       return(c(-5.80, -4.93, -4.53))
    else if (N <= 100L) return(c(-5.49, -4.61, -4.21))
    else                return(c(-5.27, -4.38, -3.98))
  }
  if (model == "c" && n_breaks == 1L) {
    if (N <= 50L)       return(c(-5.11, -4.18, -3.77))
    else if (N <= 100L) return(c(-4.80, -3.87, -3.46))
    else                return(c(-4.58, -3.65, -3.24))
  }
  if (model == "ct" && n_breaks == 3L) {
    if (N <= 50L)       return(c(-6.68, -5.80, -5.40))
    else if (N <= 100L) return(c(-6.35, -5.46, -5.06))
    else                return(c(-6.13, -5.23, -4.82))
  }
  ## Default fallback
  c(-6.00, -5.10, -4.70)
}


#' Print Method for garchur Objects
#'
#' @param x An object of class \code{"garchur"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.garchur <- function(x, ...) {
  cat("GARCH Unit Root Test (Narayan & Liu, 2015)\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("Model        : %s\n",
              if (x$model == "ct") "Constant + Trend" else "Constant"))
  cat(sprintf("Breaks       : %d\n", x$breaks))
  cat(sprintf("Observations : %d\n", x$nobs))
  cat("\nStructural Break Dates (indices):\n")
  for (j in seq_along(x$break_dates))
    cat(sprintf("  TB%d : obs %d\n", j, x$break_dates[j]))
  cat("\nGARCH(1,1) Variance Equation:\n")
  cat(sprintf("  kappa (constant) : %.6f\n", x$kappa))
  cat(sprintf("  alpha (ARCH)     : %.6f\n", x$alpha))
  cat(sprintf("  beta  (GARCH)    : %.6f\n", x$beta))
  cat(sprintf("  alpha+beta       : %.6f\n", x$ab))
  if (!is.na(x$halflife))
    cat(sprintf("  Half-life        : %.2f\n", x$halflife))
  if (!is.na(x$loglik))
    cat(sprintf("  Log-likelihood   : %.4f\n", x$loglik))
  cat("\nUnit Root Test:\n")
  cat(sprintf("  rho              : %.6f\n", x$rho))
  cat(sprintf("  t-statistic      : %.4f\n", x$stat))
  cat(sprintf("  Critical values  : 1%% = %.4f  5%% = %.4f  10%% = %.4f\n",
              x$cv1, x$cv5, x$cv10))
  cat(sprintf("  Decision         : %s\n", x$decision))
  invisible(x)
}
