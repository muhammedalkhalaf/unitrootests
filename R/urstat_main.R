#' Comprehensive Unified Unit Root and Stationarity Testing
#'
#' Runs a comprehensive battery of unit root and stationarity tests on one or
#' more time series, producing formatted summary tables and an integration-order
#' decision following the Elder and Kennedy (2001) strategy.
#'
#' @param x A numeric vector, \code{ts} object, or a named list/data frame of
#'   numeric vectors.
#' @param tests Character vector of tests to run.  Possible values:
#'   \code{"ADF"}, \code{"PP"}, \code{"KPSS"}, \code{"ERS"}, \code{"ZA"},
#'   \code{"KM"}.  Default \code{"ALL"} runs all supported tests.
#' @param max_lag Integer. Maximum lag order for ADF lag selection (default 12).
#' @param crit Character. Information criterion for ADF lag selection:
#'   \code{"BIC"} (default) or \code{"AIC"}.
#' @param pp_lag Integer. Bandwidth (Newey-West lags) for PP test (default 4).
#' @param kpss_lags Integer. Lag truncation for KPSS (default 8).
#' @param za_trim Numeric. Trimming proportion for Zivot-Andrews test
#'   (default 0.15).
#' @param level Numeric. Significance level used for decisions (default 0.05).
#' @param stars Logical. Print significance stars (default \code{TRUE}).
#' @param strategy Logical. Print Elder-Kennedy decision table
#'   (default \code{TRUE}).
#'
#' @return A list (invisibly) containing one element per series.  Each element
#'   is itself a named list with components:
#'   \describe{
#'     \item{\code{adf}}{Results from ADF tests (list).}
#'     \item{\code{pp}}{Results from PP tests (list).}
#'     \item{\code{kpss}}{Results from KPSS tests (list).}
#'     \item{\code{ers}}{Results from ERS/DF-GLS tests (list).}
#'     \item{\code{za}}{Results from Zivot-Andrews test (list).}
#'     \item{\code{decision}}{Character. Inferred integration order.}
#'     \item{\code{process}}{Character. Suggested data transformation.}
#'   }
#'
#' @references
#' Elder, J. and Kennedy, P. E. (2001). Testing for unit roots: What should
#' students be taught? \emph{Journal of Economic Education}, 32(2), 137-146.
#' \doi{10.1080/00220480109595179}
#'
#' Dickey, D. A. and Fuller, W. A. (1979). Distribution of the estimators for
#' autoregressive time series with a unit root.
#' \emph{Journal of the American Statistical Association}, 74(366), 427-431.
#' \doi{10.2307/2286348}
#'
#' Phillips, P. C. B. and Perron, P. (1988). Testing for a unit root in time
#' series regression. \emph{Biometrika}, 75(2), 335-346.
#' \doi{10.1093/biomet/75.2.335}
#'
#' Kwiatkowski, D., Phillips, P. C. B., Schmidt, P., and Shin, Y. (1992).
#' Testing the null hypothesis of stationarity against the alternative of a
#' unit root. \emph{Journal of Econometrics}, 54(1-3), 159-178.
#' \doi{10.1016/0304-4076(92)90104-Y}
#'
#' Elliott, G., Rothenberg, T. J., and Stock, J. H. (1996). Efficient tests
#' for an autoregressive unit root. \emph{Econometrica}, 64(4), 813-836.
#' \doi{10.2307/2171846}
#'
#' Zivot, E. and Andrews, D. W. K. (1992). Further evidence on the great
#' crash, the oil-price shock, and the unit-root hypothesis.
#' \emph{Journal of Business & Economic Statistics}, 10(3), 251-270.
#' \doi{10.1080/07350015.1992.10509904}
#'
#' @examples
#' set.seed(42)
#' x <- cumsum(rnorm(60))
#' res <- urstat(x, tests = c("ADF", "PP", "KPSS"), strategy = FALSE)
#'
#' @importFrom stats lm.fit setNames
#' @importFrom tseries adf.test pp.test kpss.test
#' @importFrom urca ur.ers
#' @importFrom strucchange Fstats
#' @export
urstat <- function(x,
                   tests     = "ALL",
                   max_lag   = 12L,
                   crit      = "BIC",
                   pp_lag    = 4L,
                   kpss_lags = 8L,
                   za_trim   = 0.15,
                   level     = 0.05,
                   stars     = TRUE,
                   strategy  = TRUE) {

  ## ---------- Normalise input ------------------------------------------
  series_list <- .urs_normalise_input(x)
  crit  <- toupper(crit)
  tests <- toupper(tests)
  if (identical(tests, "ALL")) {
    tests <- c("ADF", "PP", "KPSS", "ERS", "ZA", "KM")
  }
  do_adf  <- "ADF"  %in% tests
  do_pp   <- "PP"   %in% tests
  do_kpss <- "KPSS" %in% tests
  do_ers  <- "ERS"  %in% tests
  do_za   <- "ZA"   %in% tests
  do_km   <- "KM"   %in% tests

  all_results <- list()

  for (nm in names(series_list)) {
    y <- series_list[[nm]]
    if (!is.numeric(y) || length(y) < 8L) {
      message("urstat: skipping '", nm, "' (too few observations or non-numeric).")
      next
    }
    dy  <- diff(y)
    d2y <- diff(dy)

    res <- list()

    ## ---- ADF -----------------------------------------------------------
    if (do_adf) {
      res$adf <- list(
        level = list(
          nc = .urs_adf(y, det = "nc", max_lag = max_lag, crit = crit),
          c  = .urs_adf(y, det = "c",  max_lag = max_lag, crit = crit),
          ct = .urs_adf(y, det = "ct", max_lag = max_lag, crit = crit)
        ),
        d1 = list(
          nc = .urs_adf(dy, det = "nc", max_lag = max_lag, crit = crit),
          c  = .urs_adf(dy, det = "c",  max_lag = max_lag, crit = crit),
          ct = .urs_adf(dy, det = "ct", max_lag = max_lag, crit = crit)
        ),
        d2 = list(
          nc = .urs_adf(d2y, det = "nc", max_lag = max_lag, crit = crit),
          c  = .urs_adf(d2y, det = "c",  max_lag = max_lag, crit = crit),
          ct = .urs_adf(d2y, det = "ct", max_lag = max_lag, crit = crit)
        )
      )
    }

    ## ---- PP ------------------------------------------------------------
    if (do_pp) {
      res$pp <- list(
        level = list(
          c  = .urs_pp(y,   det = "constant", lags = pp_lag),
          ct = .urs_pp(y,   det = "trend",    lags = pp_lag)
        ),
        d1 = list(
          c  = .urs_pp(dy,  det = "constant", lags = pp_lag),
          ct = .urs_pp(dy,  det = "trend",    lags = pp_lag)
        ),
        d2 = list(
          c  = .urs_pp(d2y, det = "constant", lags = pp_lag),
          ct = .urs_pp(d2y, det = "trend",    lags = pp_lag)
        )
      )
    }

    ## ---- KPSS ----------------------------------------------------------
    if (do_kpss) {
      res$kpss <- list(
        level = list(
          c  = .urs_kpss(y,   type = "mu",  lags = kpss_lags),
          ct = .urs_kpss(y,   type = "tau", lags = kpss_lags)
        ),
        d1 = list(
          c  = .urs_kpss(dy,  type = "mu",  lags = kpss_lags),
          ct = .urs_kpss(dy,  type = "tau", lags = kpss_lags)
        ),
        d2 = list(
          c  = .urs_kpss(d2y, type = "mu",  lags = kpss_lags),
          ct = .urs_kpss(d2y, type = "tau", lags = kpss_lags)
        )
      )
    }

    ## ---- ERS / DF-GLS --------------------------------------------------
    if (do_ers) {
      res$ers <- list(
        level = list(
          c  = .urs_ers(y,   det = "constant"),
          ct = .urs_ers(y,   det = "trend")
        ),
        d1 = list(
          c  = .urs_ers(dy,  det = "constant"),
          ct = .urs_ers(dy,  det = "trend")
        ),
        d2 = list(
          c  = .urs_ers(d2y, det = "constant"),
          ct = .urs_ers(d2y, det = "trend")
        )
      )
    }

    ## ---- Zivot-Andrews -------------------------------------------------
    if (do_za) {
      res$za <- list(
        level = .urs_za(y,   trim = za_trim),
        d1    = .urs_za(dy,  trim = za_trim),
        d2    = .urs_za(d2y, trim = za_trim)
      )
    }

    ## ---- KM (Kobayashi-McAleer) ----------------------------------------
    if (do_km) {
      res$km <- list(
        level = .urs_km(y),
        d1    = .urs_km(dy),
        d2    = .urs_km(d2y)
      )
    }

    ## ---- Elder-Kennedy decision ----------------------------------------
    dec <- .urs_decision(res, level = level)
    res$decision <- dec$order
    res$process  <- dec$process

    all_results[[nm]] <- res

    ## ---- Print tables --------------------------------------------------
    .urs_print_table1(nm, res, do_adf = do_adf, do_pp = do_pp,
                      do_kpss = do_kpss, stars = stars)
    if (do_ers || do_za) {
      .urs_print_table2(nm, res, do_ers = do_ers, do_za = do_za, stars = stars)
    }
    if (strategy) {
      .urs_print_strategy(nm, res, stars = stars)
    }
  }

  invisible(all_results)
}


## ==========================================================================
## Internal helpers
## ==========================================================================

#' @keywords internal
.urs_normalise_input <- function(x) {
  if (is.numeric(x)) {
    nm <- deparse(substitute(x))
    if (nchar(nm) > 20) nm <- "series"
    return(stats::setNames(list(as.numeric(x)), nm))
  }
  if (is.data.frame(x)) {
    num_cols <- vapply(x, is.numeric, logical(1L))
    return(as.list(x[, num_cols, drop = FALSE]))
  }
  if (is.list(x)) {
    if (is.null(names(x))) names(x) <- paste0("series", seq_along(x))
    return(x)
  }
  stop("'x' must be a numeric vector, ts, data.frame, or named list.")
}

#' @keywords internal
.urs_adf_lag <- function(y, det, max_lag, crit) {
  if (length(y) < max_lag + 5L) max_lag <- max(0L, length(y) - 5L)
  best_ic  <- Inf
  best_lag <- 0L
  dy <- diff(y)
  n  <- length(dy)
  for (k in 0L:max_lag) {
    idx <- (k + 2L):n
    if (length(idx) < 3L) break
    y_dep  <- dy[idx]
    y_lag1 <- y[idx]            # y_{t-1} in levels
    mat <- matrix(y_lag1, ncol = 1L)
    if (det %in% c("c", "ct"))  mat <- cbind(1, mat)
    if (det == "ct") mat <- cbind(mat, seq_along(idx))
    if (k > 0L) {
      for (j in 1L:k) mat <- cbind(mat, dy[idx - j])
    }
    fit <- tryCatch(stats::lm.fit(mat, y_dep), error = function(e) NULL)
    if (is.null(fit)) next
    rss <- sum(fit$residuals^2)
    np  <- ncol(mat)
    ic  <- if (crit == "AIC") n * log(rss / n) + 2 * np else
      n * log(rss / n) + log(n) * np
    if (ic < best_ic) { best_ic <- ic; best_lag <- k }
  }
  best_lag
}

#' @keywords internal
.urs_adf <- function(y, det, max_lag, crit) {
  lag <- .urs_adf_lag(y, det, max_lag, crit)
  type <- switch(det, nc = "none", c = "drift", ct = "trend")
  out  <- tryCatch(tseries::adf.test(y, k = lag), error = function(e) NULL)
  if (is.null(out)) return(list(stat = NA_real_, pval = NA_real_, lag = lag))
  list(stat = unname(out$statistic), pval = out$p.value, lag = lag)
}

#' @keywords internal
.urs_pp <- function(y, det, lags) {
  out <- tryCatch(tseries::pp.test(y, lshort = (lags <= 4L)),
                  error = function(e) NULL)
  if (is.null(out)) return(list(stat = NA_real_, pval = NA_real_))
  list(stat = unname(out$statistic), pval = out$p.value)
}

#' @keywords internal
.urs_kpss <- function(y, type, lags) {
  null_type <- if (type == "mu") "Level" else "Trend"
  out <- tryCatch(tseries::kpss.test(y, null = null_type, lshort = TRUE),
                  error = function(e) NULL)
  if (is.null(out)) return(list(stat = NA_real_, pval = NA_real_))
  list(stat = unname(out$statistic), pval = out$p.value)
}

#' @keywords internal
.urs_ers <- function(y, det) {
  type <- if (det == "trend") "trend" else "drift"
  out  <- tryCatch(urca::ur.ers(y, type = "DF-GLS", model = type,
                                lag.max = 4L),
                   error = function(e) NULL)
  if (is.null(out)) return(list(stat = NA_real_, cv5 = NA_real_))
  stat <- out@teststat
  cv5  <- out@cval[, "5pct"]
  list(stat = stat, cv5 = cv5)
}

#' @keywords internal
.urs_za <- function(y, trim) {
  out <- tryCatch(
    strucchange::Fstats(y ~ stats::lag(y, -1L) + 1L,
                        from = trim, to = 1 - trim),
    error = function(e) NULL
  )
  if (is.null(out)) return(list(stat = NA_real_, breakpoint = NA_integer_))
  bp <- which.max(out$Fstats)
  list(stat = max(out$Fstats, na.rm = TRUE), breakpoint = bp)
}

#' @keywords internal
.urs_km <- function(y) {
  if (any(y <= 0, na.rm = TRUE)) {
    return(list(stat = NA_real_, pval = NA_real_, recommendation = "non-positive values"))
  }
  ly <- log(y)
  n  <- length(y)
  t  <- seq_len(n)
  fit_lin <- tryCatch(stats::lm(y  ~ t), error = function(e) NULL)
  fit_log <- tryCatch(stats::lm(ly ~ t), error = function(e) NULL)
  if (is.null(fit_lin) || is.null(fit_log)) {
    return(list(stat = NA_real_, pval = NA_real_, recommendation = NA_character_))
  }
  rss_lin <- sum(fit_lin$residuals^2)
  rss_log <- sum(fit_log$residuals^2)
  stat    <- (rss_lin - rss_log) / rss_log
  pval    <- 2 * stats::pt(abs(stat), df = n - 2L, lower.tail = FALSE)
  rec     <- if (!is.na(pval) && pval < 0.05) "LOGS" else "LEVELS"
  list(stat = stat, pval = pval, recommendation = rec)
}

#' @keywords internal
.urs_stars <- function(pval) {
  if (is.null(pval) || is.na(pval)) return("")
  if (pval < 0.01) return("***")
  if (pval < 0.05) return("**")
  if (pval < 0.10) return("*")
  ""
}

#' @keywords internal
.urs_fmt_stat <- function(stat, pval, stars) {
  if (is.na(stat)) return("---")
  s <- sprintf("%.4f", stat)
  if (stars) s <- paste0(s, .urs_stars(pval))
  s
}

#' @keywords internal
.urs_print_table1 <- function(nm, res, do_adf, do_pp, do_kpss, stars) {
  line <- strrep("-", 72)
  message(line)
  message("  Table 1. Unit Root Tests: ", nm)
  message(line)
  hdr <- sprintf("  %-14s | %-14s | %-14s | %-14s",
                 "Test/Det", "Level", "1st Diff", "2nd Diff")
  message(hdr)
  message(line)

  specs <- c("c (const)", "ct (c+t)")

  if (do_adf && !is.null(res$adf)) {
    message("  ADF:")
    for (sp in c("c", "ct")) {
      lab <- if (sp == "c") "  Const" else "  Const+Trend"
      get_adf <- function(lv) {
        r <- res$adf[[lv]][[sp]]
        .urs_fmt_stat(r$stat, r$pval, stars)
      }
      message(sprintf("  %-13s | %-14s | %-14s | %-14s",
                      lab, get_adf("level"), get_adf("d1"), get_adf("d2")))
    }
    message(line)
  }

  if (do_pp && !is.null(res$pp)) {
    message("  PP:")
    for (sp in c("c", "ct")) {
      lab <- if (sp == "c") "  Const" else "  Const+Trend"
      get_pp <- function(lv) {
        r <- res$pp[[lv]][[sp]]
        .urs_fmt_stat(r$stat, r$pval, stars)
      }
      message(sprintf("  %-13s | %-14s | %-14s | %-14s",
                      lab, get_pp("level"), get_pp("d1"), get_pp("d2")))
    }
    message(line)
  }

  if (do_kpss && !is.null(res$kpss)) {
    message("  KPSS (null = stationary):")
    for (sp in c("c", "ct")) {
      lab <- if (sp == "c") "  Const" else "  Const+Trend"
      get_kp <- function(lv) {
        r <- res$kpss[[lv]][[sp]]
        .urs_fmt_stat(r$stat, r$pval, stars)
      }
      message(sprintf("  %-13s | %-14s | %-14s | %-14s",
                      lab, get_kp("level"), get_kp("d1"), get_kp("d2")))
    }
    message(line)
  }

  if (stars) message("  ***, **, * = 1%, 5%, 10% significance")
  message(line)
}

#' @keywords internal
.urs_print_table2 <- function(nm, res, do_ers, do_za, stars) {
  line <- strrep("-", 72)
  message(line)
  message("  Table 2. Advanced Tests: ", nm)
  message(line)

  if (do_ers && !is.null(res$ers)) {
    message("  ERS/DF-GLS (null = unit root):")
    for (sp in c("c", "ct")) {
      lab <- if (sp == "c") "  Const" else "  Trend"
      get_ers <- function(lv) {
        r <- res$ers[[lv]][[sp]]
        if (is.na(r$stat)) return("---")
        sprintf("%.4f", r$stat)
      }
      message(sprintf("  %-13s | Level: %-8s | d1: %-8s | d2: %-8s",
                      lab, get_ers("level"), get_ers("d1"), get_ers("d2")))
    }
    message(line)
  }

  if (do_za && !is.null(res$za)) {
    message("  Zivot-Andrews (null = unit root):")
    for (lv in c("level", "d1", "d2")) {
      r <- res$za[[lv]]
      lbl <- switch(lv, level = "Level", d1 = "1st Diff", d2 = "2nd Diff")
      stat_str <- if (is.na(r$stat)) "---" else sprintf("%.4f", r$stat)
      bp_str   <- if (is.na(r$breakpoint)) "---" else as.character(r$breakpoint)
      message(sprintf("  %-14s  F-stat: %-10s  Break obs: %s",
                      lbl, stat_str, bp_str))
    }
    message(line)
  }
}

#' @keywords internal
.urs_decision <- function(res, level) {
  get_p <- function(lst, path) {
    tryCatch({
      r <- lst
      for (nm in path) {
        r <- r[[nm]]
        if (is.null(r)) return(NA_real_)
      }
      pv <- r$pval
      if (is.null(pv)) NA_real_ else pv
    }, error = function(e) NA_real_)
  }

  p_lct <- get_p(res, c("adf", "level", "ct"))
  p_dc  <- get_p(res, c("adf", "d1",    "c"))
  p_d2c <- get_p(res, c("adf", "d2",    "c"))

  trend_sig <- FALSE   # simplified: no trend test here

  if (!is.na(p_lct) && !is.null(p_lct) && p_lct < level) {
    return(list(order = "I(0)", process = "Stationary"))
  }
  if (!is.na(p_dc) && !is.null(p_dc) && p_dc < level) {
    return(list(order = "I(1)", process = "Difference-stationary"))
  }
  if (!is.na(p_d2c) && !is.null(p_d2c) && p_d2c < level) {
    return(list(order = "I(2)", process = "Difference-stationary (2x)"))
  }
  list(order = "I(>2)", process = "Possibly non-stationary")
}

#' @keywords internal
.urs_print_strategy <- function(nm, res, stars) {
  line <- strrep("-", 72)
  message(line)
  message("  Elder-Kennedy Decision: ", nm)
  message("  Integration order : ", res$decision)
  message("  Suggested process : ", res$process)
  message(line)
}
