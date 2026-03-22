# ==============================================================================
#  static_modelling.R
#  Static (single time-point) survival models for mortgage default.
#  Covers data preparation for AFT/PH models, model fitting (Cox, Weibull PH,
#  Weibull AFT via flexsurv), and cumulative-hazard visualisation.
# ==============================================================================


#' Prepare a one-row-per-individual dataset for static survival models
#'
#' Converts the wide-format simulation output into the cross-sectional format
#' required by Cox, Weibull PH, and Weibull AFT models. Each individual
#' contributes exactly one row: either the first default event observed, or
#' (if no event) the covariate values from the last landmark period.
#'
#' @param time_matrix Numeric matrix (n × m) of simulated event times.
#' @param status_matrix Integer matrix (n × m) of event statuses
#'   (1 = default, 0 = censored, \code{NA} = already defaulted).
#' @param x1_matrix Numeric matrix of interest-rate covariate values.
#' @param x2_matrix Numeric matrix of inflation covariate values.
#' @param x3_matrix Numeric matrix of loan-to-income covariate values.
#' @param x4_matrix Numeric matrix of borrower-age covariate values.
#' @param obsdate_matrix Date matrix of observation dates.
#' @param credit_start_dates Date vector of mortgage origination dates (length n).
#' @param landmarks Numeric vector of landmark time-points (months).
#'   Default \code{c(60, 62, 64, 66, 68, 70, 72)}.
#' @param censor_time Numeric. Time (months) assigned to censored individuals
#'   with no event. Default \code{72}.
#'
#' @return A \code{data.frame} with one row per borrower and columns:
#'   \code{id}, \code{status}, \code{time}, \code{iteration},
#'   \code{x1}–\code{x4}, \code{month_final}.
#'
#' @details
#' For borrowers who defaulted, \code{time} is the observed time (in months
#' from origination) at the iteration of first default, and \code{status = 1}.
#' For censored borrowers, \code{time = censor_time} and \code{status = 0},
#' with covariates taken from the final landmark iteration.
#'
#' @examples
#' \dontrun{
#' sim    <- simulate_mortgage_data(n = 500, seed = 42)
#' df_aft <- prepare_aft_data(
#'   time_matrix        = sim$time_matrix,
#'   status_matrix      = sim$status_matrix,
#'   x1_matrix          = sim$x1_matrix,
#'   x2_matrix          = sim$x2_matrix,
#'   x3_matrix          = sim$x3_matrix,
#'   x4_matrix          = sim$x4_matrix,
#'   obsdate_matrix     = sim$obsdate_matrix,
#'   credit_start_dates = sim$credit_start_dates
#' )
#' head(df_aft)
#' }
#'
#' @importFrom lubridate as_date
#' @export
prepare_aft_data <- function(time_matrix,
                              status_matrix,
                              x1_matrix,
                              x2_matrix,
                              x3_matrix,
                              x4_matrix,
                              obsdate_matrix,
                              credit_start_dates,
                              landmarks    = c(60, 62, 64, 66, 68, 70, 72),
                              censor_time  = 72) {

  n <- nrow(time_matrix)
  m <- ncol(time_matrix)

  # ── Build long-format data frame ────────────────────────────────────────────
  df_long <- data.frame()
  for (i in seq_len(m)) {
    df_long <- rbind(df_long,
                     data.frame(
                       id        = seq_len(n),
                       iteration = i,
                       time      = time_matrix[, i],
                       status    = status_matrix[, i],
                       x1        = x1_matrix[, i],
                       x2        = x2_matrix[, i],
                       x3        = x3_matrix[, i],
                       x4        = x4_matrix[, i],
                       obs_date  = obsdate_matrix[, i]
                     ))
  }

  df_long$obs_date   <- as.Date(df_long$obs_date)
  df_long$obs_time   <- as.numeric(
    difftime(df_long$obs_date, credit_start_dates[df_long$id], units = "days")
  ) / 30.44
  df_long$month_final <- landmarks[df_long$iteration]

  # ── Collapse to one row per individual ──────────────────────────────────────
  ids  <- seq_len(n)
  out  <- lapply(ids, function(pid) {
    sub         <- df_long[df_long$id == pid, ]
    event_rows  <- which(!is.na(sub$status) & sub$status == 1)
    last_iter   <- max(sub$iteration)

    if (length(event_rows) > 0) {
      first_ev <- event_rows[1]
      data.frame(
        id          = pid,
        status      = 1L,
        time        = sub$obs_time[first_ev],
        iteration   = sub$iteration[first_ev],
        x1          = sub$x1[first_ev],
        x2          = sub$x2[first_ev],
        x3          = sub$x3[first_ev],
        x4          = sub$x4[first_ev],
        month_final = sub$month_final[first_ev]
      )
    } else {
      last_sub <- sub[sub$iteration == last_iter, ][1, ]
      data.frame(
        id          = pid,
        status      = 0L,
        time        = censor_time,
        iteration   = last_iter,
        x1          = last_sub$x1,
        x2          = last_sub$x2,
        x3          = last_sub$x3,
        x4          = last_sub$x4,
        month_final = censor_time
      )
    }
  })

  do.call(rbind, out)
}


#' Fit static survival models for mortgage default
#'
#' Fits three complementary static survival models on a one-row-per-individual
#' dataset:
#' \enumerate{
#'   \item **Cox** semi-parametric proportional hazards (\code{survival::coxph}).
#'   \item **Weibull PH** fully parametric PH model (\code{flexsurv::flexsurvreg},
#'     \code{dist = "weibullPH"}).
#'   \item **Weibull AFT** accelerated failure time model
#'     (\code{flexsurv::flexsurvreg}, \code{dist = "weibull"}).
#' }
#'
#' @param df_aft A \code{data.frame} as returned by \code{\link{prepare_aft_data}},
#'   with columns \code{time}, \code{status}, \code{x1}, \code{x2}, \code{x3},
#'   \code{x4}.
#' @param formula A \code{formula} object for the survival model.
#'   Default \code{Surv(time, status) ~ x1 + x2 + x3 + x4}.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{cox}{Fitted \code{coxph} object.}
#'   \item{weibull_ph}{Fitted \code{flexsurvreg} object (Weibull PH).}
#'   \item{weibull_aft}{Fitted \code{flexsurvreg} object (Weibull AFT).}
#' }
#'
#' @examples
#' \dontrun{
#' sim    <- simulate_mortgage_data(n = 500, seed = 42)
#' df_aft <- prepare_aft_data(
#'   time_matrix        = sim$time_matrix,
#'   status_matrix      = sim$status_matrix,
#'   x1_matrix          = sim$x1_matrix,
#'   x2_matrix          = sim$x2_matrix,
#'   x3_matrix          = sim$x3_matrix,
#'   x4_matrix          = sim$x4_matrix,
#'   obsdate_matrix     = sim$obsdate_matrix,
#'   credit_start_dates = sim$credit_start_dates
#' )
#' fits <- fit_static_models(df_aft)
#' summary(fits$weibull_aft)
#' }
#'
#' @importFrom survival coxph Surv
#' @importFrom flexsurv flexsurvreg
#' @export
fit_static_models <- function(df_aft,
                               formula = survival::Surv(time, status) ~ x1 + x2 + x3 + x4) {

  cox_fit        <- survival::coxph(formula, data = df_aft)
  fit_weibull_ph  <- flexsurv::flexsurvreg(formula, data = df_aft,
                                            dist = "weibullPH")
  fit_weibull_aft <- flexsurv::flexsurvreg(formula, data = df_aft,
                                            dist = "weibull")

  list(
    cox         = cox_fit,
    weibull_ph  = fit_weibull_ph,
    weibull_aft = fit_weibull_aft
  )
}


#' Plot cumulative hazard functions for static models
#'
#' Overlays the Kaplan–Meier, Weibull PH, and Weibull AFT cumulative hazard
#' curves on a single plot, with optional vertical lines at each landmark.
#'
#' @param fit_weibull_ph A fitted Weibull PH model (from
#'   \code{\link{fit_static_models}} or \code{flexsurv::flexsurvreg}).
#' @param fit_weibull_aft A fitted Weibull AFT model.
#' @param landmarks Numeric vector of landmark time-points to draw as vertical
#'   lines. Default \code{c(60, 62, 64, 66, 68, 70, 72)}.
#' @param xlim Numeric vector of length 2 giving the time axis limits.
#'   Default \code{c(50, 72)}.
#' @param landmark_labels Character vector of labels for the landmark lines.
#'   Default \code{c("t1","t2","t3","t4","t5","t6","t7")}.
#' @param landmark_colors Character vector of colours for the landmark lines.
#'   Defaults to a preset palette.
#'
#' @return Invisible \code{NULL}. Draws to the active graphics device.
#'
#' @examples
#' \dontrun{
#' fits <- fit_static_models(df_aft)
#' plot_cumhaz_static(fits$weibull_ph, fits$weibull_aft)
#' }
#'
#' @importFrom graphics abline legend
#' @export
plot_cumhaz_static <- function(fit_weibull_ph,
                                fit_weibull_aft,
                                landmarks        = c(60, 62, 64, 66, 68, 70, 72),
                                xlim             = c(50, 72),
                                landmark_labels  = paste0("t", seq_along(landmarks)),
                                landmark_colors  = c("#FF9999", "#FFCC33",
                                                     "#2EB67D", "#57C9BE",
                                                     "#7FDBFF", "#9999FF",
                                                     "#FF99B8")) {

  plot(fit_weibull_ph,
       type = "cumhaz",
       col  = "blue",
       lwd  = 2,
       xlab = "Time (months from origination)",
       ylab = "Cumulative Hazard",
       main = "Cumulative Hazard Function – Static Models",
       xlim = xlim)

  plot(fit_weibull_aft,
       type = "cumhaz",
       col  = "green3",
       lwd  = 2,
       add  = TRUE)

  graphics::legend("topleft",
                   legend = c("KM", "Weibull PH", "Weibull AFT"),
                   col    = c("black", "blue", "green3"),
                   lwd    = 2,
                   bty    = "n")

  # Landmark lines
  n_lm <- min(length(landmarks), length(landmark_labels),
               length(landmark_colors))
  for (j in seq_len(n_lm)) {
    graphics::abline(v     = landmarks[j],
                     col   = landmark_colors[j],
                     lty   = 2,
                     lwd   = 1.5)
  }

  invisible(NULL)
}
