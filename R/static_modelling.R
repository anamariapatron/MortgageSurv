# ==============================================================================
#  static_modelling.R
#  Mirrors 2_StaticModelling.R EXACTLY — including obsdate_matrix[n, i].
# ==============================================================================


#' Prepare one-row-per-borrower dataset for static survival models
#'
#' Mirrors \code{2_StaticModelling.R} exactly.
#' Uses \code{obsdate_matrix[n, i]} (last row) as observation date for all
#' individuals per iteration — same as the original script.
#'
#' @param time_matrix Numeric matrix (n x m).
#' @param status_matrix Integer matrix (n x m).
#' @param x1_matrix Numeric matrix of interest-rate values.
#' @param x2_matrix Numeric matrix of inflation values.
#' @param x3_matrix Numeric matrix of LTI values.
#' @param x4_matrix Numeric matrix of borrower age values.
#' @param obsdate_matrix Date matrix of observation dates.
#' @param credit_start_dates Date vector of origination dates (length n).
#' @param landmarks Numeric vector of landmark months.
#'   Default \code{c(60, 62, 64, 66, 68, 70, 72)}.
#' @param censor_time Numeric. Time assigned to censored borrowers. Default \code{72}.
#' @param last_iteration Integer. Last iteration index for censored borrowers.
#'   Default \code{7}.
#'
#' @return A \code{data.frame} with one row per borrower.
#'
#' @importFrom dplyr group_by summarise mutate ungroup rowwise
#' @export
prepare_aft_data <- function(time_matrix,
                              status_matrix,
                              x1_matrix,
                              x2_matrix,
                              x3_matrix,
                              x4_matrix,
                              obsdate_matrix,
                              credit_start_dates,
                              landmarks     = c(60, 62, 64, 66, 68, 70, 72),
                              censor_time   = 72,
                              last_iteration = 7L) {

  n <- nrow(time_matrix)
  m <- ncol(time_matrix)

  # ── Build long format — uses obsdate_matrix[n, i] for all individuals ─────
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
                       obs_date  = obsdate_matrix[n, i]   # last row, all individuals
                     ))
  }

  df_long$obs_date    <- as.Date(df_long$obs_date)
  df_long$obs_time    <- as.numeric(
    difftime(df_long$obs_date,
             credit_start_dates[df_long$id],
             units = "days")) / 30.44
  df_long$month_final <- landmarks[df_long$iteration]

  # ── Collapse to one row per individual — mirrors original summarise/mutate ─
  df_aft <- df_long |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      event_index = which(status == 1)[1],
      .groups     = "drop"
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      status = ifelse(!is.na(event_index), 1, 0),
      time   = ifelse(!is.na(event_index),
                      df_long$obs_time[df_long$id == id][event_index],
                      censor_time),
      iteration = ifelse(!is.na(event_index),
                         df_long$iteration[df_long$id == id][event_index],
                         last_iteration),
      x1 = ifelse(!is.na(event_index),
                  df_long$x1[df_long$id == id][event_index],
                  df_long$x1[df_long$id == id &
                               df_long$iteration == last_iteration]),
      x2 = ifelse(!is.na(event_index),
                  df_long$x2[df_long$id == id][event_index],
                  df_long$x2[df_long$id == id &
                               df_long$iteration == last_iteration]),
      x3 = ifelse(!is.na(event_index),
                  df_long$x3[df_long$id == id][event_index],
                  df_long$x3[df_long$id == id &
                               df_long$iteration == last_iteration]),
      x4 = ifelse(!is.na(event_index),
                  df_long$x4[df_long$id == id][event_index],
                  df_long$x4[df_long$id == id &
                               df_long$iteration == last_iteration]),
      month_final = ifelse(!is.na(event_index),
                           df_long$month_final[df_long$id == id][event_index],
                           censor_time)
    ) |>
    dplyr::ungroup()

  as.data.frame(df_aft)
}


#' Fit static survival models (Cox, Weibull PH, Weibull AFT)
#'
#' Mirrors \code{2_StaticModelling.R} exactly.
#'
#' @param df_aft A \code{data.frame} from \code{\link{prepare_aft_data}}.
#' @param formula A \code{formula}. Default
#'   \code{Surv(time, status) ~ x1 + x2 + x3 + x4}.
#'
#' @return Named list: \code{cox}, \code{weibull_ph}, \code{weibull_aft}.
#'
#' @importFrom survival coxph Surv
#' @importFrom flexsurv flexsurvreg
#' @export
fit_static_models <- function(df_aft,
                               formula = survival::Surv(time, status) ~
                                           x1 + x2 + x3 + x4) {
  list(
    cox         = survival::coxph(formula, data = df_aft),
    weibull_ph  = flexsurv::flexsurvreg(formula, data = df_aft,
                                         dist = "weibullPH"),
    weibull_aft = flexsurv::flexsurvreg(formula, data = df_aft,
                                         dist = "weibull")
  )
}


#' Plot cumulative hazard from static models
#'
#' @param fit_weibull_ph Fitted Weibull PH model.
#' @param fit_weibull_aft Fitted Weibull AFT model.
#' @param landmarks Numeric vector of landmark lines to draw.
#' @param xlim Numeric vector of length 2. Default \code{c(50, 72)}.
#'
#' @return Invisible \code{NULL}.
#'
#' @importFrom graphics plot legend abline
#' @export
plot_cumhaz_static <- function(fit_weibull_ph,
                                fit_weibull_aft,
                                landmarks = c(60, 62, 64, 66, 68, 70, 72),
                                xlim      = c(50, 72)) {
  colors <- c("#FF9999", "#FFCC33", "#2EB67D",
              "#57C9BE", "#7FDBFF", "#9999FF", "#FF99B8")

  plot(fit_weibull_ph,
       type = "cumhaz", col = "blue", lwd = 2,
       xlab = "Time", ylab = "Cumulative Hazard",
       main = "Cumulative Hazard Function - Static Modelling",
       xlim = xlim)

  plot(fit_weibull_aft,
       type = "cumhaz", col = "green", lwd = 2, add = TRUE)

  graphics::legend("topleft",
                   legend = c("KM", "PH", "AFT"),
                   col    = c("black", "blue", "green"),
                   lwd    = 2, bty = "n")

  n_lm <- min(length(landmarks), length(colors))
  for (j in seq_len(n_lm)) {
    graphics::abline(v = landmarks[j], col = colors[j], lty = 2, lwd = 1.5)
  }
  invisible(NULL)
}
