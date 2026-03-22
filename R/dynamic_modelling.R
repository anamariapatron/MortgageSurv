# ==============================================================================
#  dynamic_modelling.R
#  Landmarking (LOCF) model fitting for mortgage default prediction.
#  Wraps the Landmarking package within the mortgage-default context.
# ==============================================================================


#' Prepare longitudinal data for landmark modelling
#'
#' Converts the wide-format simulation output into the long format required by
#' \code{Landmarking::return_ids_with_LOCF()} and
#' \code{Landmarking::fit_LOCF_landmark()}.
#'
#' @param time_matrix Numeric matrix (n x m) of simulated event times.
#' @param status_matrix Integer matrix (n x m) of event statuses.
#' @param x1_matrix Numeric matrix of interest-rate covariate values.
#' @param x2_matrix Numeric matrix of inflation covariate values.
#' @param x3_matrix Numeric matrix of loan-to-income covariate values.
#' @param x4_matrix Numeric matrix of borrower-age covariate values.
#' @param obsdate_matrix Date matrix of observation dates.
#' @param credit_start_dates Date vector of mortgage origination dates (length n).
#' @param landmarks Numeric vector of landmark time-points (in months from
#'   origination). Must have the same length as \code{ncol(time_matrix)}.
#'   Example: \code{c(54, 56, 58, 60, 62, 64)}.
#'
#' @return A \code{data.frame} in long format containing columns
#'   \code{id}, \code{iteration}, \code{time_raw}, \code{status_raw},
#'   \code{x1}–\code{x4}, \code{obs_date}, \code{obs_time},
#'   \code{month_final}, \code{event_time}, \code{event_status}.
#'
#' @importFrom lubridate interval months
#' @importFrom dplyr group_by mutate ungroup filter if_else any
#' @export
prepare_landmark_data <- function(time_matrix,
                                  status_matrix,
                                  x1_matrix,
                                  x2_matrix,
                                  x3_matrix,
                                  x4_matrix,
                                  obsdate_matrix,
                                  credit_start_dates,
                                  landmarks) {

  n <- nrow(time_matrix)
  m <- ncol(time_matrix)

  stopifnot(length(landmarks) == m)

  df_long <- data.frame()

  for (i in seq_len(m)) {
    df_long <- rbind(df_long,
                     data.frame(
                       id          = seq_len(n),
                       iteration   = i,
                       time_raw    = time_matrix[, i],
                       status_raw  = status_matrix[, i],
                       x1          = x1_matrix[, i],
                       x2          = x2_matrix[, i],
                       x3          = x3_matrix[, i],
                       x4          = x4_matrix[, i],
                       obs_date    = obsdate_matrix[, i]
                     ))
  }

  df_long$obs_date <- as.Date(df_long$obs_date)
  df_long$obs_time <- as.numeric(
    difftime(df_long$obs_date, credit_start_dates[df_long$id], units = "days")
  ) / 30.44
  df_long$month_final <- landmarks[df_long$iteration]

  # Remove rows where individual has already defaulted (NA status)
  df_long <- df_long[!is.na(df_long$status_raw), ]

  # Compute per-individual event time (first default) and status flag
  df_long <- df_long |>
    (\(d) {
      split_d <- split(d, d$id)
      lapply(split_d, function(sub) {
        if (any(sub$status_raw == 1)) {
          sub$event_time   <- sub$time_raw[sub$status_raw == 1][1]
          sub$event_status <- 1L
        } else {
          sub$event_time   <- max(sub$time_raw, na.rm = TRUE)
          sub$event_status <- 0L
        }
        sub
      }) |> do.call(what = rbind)
    })()

  rownames(df_long) <- NULL
  df_long
}


#' Fit a dynamic Landmarking model (LOCF) to mortgage default data
#'
#' Wrapper around \code{Landmarking::return_ids_with_LOCF()} and
#' \code{Landmarking::fit_LOCF_landmark()} using a standard Cox
#' proportional-hazards sub-model.
#'
#' @param df_long A \code{data.frame} as returned by
#'   \code{\link{prepare_landmark_data}}.
#' @param x_L Numeric vector of landmark time-points (months).
#' @param x_hor Numeric vector of horizon time-points (one per landmark).
#'   Typically \code{x_L + 12} for a 12-month prediction horizon.
#' @param covariates Character vector of covariate column names.
#'   Default \code{c("x1", "x2", "x3", "x4")}.
#' @param survival_submodel Character. Sub-model type passed to
#'   \code{fit_LOCF_landmark}. Default \code{"standard_cox"}.
#' @param k Integer. Minimum number of observations required per landmark.
#'   Default \code{1}.
#'
#' @return A named list (one element per landmark) as returned by
#'   \code{Landmarking::fit_LOCF_landmark()}, where each element contains
#'   the fitted model and individual-level predicted default probabilities.
#'
#' @examples
#' \dontrun{
#' sim     <- simulate_mortgage_data(n = 500, seed = 42)
#' df_long <- prepare_landmark_data(
#'   time_matrix        = sim$time_matrix,
#'   status_matrix      = sim$status_matrix,
#'   x1_matrix          = sim$x1_matrix,
#'   x2_matrix          = sim$x2_matrix,
#'   x3_matrix          = sim$x3_matrix,
#'   x4_matrix          = sim$x4_matrix,
#'   obsdate_matrix     = sim$obsdate_matrix,
#'   credit_start_dates = sim$credit_start_dates,
#'   landmarks          = c(54, 56, 58, 60, 62, 64)
#' )
#' fit <- fit_landmark_model(df_long, x_L = c(54, 56, 58, 60, 62, 64))
#' }
#'
#' @importFrom Landmarking return_ids_with_LOCF fit_LOCF_landmark
#' @export
fit_landmark_model <- function(df_long,
                                x_L,
                                x_hor           = x_L + 12,
                                covariates       = c("x1", "x2", "x3", "x4"),
                                survival_submodel = "standard_cox",
                                k                = 1L) {

  covariates_time <- rep("obs_time", length(covariates))

  df_locf <- Landmarking::return_ids_with_LOCF(
    data_long       = df_long,
    individual_id   = "id",
    covariates      = covariates,
    covariates_time = covariates_time,
    x_L             = x_L
  )

  Landmarking::fit_LOCF_landmark(
    data_long         = df_locf,
    x_L               = x_L,
    x_hor             = x_hor,
    covariates        = covariates,
    covariates_time   = covariates_time,
    k                 = k,
    individual_id     = "id",
    event_time        = "event_time",
    event_status      = "event_status",
    survival_submodel = survival_submodel
  )
}


#' Plot the cumulative hazard trajectory for a specific borrower
#'
#' For a given individual, extracts the predicted default probability at each
#' landmark from the fitted landmark model and converts it to a cumulative
#' hazard via \eqn{H = -\log(1 - p)}.
#'
#' @param landmark_model A fitted landmark model list as returned by
#'   \code{\link{fit_landmark_model}}.
#' @param person_id Integer. The individual identifier to plot.
#'
#' @return A \code{ggplot2} object (invisible). Also printed to the current device.
#'
#' @examples
#' \dontrun{
#' plot_cumulative_hazard(fit, person_id = 3)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_minimal
#' @importFrom dplyr filter bind_rows
#' @export
plot_cumulative_hazard <- function(landmark_model, person_id = 1L) {
  LMs <- names(landmark_model)

  cumhaz_person <- lapply(LMs, function(LM) {
    df_lm    <- landmark_model[[LM]]$data
    df_p     <- df_lm[df_lm$id == person_id, , drop = FALSE]
    if (nrow(df_p) == 0) return(NULL)
    data.frame(
      LM     = as.numeric(LM),
      cumhaz = -log(1 - df_p$event_prediction)
    )
  })
  cumhaz_person <- do.call(rbind, Filter(Negate(is.null), cumhaz_person))

  p <- ggplot2::ggplot(cumhaz_person, ggplot2::aes(x = LM, y = cumhaz)) +
    ggplot2::geom_point(size = 3, color = "blue") +
    ggplot2::geom_line(color = "blue", linewidth = 1) +
    ggplot2::labs(
      x     = "Landmark time (months)",
      y     = "Cumulative hazard",
      title = paste("Cumulative hazard – borrower", person_id)
    ) +
    ggplot2::theme_minimal()
  print(p)
  invisible(p)
}
