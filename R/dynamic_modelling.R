# ==============================================================================
#  dynamic_modelling.R
#  Mirrors 2_DynamicModelling.R EXACTLY — including obsdate_matrix[n, i]
#  (last row for all individuals), which is what the original does because
#  j = n = 1000 is left over from the simulation loop.
# ==============================================================================


#' Prepare longitudinal data for landmark modelling
#'
#' Mirrors \code{2_DynamicModelling.R} line by line.
#' Uses \code{obsdate_matrix[n, i]} (last row) as the common observation
#' date for all individuals per iteration — exactly as the original script.
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
#'   Must have same length as \code{ncol(time_matrix)}.
#'
#' @return A \code{data.frame} ready for \code{\link{fit_landmark_model}}.
#'
#' @importFrom dplyr group_by mutate ungroup if_else filter
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

  # ── Build long format ─────────────────────────────────────────────────────
  # Uses obsdate_matrix[n, i] — last row — for ALL individuals per iteration.
  # This mirrors the original: obs_date = obsdate_matrix[j, i] where j = n.
  df_long <- data.frame()
  for (i in seq_len(m)) {
    df_long <- rbind(df_long,
                     data.frame(
                       id         = seq_len(n),
                       iteration  = i,
                       time_raw   = time_matrix[, i],
                       status_raw = status_matrix[, i],
                       x1         = x1_matrix[, i],
                       x2         = x2_matrix[, i],
                       x3         = x3_matrix[, i],
                       x4         = x4_matrix[, i],
                       obs_date   = obsdate_matrix[n, i]   # last row, all individuals
                     ))
  }

  df_long$obs_date <- as.Date(df_long$obs_date)

  # obs_time = (common_date_i - credit_start_k) / 30.44 months
  # Varies across individuals because credit_start_dates differ
  df_long$obs_time <- as.numeric(
    difftime(df_long$obs_date,
             credit_start_dates[df_long$id],
             units = "days")) / 30.44

  df_long$month_final <- landmarks[df_long$iteration]

  # filter(!is.na(status_raw))
  df_long <- df_long[!is.na(df_long$status_raw), ]

  # event_time
  df_long <- dplyr::group_by(df_long, id)
  df_long <- dplyr::mutate(df_long,
    event_time = if (any(status_raw == 1)) {
      time_raw[status_raw == 1][1]
    } else {
      max(time_raw, na.rm = TRUE)
    }
  )
  df_long <- dplyr::ungroup(df_long)

  # event_status as double (1/0) — Landmarking requires numeric not integer
  df_long <- dplyr::group_by(df_long, id)
  df_long <- dplyr::mutate(df_long,
    event_status = dplyr::if_else(any(status_raw == 1), 1, 0)
  )
  df_long <- dplyr::ungroup(df_long)
  df_long <- as.data.frame(df_long)

  rownames(df_long) <- NULL
  df_long
}


#' Fit a dynamic Landmarking model (LOCF)
#'
#' Mirrors \code{2_DynamicModelling.R} exactly.
#'
#' @param df_long A \code{data.frame} from \code{\link{prepare_landmark_data}}.
#' @param x_L Numeric vector of landmark time-points (months).
#' @param x_hor Numeric vector of horizons. Default \code{x_L + 12}.
#' @param covariates Character vector of covariate names.
#' @param survival_submodel Character. Default \code{"standard_cox"}.
#' @param k Integer. Min observations per landmark. Default \code{1}.
#'
#' @return Named list from \code{Landmarking::fit_LOCF_landmark()}.
#'
#' @examples
#' \dontrun{
#' sim       <- simulate_mortgage_data(n = 1000, n_iterations = 6, seed = 1234)
#' landmarks <- c(54, 56, 58, 60, 62, 64)
#' df_long   <- prepare_landmark_data(
#'   sim$time_matrix, sim$status_matrix,
#'   sim$x1_matrix,   sim$x2_matrix,
#'   sim$x3_matrix,   sim$x4_matrix,
#'   sim$obsdate_matrix, sim$credit_start_dates, landmarks
#' )
#' lm_fit <- fit_landmark_model(df_long, x_L = landmarks)
#' }
#'
#' @importFrom Landmarking return_ids_with_LOCF fit_LOCF_landmark
#' @export
fit_landmark_model <- function(df_long,
                                x_L,
                                x_hor             = x_L + 12,
                                covariates        = c("x1", "x2", "x3", "x4"),
                                survival_submodel = "standard_cox",
                                k                 = 1L) {

  # Igual que el código original — obs_time, no month_final
  covariates_time <- rep("obs_time", length(covariates))

  # Step 1: LOCF
  df_long2 <- Landmarking::return_ids_with_LOCF(
    data_long       = df_long,
    individual_id   = "id",
    covariates      = covariates,
    covariates_time = covariates_time,
    x_L             = x_L
  )

  # Step 2: fit — sin tocar df_long2
  Landmarking::fit_LOCF_landmark(
    data_long         = df_long2,
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


#' Plot cumulative hazard for a specific borrower
#'
#' @param landmark_model Fitted model from \code{\link{fit_landmark_model}}.
#' @param person_id Integer. Borrower ID. Default \code{1}.
#' @return A \code{ggplot2} object (invisible).
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_minimal
#' @export
plot_cumulative_hazard <- function(landmark_model, person_id = 1L) {
  LMs <- names(landmark_model)
  cumhaz_person <- lapply(LMs, function(LM) {
    df_lm     <- landmark_model[[LM]]$data
    df_person <- df_lm[df_lm$id == person_id, , drop = FALSE]
    if (nrow(df_person) == 0) return(NULL)
    data.frame(LM = as.numeric(LM),
               cumhaz = -log(1 - df_person$event_prediction))
  })
  cumhaz_person <- do.call(rbind, Filter(Negate(is.null), cumhaz_person))

  p <- ggplot2::ggplot(cumhaz_person, ggplot2::aes(x = LM, y = cumhaz)) +
    ggplot2::geom_point(size = 3, color = "blue") +
    ggplot2::geom_line(color = "blue", linewidth = 1) +
    ggplot2::labs(x = "Landmark time", y = "Cumulative hazard",
                  title = paste("Cumulative hazard for person", person_id)) +
    ggplot2::theme_minimal()
  print(p)
  invisible(p)
}
