# ==============================================================================
#  simulation.R
#  Functions for simulating UK mortgage loan data under a Landmark PH Weibull
#  framework. Mortgage origination defines entry into the sample and default is
#  the event of interest.
# ==============================================================================


#' Find the Weibull scale parameter (lambda) for a target default probability
#'
#' Given a Weibull shape parameter \code{k} and an interval \eqn{[t0, t1]},
#' solves numerically for the scale parameter \eqn{\lambda} such that the
#' probability of default within that interval equals \code{p_target}.
#'
#' @param t0 Numeric. Start of the interval (in months).
#' @param t1 Numeric. End of the interval (in months).
#' @param k  Numeric. Weibull shape parameter.
#' @param p_target Numeric. Target default probability within \eqn{[t0, t1]}.
#'   Default is \code{0.2} (20\%).
#'
#' @return Numeric scalar: the scale parameter \eqn{\lambda}.
#'
#' @examples
#' lambda <- find_lambda(t0 = 60, t1 = 62, k = 1.5, p_target = 0.2)
#' lambda
#'
#' @export
find_lambda <- function(t0, t1, k, p_target = 0.2) {
  f <- function(lambda) {
    ratio <- exp(-((t1 / lambda)^k - (t0 / lambda)^k))
    ratio - (1 - p_target)
  }
  stats::uniroot(f, interval = c(1e-6, 1e6))$root
}


#' Generate borrower-specific landmark dates from a credit start date
#'
#' Creates a sequence of landmark dates for a borrower, starting five years
#' after mortgage origination and then spaced every \code{spacing_weeks} weeks.
#' The returned vector has length \code{n_intervals + 1}, so it can be used both
#' as the observation grid (first \code{n_intervals} dates) and as the interval
#' boundaries required by the simulation routine.
#'
#' @param credit_start A \code{Date} object. The origination date of the mortgage.
#' @param n_intervals Integer. Number of landmark intervals. Default \code{7}.
#' @param spacing_weeks Integer. Weeks between consecutive landmark dates after
#'   the first landmark. Default \code{6}.
#'
#' @return A \code{Date} vector of length \code{n_intervals + 1}.
#'
#' @examples
#' lms <- generate_landmarks(as.Date("2015-03-15"), n_intervals = 7)
#' lms
#'
#' @importFrom lubridate years weeks as_date
#' @export
generate_landmarks <- function(credit_start,
                               n_intervals = 7L,
                               spacing_weeks = 6L) {
  if (!inherits(credit_start, "Date")) {
    stop("`credit_start` must be a Date object.")
  }
  if (length(n_intervals) != 1L || !is.numeric(n_intervals) || n_intervals < 1) {
    stop("`n_intervals` must be a positive integer.")
  }
  if (length(spacing_weeks) != 1L || !is.numeric(spacing_weeks) || spacing_weeks <= 0) {
    stop("`spacing_weeks` must be a positive number.")
  }

  first_landmark <- credit_start + lubridate::years(5)
  week_offsets   <- c(0, cumsum(rep(spacing_weeks, as.integer(n_intervals))))
  landmark_dates <- first_landmark + lubridate::weeks(week_offsets)

  lubridate::as_date(landmark_dates)
}


# Convert borrower-specific landmark dates to elapsed months since origination.
landmark_times_in_months <- function(credit_start, landmark_dates) {
  as.numeric(
    lubridate::time_length(
      lubridate::interval(credit_start, landmark_dates),
      unit = "month"
    )
  )
}


#' Simulate event times under a Landmark Proportional Hazards Weibull model
#'
#' For each individual (row in \code{x}), simulates the time to mortgage default
#' using a piecewise Weibull baseline hazard defined over a sequence of numeric
#' landmark times measured in months since origination. The linear predictor
#' scales the baseline via the proportional hazards assumption.
#'
#' @param seed Integer or \code{NULL}. Optional random seed for reproducibility.
#'   If \code{NULL}, the current RNG state is used.
#' @param x Numeric matrix of covariates, one row per individual.
#' @param Betas Numeric matrix of regression coefficients, one row per
#'   landmark interval (\eqn{K} rows, one column per covariate).
#' @param Thetas Numeric matrix of Weibull parameters, one row per interval.
#'   Column 1 is the shape and column 2 is the scale.
#' @param LMs Numeric vector of landmark time-points of length \eqn{K + 1},
#'   measured in months since origination.
#' @param dist Character. Only \code{"W"} (Weibull) is currently supported.
#'
#' @return Numeric vector of length \code{nrow(x)} with simulated event times in
#'   months since origination.
#'
#' @examples
#' \dontrun{
#' x      <- matrix(c(0.03, 0.02, 3.5, 40), nrow = 1)
#' Betas  <- matrix(c(50, 0.9, 0.8, -0.2), nrow = 1)
#' Thetas <- matrix(c(0.5, 22.25), nrow = 1)
#' LMs    <- c(60, 62)
#' t      <- simLMPH(seed = 42, x = x, Betas = Betas,
#'                   Thetas = Thetas, LMs = LMs, dist = "W")
#' }
#'
#' @export
simLMPH <- function(seed = NULL, x, Betas, Thetas, LMs, dist = "W") {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  x <- as.matrix(x)
  LMs <- as.numeric(LMs)

  if (anyNA(LMs) || length(LMs) < 2L) {
    stop("`LMs` must be a numeric vector of landmark times with length at least 2.")
  }
  if (any(diff(LMs) <= 0)) {
    stop("`LMs` must be strictly increasing.")
  }
  if (ncol(Betas) != ncol(x)) {
    stop("`Betas` must have the same number of columns as `x`.")
  }
  if (nrow(Thetas) != nrow(Betas)) {
    stop("`Thetas` and `Betas` must have the same number of rows.")
  }

  K <- length(LMs) - 1L
  if (nrow(Betas) != K) {
    stop("`LMs` must have length nrow(Betas) + 1.")
  }

  if (dist != "W") {
    stop("Only Weibull distribution is currently implemented (dist = 'W').")
  }

  n        <- nrow(x)
  tpoints  <- as.vector(LMs)
  out_time <- rep(NA_real_, n)

  for (k in seq_len(K)) {
    beta_k  <- Betas[k, , drop = FALSE]
    theta_k <- Thetas[k, ]
    linpred <- x %*% t(beta_k)
    mlambda <- exp(drop(linpred))

    S0f    <- function(t) stats::pweibull(t, shape = theta_k[1], scale = theta_k[2],
                                          lower.tail = FALSE)
    quantf <- function(p) stats::qweibull(p, shape = theta_k[1], scale = theta_k[2])

    for (i in seq_len(n)) {
      S_L    <- S0f(tpoints[k])
      u      <- stats::runif(1)
      p_draw <- max(1 - S_L * u^(1 / mlambda[i]), .Machine$double.eps)
      t_rel  <- quantf(p_draw)

      if (is.na(out_time[i]) && t_rel <= tpoints[k + 1]) {
        out_time[i] <- t_rel
      }
    }
  }

  # Extrapolate beyond the last landmark boundary for individuals without an event
  no_event <- which(is.na(out_time))
  if (length(no_event) > 0) {
    beta_k  <- Betas[K, , drop = FALSE]
    theta_k <- Thetas[K, ]
    linpred <- x[no_event, , drop = FALSE] %*% t(beta_k)
    mlambda <- exp(drop(linpred))

    S0f    <- function(t) stats::pweibull(t, shape = theta_k[1], scale = theta_k[2],
                                          lower.tail = FALSE)
    quantf <- function(p) stats::qweibull(p, shape = theta_k[1], scale = theta_k[2])
    S_L    <- S0f(tpoints[K])

    for (j in seq_along(no_event)) {
      u      <- stats::runif(1)
      p_draw <- max(1 - S_L * u^(1 / mlambda[j]), .Machine$double.eps)
      out_time[no_event[j]] <- quantf(p_draw)
    }
  }

  out_time
}


#' Simulate a full UK mortgage portfolio dataset
#'
#' Generates simulated data for \code{n} borrowers over \code{n_iterations}
#' landmark periods. Covariates (interest rate, inflation, LTI, age) are drawn
#' from truncated normal distributions calibrated to the UK mortgage market.
#' Default times are simulated via \code{\link{simLMPH}}.
#'
#' @param n Integer. Number of borrowers to simulate. Default \code{1000}.
#' @param n_iterations Integer. Number of landmark periods. Default \code{7}.
#' @param seed Integer. Base random seed. Default \code{1234}.
#' @param betas Numeric matrix (\eqn{n\_iterations \times 4}) of regression
#'   coefficients for (interest rate, inflation, LTI, age). Defaults to the
#'   values from the original study.
#' @param thetas Numeric matrix (\eqn{n\_iterations \times 2}) of Weibull
#'   (shape, scale) parameters per landmark interval.
#' @param credit_start_year Integer. Year in which mortgages are originated.
#'   Default \code{2015}.
#' @param spacing_weeks Integer. Weeks between consecutive borrower-specific
#'   landmark dates after the first landmark. Default \code{6}.
#'
#' @return A named list with:
#' \describe{
#'   \item{time_matrix}{n x n_iterations matrix of simulated event times.}
#'   \item{status_matrix}{n x n_iterations integer matrix
#'     (1 = default, 0 = censored, NA = already defaulted).}
#'   \item{obsdate_matrix}{n x n_iterations \code{Date} matrix of observation dates.}
#'   \item{x1_matrix}{Interest rate covariate matrix.}
#'   \item{x2_matrix}{Inflation covariate matrix.}
#'   \item{x3_matrix}{Loan-to-income covariate matrix.}
#'   \item{x4_matrix}{Borrower age covariate matrix.}
#'   \item{credit_start_dates}{Length-n \code{Date} vector of origination dates.}
#'   \item{landmarks_list}{List of borrower-specific landmark date vectors of
#'     length \code{n_iterations + 1}.}
#'   \item{landmark_months_list}{List of borrower-specific landmark times in months
#'     since origination.}
#' }
#'
#' @examples
#' \dontrun{
#' sim <- simulate_mortgage_data(n = 200, n_iterations = 7, seed = 42)
#' str(sim)
#' }
#'
#' @importFrom lubridate as_date days years weeks interval time_length months
#' @importFrom truncnorm rtruncnorm
#' @export
simulate_mortgage_data <- function(n                 = 1000L,
                                   n_iterations      = 7L,
                                   seed              = 1234L,
                                   betas             = NULL,
                                   thetas            = NULL,
                                   credit_start_year = 2015L,
                                   spacing_weeks     = 6L) {

  set.seed(seed)

  if (is.null(betas)) {
    betas <- matrix(c(
       50,  0.90,  0.80, -0.20,
       51,  0.92,  0.81, -0.21,
       49,  0.88,  0.79, -0.19,
       52,  0.95,  0.82, -0.22,
       48,  0.89,  0.78, -0.18,
        5,  0.90,  0.80, -0.20,
       51,  0.91,  0.81, -0.21
    ), nrow = 7, byrow = TRUE)
  }

  if (is.null(thetas)) {
    thetas <- matrix(c(
      0.5, 22.25000,
      0.5, 22.49000,
      0.5, 22.49756,
      0.5, 22.73315,
      0.5, 22.96397,
      0.5, 23.19023,
      0.5, 23.41217
    ), nrow = 7, byrow = TRUE)
  }

  if (nrow(betas) != n_iterations) {
    stop("`betas` must have `n_iterations` rows.")
  }
  if (nrow(thetas) != n_iterations) {
    stop("`thetas` must have `n_iterations` rows.")
  }

  base_date <- lubridate::as_date(paste0(credit_start_year, "-01-01"))
  credit_start_dates <- base_date + lubridate::days(sample(0:364, n, replace = TRUE))

  landmarks_list <- lapply(
    credit_start_dates,
    generate_landmarks,
    n_intervals = n_iterations,
    spacing_weeks = spacing_weeks
  )
  landmark_months_list <- Map(
    landmark_times_in_months,
    credit_start = credit_start_dates,
    landmark_dates = landmarks_list
  )

  x1_matrix <- matrix(
    truncnorm::rtruncnorm(n * n_iterations, a = 0.001, b = 0.06, mean = 0.03, sd = 0.01),
    nrow = n, ncol = n_iterations
  )
  x2_matrix <- matrix(
    truncnorm::rtruncnorm(n * n_iterations, a = 0, b = 0.10, mean = 0.02, sd = 0.015),
    nrow = n, ncol = n_iterations
  )
  x3_matrix <- matrix(
    truncnorm::rtruncnorm(n * n_iterations, a = 1, b = 6, mean = 3.5, sd = 2),
    nrow = n, ncol = n_iterations
  )
  x4_matrix <- matrix(
    truncnorm::rtruncnorm(n * n_iterations, a = 20, b = 70, mean = 40, sd = 10),
    nrow = n, ncol = n_iterations
  )

  time_matrix    <- matrix(NA_real_,    nrow = n, ncol = n_iterations)
  status_matrix  <- matrix(NA_integer_, nrow = n, ncol = n_iterations)
  obsdate_matrix <- matrix(as.Date(NA), nrow = n, ncol = n_iterations)

  alive <- rep(TRUE, n)

  for (i in seq_len(n_iterations)) {
    message("Simulating iteration ", i, " of ", n_iterations, " ...")

    covariates_i <- cbind(
      x1_matrix[, i],
      x2_matrix[, i],
      x3_matrix[, i],
      x4_matrix[, i]
    )

    sim_times <- numeric(n)
    obs_dates <- as.Date(rep(NA, n))

    for (j in seq_len(n)) {
      sim_times[j] <- simLMPH(
        seed   = NULL,
        x      = matrix(covariates_i[j, ], nrow = 1),
        Betas  = betas,
        Thetas = thetas,
        LMs    = landmark_months_list[[j]],
        dist   = "W"
      )

      obs_dates[j] <- landmarks_list[[j]][i] + sample(0:3, 1)
    }

    time_matrix[, i]    <- sim_times
    obsdate_matrix[, i] <- obs_dates

    months_elapsed <- as.numeric(
      lubridate::time_length(
        lubridate::interval(credit_start_dates, obs_dates),
        unit = "month"
      )
    )

    status_vals         <- as.integer(months_elapsed >= sim_times)
    status_vals[!alive] <- NA
    status_matrix[, i]  <- status_vals

    alive <- ifelse(is.na(status_vals), FALSE, status_vals == 0)
  }

  list(
    time_matrix         = time_matrix,
    status_matrix       = status_matrix,
    obsdate_matrix      = obsdate_matrix,
    x1_matrix           = x1_matrix,
    x2_matrix           = x2_matrix,
    x3_matrix           = x3_matrix,
    x4_matrix           = x4_matrix,
    credit_start_dates  = credit_start_dates,
    landmarks_list      = landmarks_list,
    landmark_months_list = landmark_months_list
  )
}


# ── Plotting helpers ──────────────────────────────────────────────────────────


#' Plot cumulative mortgage defaults over landmark time-points
#'
#' @param status_matrix Integer matrix of status values from
#'   \code{\link{simulate_mortgage_data}}.
#' @param LMs Optional numeric vector of landmark time-points in months. If
#'   \code{NULL}, a regular sequence starting at month 60 with step 2 is used.
#'
#' @return A \code{ggplot2} object (invisible).
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal
#' @export
plot_cumulative_deaths <- function(status_matrix, LMs = NULL) {
  if (is.null(LMs)) {
    LMs <- seq(from = 60, by = 2, length.out = ncol(status_matrix))
  }
  if (length(LMs) != ncol(status_matrix)) {
    stop("`LMs` must have length equal to ncol(status_matrix).")
  }

  cum_deaths <- cumsum(colSums(status_matrix, na.rm = TRUE))
  df_cum <- data.frame(time = LMs, cum_deaths = cum_deaths)

  p <- ggplot2::ggplot(df_cum, ggplot2::aes(x = time, y = cum_deaths)) +
    ggplot2::geom_line(color = "blue", linewidth = 1) +
    ggplot2::geom_point(color = "darkblue") +
    ggplot2::labs(
      title = "Cumulative Defaults Over Landmark Time-Points",
      x     = "Landmark (months since origination)",
      y     = "Cumulative Defaults"
    ) +
    ggplot2::theme_minimal()

  print(p)
  invisible(p)
}


#' Plot default times per individual
#'
#' @param time_matrix Numeric matrix of event times.
#' @param status_matrix Integer matrix of event statuses.
#'
#' @return A \code{ggplot2} object (invisible).
#'
#' @importFrom ggplot2 ggplot aes geom_rect geom_point labs theme_minimal
#'   element_text theme
#' @export
plot_death_times <- function(time_matrix, status_matrix) {
  df_events <- data.frame(
    individual = rep(seq_len(nrow(status_matrix)), times = ncol(status_matrix)),
    time       = as.vector(time_matrix),
    status     = as.vector(status_matrix)
  )
  df_deaths <- subset(df_events, status == 1 & !is.na(time))

  time_min       <- floor(min(df_deaths$time, na.rm = TRUE))
  time_max       <- ceiling(max(df_deaths$time, na.rm = TRUE))
  quarter_starts <- seq(time_min, time_max, by = 3)
  quarters <- data.frame(start = quarter_starts, end = quarter_starts + 3)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = quarters,
      ggplot2::aes(xmin = start, xmax = end, ymin = 0,
                   ymax = max(df_deaths$individual)),
      fill = "grey90", alpha = 0.3
    ) +
    ggplot2::geom_point(
      data = df_deaths,
      ggplot2::aes(x = time, y = individual),
      color = "red", size = 2, alpha = 0.7
    ) +
    ggplot2::labs(
      title    = "Default Times per Individual",
      subtitle = "Each point represents the time an individual defaulted",
      x        = "Time elapsed since origination (months)",
      y        = "Individual"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 13),
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1)
    )

  print(p)
  invisible(p)
}


#' Plot default dates per individual with landmark lines
#'
#' @param obsdate_matrix Date matrix of observation dates.
#' @param status_matrix Integer matrix of event statuses.
#' @param LM_start_date A \code{Date}. The reference start date for landmark
#'   labels. Default \code{as.Date("2015-01-01")}.
#' @param LMs Optional numeric vector of landmark months. If \code{NULL}, a
#'   regular sequence starting at month 60 with step 2 is used.
#'
#' @return A \code{ggplot2} object (invisible).
#'
#' @importFrom ggplot2 ggplot aes geom_rect geom_point geom_vline geom_text
#'   scale_x_date labs theme_minimal theme element_text guides
#' @importFrom lubridate floor_date ceiling_date "%m+%" months
#' @export
plot_default_dates <- function(obsdate_matrix,
                               status_matrix,
                               LM_start_date = as.Date("2015-01-01"),
                               LMs = NULL) {
  if (is.null(LMs)) {
    LMs <- seq(from = 60, by = 2, length.out = ncol(obsdate_matrix))
  }
  if (length(LMs) != ncol(obsdate_matrix)) {
    stop("`LMs` must have length equal to ncol(obsdate_matrix).")
  }

  df_events <- data.frame(
    individual = rep(seq_len(nrow(status_matrix)), times = ncol(status_matrix)),
    date       = as.Date(as.vector(obsdate_matrix)),
    status     = as.vector(status_matrix)
  )
  df_deaths <- subset(df_events, status == 1 & !is.na(date))

  quarter_starts <- seq(
    lubridate::floor_date(min(df_deaths$date), unit = "quarter"),
    lubridate::ceiling_date(max(df_deaths$date), unit = "quarter"),
    by = "3 months"
  )
  quarters <- data.frame(
    start = quarter_starts,
    end   = quarter_starts + lubridate::months(3)
  )

  mean_dates <- LM_start_date %m+% lubridate::months(LMs)
  labels_df  <- data.frame(
    mean_dates = mean_dates,
    label      = paste0("t", seq_along(mean_dates))
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = quarters,
      ggplot2::aes(xmin = as.Date(start), xmax = as.Date(end),
                   ymin = 0, ymax = max(df_deaths$individual)),
      fill = "grey90", alpha = 0.3
    ) +
    ggplot2::geom_point(
      data = df_deaths,
      ggplot2::aes(x = as.Date(date), y = individual),
      color = "red", size = 2, alpha = 0.7
    ) +
    ggplot2::geom_vline(
      data = labels_df,
      ggplot2::aes(xintercept = mean_dates, color = label),
      linetype = "solid", linewidth = 1
    ) +
    ggplot2::geom_text(
      data = labels_df,
      ggplot2::aes(x = mean_dates, y = max(df_deaths$individual),
                   label = label, color = label),
      vjust = 0.2, hjust = 1, size = 5, fontface = "bold"
    ) +
    ggplot2::scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
    ggplot2::labs(
      title    = "Default Dates per Individual",
      subtitle = "Each point shows the date an individual defaulted",
      x        = "Month",
      y        = "Individual"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle   = ggplot2::element_text(hjust = 0.5, size = 13),
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::guides(color = "none")

  print(p)
  invisible(p)
}
