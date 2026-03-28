# ==============================================================================
# simulation.R
# ==============================================================================

#' Find the Weibull scale parameter (lambda) for a target default probability
#'
#' @param t0 Numeric. Start of the interval (in months).
#' @param t1 Numeric. End of the interval (in months).
#' @param k Numeric. Weibull shape parameter.
#' @param p_target Numeric. Target default probability in [t0, t1].
#'
#' @return Numeric scalar.
#' @export
find_lambda <- function(t0, t1, k, p_target = 0.2) {
  f <- function(lambda) {
    ratio <- exp(-((t1 / lambda)^k - (t0 / lambda)^k))
    ratio - (1 - p_target)
  }
  stats::uniroot(f, interval = c(1e-6, 1e6))$root
}


#' Generate individual landmark dates from a credit start date
#'
#' @param credit_start Date.
#' @param n_iterations Integer. Number of observation landmarks.
#' @param spacing_weeks Integer. Weeks between landmarks after the first.
#'
#' @return Date vector of length n_iterations + 1.
#' @importFrom lubridate years weeks as_date
#' @export
generate_landmarks <- function(credit_start,
                               n_iterations = 7L,
                               spacing_weeks = 6L) {
  lm1 <- credit_start + lubridate::years(5)
  offsets <- c(0, cumsum(rep(spacing_weeks, n_iterations)))
  lubridate::as_date(lm1 + lubridate::weeks(offsets))
}


# Internal helper: convert landmark dates to months since origination
landmark_times_in_months <- function(credit_start, landmark_dates) {
  as.numeric(
    lubridate::time_length(
      lubridate::interval(credit_start, landmark_dates),
      unit = "month"
    )
  )
}


#' Simulate event times under a Landmark PH Weibull model
#'
#' @param seed Integer or NULL.
#' @param x Numeric matrix of covariates.
#' @param Betas Numeric matrix of coefficients.
#' @param Thetas Numeric matrix of Weibull parameters.
#' @param LMs Numeric vector of landmark times in months, length K + 1.
#' @param dist Character. Only "W".
#'
#' @return Numeric vector of event times in months.
#' @export
simLMPH <- function(seed = NULL, x, Betas, Thetas, LMs, dist = "W") {
  if (!is.null(seed)) set.seed(seed)
  
  x <- as.matrix(x)
  LMs <- as.numeric(LMs)
  
  if (anyNA(LMs) || length(LMs) < 2L) {
    stop("`LMs` must be a numeric vector of length at least 2.")
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
    stop("Only Weibull implemented.")
  }
  
  n <- nrow(x)
  tpoints <- as.vector(LMs)
  out_time <- rep(NA_real_, n)
  
  for (k in seq_len(K)) {
    beta_k  <- Betas[k, , drop = FALSE]
    theta_k <- Thetas[k, ]
    linpred <- x %*% t(beta_k)
    mlambda <- exp(drop(linpred))
    
    S0f <- function(t) {
      stats::pweibull(t, shape = theta_k[1], scale = theta_k[2], lower.tail = FALSE)
    }
    quantf <- function(p) {
      stats::qweibull(p, shape = theta_k[1], scale = theta_k[2])
    }
    
    for (i in seq_len(n)) {
      S_L <- S0f(tpoints[k])
      u <- stats::runif(1)
      p_draw <- max(1 - S_L * u^(1 / mlambda[i]), .Machine$double.eps)
      t_rel <- quantf(p_draw)
      
      if (is.na(out_time[i]) && t_rel <= tpoints[k + 1]) {
        out_time[i] <- t_rel
      }
    }
  }
  
  no_event <- which(is.na(out_time))
  if (length(no_event) > 0) {
    beta_k  <- Betas[K, , drop = FALSE]
    theta_k <- Thetas[K, ]
    linpred <- x[no_event, , drop = FALSE] %*% t(beta_k)
    mlambda <- exp(drop(linpred))
    
    S0f <- function(t) {
      stats::pweibull(t, shape = theta_k[1], scale = theta_k[2], lower.tail = FALSE)
    }
    quantf <- function(p) {
      stats::qweibull(p, shape = theta_k[1], scale = theta_k[2])
    }
    S_L <- S0f(tpoints[K])
    
    for (j in seq_along(no_event)) {
      u <- stats::runif(1)
      p_draw <- max(1 - S_L * u^(1 / mlambda[j]), .Machine$double.eps)
      out_time[no_event[j]] <- quantf(p_draw)
    }
  }
  
  out_time
}


#' Simulate a full mortgage dataset
#'
#' @param n Integer. Number of borrowers.
#' @param n_iterations Integer. Number of observation landmarks.
#' @param seed Integer.
#' @param betas Numeric matrix with n_iterations rows and 4 columns.
#' @param thetas Numeric matrix with n_iterations rows and 2 columns.
#' @param credit_start_year Integer.
#' @param spacing_weeks Integer.
#'
#' @return Named list with matrices and dates.
#' @importFrom lubridate as_date days interval time_length
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
      50, 0.90, 0.80, -0.20,
      51, 0.92, 0.81, -0.21,
      49, 0.88, 0.79, -0.19,
      52, 0.95, 0.82, -0.22,
      48, 0.89, 0.78, -0.18,
      5, 0.90, 0.80, -0.20,
      51, 0.91, 0.81, -0.21
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
  
  if (nrow(betas) != n_iterations || ncol(betas) != 4) {
    stop("`betas` must have n_iterations rows and 4 columns.")
  }
  if (nrow(thetas) != n_iterations || ncol(thetas) != 2) {
    stop("`thetas` must have n_iterations rows and 2 columns.")
  }
  
  base_date <- lubridate::as_date(paste0(credit_start_year, "-01-01"))
  credit_start_dates <- base_date + lubridate::days(sample(0:364, n, replace = TRUE))
  
  landmarks_list <- lapply(
    credit_start_dates,
    generate_landmarks,
    n_iterations = n_iterations,
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
    
    status_vals <- as.integer(months_elapsed >= sim_times)
    status_vals[!alive] <- NA
    status_matrix[, i] <- status_vals
    
    alive <- ifelse(is.na(status_vals), FALSE, status_vals == 0)
  }
  
  list(
    time_matrix          = time_matrix,
    status_matrix        = status_matrix,
    obsdate_matrix       = obsdate_matrix,
    x1_matrix            = x1_matrix,
    x2_matrix            = x2_matrix,
    x3_matrix            = x3_matrix,
    x4_matrix            = x4_matrix,
    credit_start_dates   = credit_start_dates,
    landmarks_list       = landmarks_list,
    landmark_months_list = landmark_months_list
  )
}


#' Plot cumulative mortgage defaults
#'
#' @param status_matrix Integer matrix.
#' @param LMs Numeric vector of landmark times.
#' @return ggplot object.
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal
#' @export
plot_cumulative_deaths <- function(status_matrix,
                                   LMs = c(60, 62, 64, 66, 68, 70, 72)) {
  cum_deaths <- cumsum(colSums(status_matrix, na.rm = TRUE))
  df_cum <- data.frame(time = LMs, cum_deaths = cum_deaths)
  
  p <- ggplot2::ggplot(df_cum, ggplot2::aes(x = time, y = cum_deaths)) +
    ggplot2::geom_line(color = "blue", linewidth = 1) +
    ggplot2::geom_point(color = "darkblue") +
    ggplot2::labs(
      title = "Cumulative Defaults Over Time",
      x = "Landmark time (months)",
      y = "Cumulative Defaults"
    ) +
    ggplot2::theme_minimal()
  
  print(p)
  invisible(p)
}


#' Plot default times per individual
#'
#' @param time_matrix Numeric matrix.
#' @param status_matrix Integer matrix.
#' @return ggplot object.
#' @importFrom ggplot2 ggplot aes geom_rect geom_point labs theme_minimal theme element_text
#' @export
plot_death_times <- function(time_matrix, status_matrix) {
  df_events <- data.frame(
    individual = rep(seq_len(nrow(status_matrix)), times = ncol(status_matrix)),
    time       = as.vector(time_matrix),
    status     = as.vector(status_matrix)
  )
  
  df_deaths <- subset(df_events, status == 1 & !is.na(time))
  
  time_min <- floor(min(df_deaths$time, na.rm = TRUE))
  time_max <- ceiling(max(df_deaths$time, na.rm = TRUE))
  quarter_starts <- seq(time_min, time_max, by = 3)
  
  quarters <- data.frame(
    start = quarter_starts,
    end   = quarter_starts + 3
  )
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = quarters,
      ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = max(df_deaths$individual)),
      fill = "grey90", alpha = 0.3
    ) +
    ggplot2::geom_point(
      data = df_deaths,
      ggplot2::aes(x = time, y = individual),
      color = "red", size = 2, alpha = 0.7
    ) +
    ggplot2::labs(
      title = "Default times per individual",
      subtitle = "Each point represents the time an individual defaulted",
      x = "Time elapsed since origination",
      y = "Individual"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 13),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  print(p)
  invisible(p)
}


#' Plot default dates per individual
#'
#' @param obsdate_matrix Date matrix.
#' @param status_matrix Integer matrix.
#' @param LM_start_date Date.
#' @param LMs Numeric vector.
#' @return ggplot object.
#' @importFrom ggplot2 ggplot aes geom_rect geom_point geom_vline geom_text scale_x_date labs theme_minimal theme element_text guides
#' @importFrom lubridate floor_date ceiling_date "%m+%" period
#' @export
plot_default_dates <- function(obsdate_matrix,
                               status_matrix,
                               LM_start_date = as.Date("2015-01-01"),
                               LMs = c(60, 62, 64, 66, 68, 70, 72)) {
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
    end   = quarter_starts + lubridate::period(3, "month")
  )
  
  labels_df <- data.frame(
    mean_dates = LM_start_date %m+% lubridate::period(LMs, "month"),
    label = paste0("t", seq_along(LMs))
  )
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = quarters,
      ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = max(df_deaths$individual)),
      fill = "grey90", alpha = 0.3
    ) +
    ggplot2::geom_point(
      data = df_deaths,
      ggplot2::aes(x = date, y = individual),
      color = "red", size = 2, alpha = 0.7
    ) +
    ggplot2::geom_vline(
      data = labels_df,
      ggplot2::aes(xintercept = mean_dates, color = label),
      linetype = "solid",
      linewidth = 1
    ) +
    ggplot2::geom_text(
      data = labels_df,
      ggplot2::aes(
        x = mean_dates,
        y = max(df_deaths$individual),
        label = label,
        color = label
      ),
      vjust = 0.2, hjust = 1, size = 5, fontface = "bold"
    ) +
    ggplot2::scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
    ggplot2::labs(
      title = "Default Dates per Individual",
      subtitle = "Each point represents the date an individual defaulted",
      x = "Month",
      y = "Individual"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 13),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::guides(color = "none")
  
  print(p)
  invisible(p)
}
