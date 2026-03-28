# ==============================================================================
#  real_data.R
#  Functions to use real mortgage data instead of simulated data.
#  Converts a user-supplied data frame into the matrix format expected by
#  the rest of the MortgageSurv pipeline.
# ==============================================================================


#' Convert a real mortgage dataset into MortgageSurv format
#'
#' Takes a user-supplied \code{data.frame} in long format (one row per
#' borrower per landmark period) and converts it into the list of matrices
#' that all other \pkg{MortgageSurv} functions expect. This is the entry
#' point for working with \strong{real data} instead of simulated data.
#'
#' @param data A \code{data.frame} in long format. Must contain one row per
#'   borrower per landmark period.
#' @param id_col Character. Name of the column identifying each borrower.
#'   Default \code{"id"}.
#' @param time_col Character. Name of the column with the event or censoring
#'   time (numeric, in months from origination). Default \code{"time"}.
#' @param status_col Character. Name of the column with the event indicator
#'   (1 = default, 0 = censored). Default \code{"status"}.
#' @param covariate_cols Character vector. Names of the covariate columns
#'   to include. Default \code{c("x1", "x2", "x3", "x4")}.
#' @param date_col Character or \code{NULL}. Name of the observation date
#'   column (\code{Date} type). If \code{NULL}, observation dates are set to
#'   \code{NA}. Default \code{"obs_date"}.
#' @param start_date_col Character or \code{NULL}. Name of the mortgage
#'   origination date column (\code{Date} type). If \code{NULL}, credit start
#'   dates are set to \code{NA}. Default \code{"credit_start"}.
#' @param iteration_col Character. Name of the column indicating the landmark
#'   period number (integer, starting at 1). Default \code{"iteration"}.
#'
#' @return A named list with the same structure as
#'   \code{\link{simulate_mortgage_data}}, ready to be passed to
#'   \code{\link{prepare_aft_data}} or \code{\link{prepare_landmark_data}}:
#' \describe{
#'   \item{time_matrix}{n x m numeric matrix of event/censoring times.}
#'   \item{status_matrix}{n x m integer matrix of event indicators.}
#'   \item{obsdate_matrix}{n x m \code{Date} matrix of observation dates.}
#'   \item{x1_matrix ... xk_matrix}{One n x m matrix per covariate.}
#'   \item{credit_start_dates}{Length-n \code{Date} vector of origination dates.}
#'   \item{n_individuals}{Number of unique borrowers.}
#'   \item{n_iterations}{Number of landmark periods.}
#'   \item{covariate_names}{Names of the covariate columns used.}
#' }
#'
#' @details
#' ## Expected data format
#'
#' Your \code{data.frame} should look like this:
#'
#' | id | iteration | time | status | x1   | x2   | obs_date   | credit_start |
#' |----|-----------|------|--------|------|------|------------|--------------|
#' | 1  | 1         | 62.3 | 0      | 0.03 | 0.02 | 2020-01-15 | 2015-01-10   |
#' | 1  | 2         | 62.3 | 0      | 0.04 | 0.02 | 2020-03-20 | 2015-01-10   |
#' | 1  | 3         | 62.3 | 1      | 0.04 | 0.03 | 2020-05-25 | 2015-01-10   |
#' | 2  | 1         | 58.1 | 0      | 0.02 | 0.01 | 2020-01-15 | 2015-03-22   |
#' | ...| ...       | ...  | ...    | ...  | ...  | ...        | ...          |
#'
#' Each borrower (\code{id}) should appear once per landmark period
#' (\code{iteration}), until they default or are censored.
#'
#' ## Handling missing landmark periods
#' If a borrower defaulted before the last period, they will have \code{NA}
#' in the status matrix for subsequent periods — consistent with how
#' \code{\link{simulate_mortgage_data}} works.
#'
#' @examples
#' \dontrun{
#' # Example with a real CSV file
#' my_data <- read.csv("mortgage_data.csv")
#'
#' sim <- from_real_data(
#'   data             = my_data,
#'   id_col           = "loan_id",
#'   time_col         = "months_to_event",
#'   status_col       = "defaulted",
#'   covariate_cols   = c("interest_rate", "inflation", "lti", "age"),
#'   date_col         = "obs_date",
#'   start_date_col   = "origination_date",
#'   iteration_col    = "period"
#' )
#'
#' # Now use it exactly like simulate_mortgage_data() output:
#' df_aft  <- prepare_aft_data(sim$time_matrix, sim$status_matrix, ...)
#' df_long <- prepare_landmark_data(...)
#' }
#'
#' @export
from_real_data <- function(data,
                            id_col           = "id",
                            time_col         = "time",
                            status_col       = "status",
                            covariate_cols   = c("x1", "x2", "x3", "x4"),
                            date_col         = "obs_date",
                            start_date_col   = "credit_start",
                            iteration_col    = "iteration") {

  # ── Validate inputs ──────────────────────────────────────────────────────────
  required_cols <- c(id_col, time_col, status_col, iteration_col, covariate_cols)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing from your data: ",
         paste(missing_cols, collapse = ", "), ".\n",
         "Please check the column names or use the arguments ",
         "(id_col, time_col, status_col, etc.) to map your column names.")
  }

  # ── Extract dimensions ───────────────────────────────────────────────────────
  ids         <- sort(unique(data[[id_col]]))
  iterations  <- sort(unique(data[[iteration_col]]))
  n           <- length(ids)
  m           <- length(iterations)

  id_index    <- stats::setNames(seq_along(ids), as.character(ids))
  iter_index  <- stats::setNames(seq_along(iterations), as.character(iterations))

  # ── Initialise output matrices ───────────────────────────────────────────────
  time_matrix    <- matrix(NA_real_,    nrow = n, ncol = m)
  status_matrix  <- matrix(NA_integer_, nrow = n, ncol = m)
  obsdate_matrix <- matrix(as.Date(NA), nrow = n, ncol = m)

  cov_matrices   <- lapply(covariate_cols, function(cv)
    matrix(NA_real_, nrow = n, ncol = m))
  names(cov_matrices) <- covariate_cols

  credit_start_dates <- rep(as.Date(NA), n)

  # ── Fill matrices row by row ─────────────────────────────────────────────────
  for (row_i in seq_len(nrow(data))) {
    row      <- data[row_i, ]
    i        <- id_index[as.character(row[[id_col]])]
    j        <- iter_index[as.character(row[[iteration_col]])]

    time_matrix[i, j]   <- as.numeric(row[[time_col]])
    status_matrix[i, j] <- as.integer(row[[status_col]])

    if (!is.null(date_col) && date_col %in% names(data)) {
      obsdate_matrix[i, j] <- as.Date(row[[date_col]])
    }

    if (!is.null(start_date_col) && start_date_col %in% names(data)) {
      credit_start_dates[i] <- as.Date(row[[start_date_col]])
    }

    for (cv in covariate_cols) {
      cov_matrices[[cv]][i, j] <- as.numeric(row[[cv]])
    }
  }

  # ── Mark post-default periods as NA (consistent with simulate_mortgage_data)
  for (i in seq_len(n)) {
    event_iter <- which(status_matrix[i, ] == 1)
    if (length(event_iter) > 0) {
      first_ev <- min(event_iter)
      if (first_ev < m) {
        status_matrix[i, (first_ev + 1):m] <- NA_integer_
      }
    }
  }

  # ── Build output list ────────────────────────────────────────────────────────
  result <- list(
    time_matrix        = time_matrix,
    status_matrix      = status_matrix,
    obsdate_matrix     = obsdate_matrix,
    credit_start_dates = credit_start_dates,
    n_individuals      = n,
    n_iterations       = m,
    covariate_names    = covariate_cols
  )

  # Add one named matrix per covariate (x1_matrix, x2_matrix, etc.)
  for (cv in covariate_cols) {
    result[[paste0(cv, "_matrix")]] <- cov_matrices[[cv]]
  }

  message("Data loaded: ", n, " borrowers x ", m, " landmark periods.")
  message("Covariates: ", paste(covariate_cols, collapse = ", "))
  message("Events:     ", sum(status_matrix == 1, na.rm = TRUE),
          " defaults out of ", n, " borrowers (",
          round(mean(rowSums(status_matrix == 1, na.rm = TRUE) > 0) * 100, 1),
          "%)")

  result
}


#' Validate that a MortgageSurv data object is well-formed
#'
#' Checks that a list produced by either \code{\link{simulate_mortgage_data}}
#' or \code{\link{from_real_data}} has the correct structure and dimensions
#' before passing it to modelling functions.
#'
#' @param sim A list as returned by \code{\link{simulate_mortgage_data}} or
#'   \code{\link{from_real_data}}.
#' @param covariate_cols Character vector of expected covariate matrix names
#'   (without the \code{_matrix} suffix). Default \code{c("x1","x2","x3","x4")}.
#'
#' @return Invisible \code{TRUE} if all checks pass. Throws an informative
#'   error if any check fails.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_mortgage_data(n = 200, seed = 42)
#' validate_mortgage_data(sim)   # should pass silently
#' }
#'
#' @export
validate_mortgage_data <- function(sim,
                                    covariate_cols = c("x1", "x2", "x3", "x4")) {

  required <- c("time_matrix", "status_matrix", "obsdate_matrix",
                 "credit_start_dates",
                 paste0(covariate_cols, "_matrix"))

  missing <- setdiff(required, names(sim))
  if (length(missing) > 0)
    stop("Missing elements in data object: ", paste(missing, collapse = ", "))

  n <- nrow(sim$time_matrix)
  m <- ncol(sim$time_matrix)

  for (nm in paste0(covariate_cols, "_matrix")) {
    if (!identical(dim(sim[[nm]]), c(n, m)))
      stop(nm, " has wrong dimensions. Expected (", n, " x ", m, ").")
  }

  if (length(sim$credit_start_dates) != n)
    stop("credit_start_dates must have length ", n, ".")

  message("Validation passed: ", n, " borrowers x ", m, " periods, ",
          length(covariate_cols), " covariates.")
  invisible(TRUE)
}


#' Create a minimal example dataset of real mortgage data
#'
#' Generates a small \code{data.frame} in the format expected by
#' \code{\link{from_real_data}}. Useful for testing and as a template
#' for formatting your own data.
#'
#' @param n Integer. Number of borrowers. Default \code{50}.
#' @param m Integer. Number of landmark periods. Default \code{4}.
#' @param seed Integer. Random seed. Default \code{42}.
#'
#' @return A \code{data.frame} with columns \code{id}, \code{iteration},
#'   \code{time}, \code{status}, \code{x1}, \code{x2}, \code{x3}, \code{x4},
#'   \code{obs_date}, \code{credit_start}.
#'
#' @examples
#' template <- example_real_data(n = 10, m = 3)
#' head(template)
#'
#' # Use it with from_real_data():
#' sim <- from_real_data(template)
#'
#' @importFrom stats runif rnorm
#' @export
example_real_data <- function(n = 50L, m = 4L, seed = 42L) {
  set.seed(seed)

  rows <- vector("list", n * m)
  k    <- 1L

  for (i in seq_len(n)) {
    start_date   <- as.Date("2015-01-01") + sample(0:364, 1)
    default_iter <- if (stats::runif(1) < 0.3) sample(seq_len(m), 1) else NA
    base_time    <- 54 + (i %% 7) * 2 + stats::rnorm(1, 0, 0.5)

    for (j in seq_len(m)) {
      if (!is.na(default_iter) && j > default_iter) {
        # Post-default: skip (NA will be handled by from_real_data)
        next
      }
      obs_date <- start_date + round(j * 30.44 * 2)
      rows[[k]] <- data.frame(
        id           = i,
        iteration    = j,
        time         = round(base_time + j * 0.5, 2),
        status       = as.integer(!is.na(default_iter) && j == default_iter),
        x1           = round(stats::runif(1, 0.01, 0.06), 4),   # interest rate
        x2           = round(stats::runif(1, 0.01, 0.08), 4),   # inflation
        x3           = round(stats::runif(1, 1.5,  5.0),  2),   # LTI
        x4           = round(stats::runif(1, 25,   65),   1),   # age
        obs_date     = obs_date,
        credit_start = start_date
      )
      k <- k + 1L
    }
  }

  do.call(rbind, rows[!sapply(rows, is.null)])
}
