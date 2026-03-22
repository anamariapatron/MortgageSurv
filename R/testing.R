# ==============================================================================
#  testing.R
#  Pseudo-evaluation metric for the AFT and Landmarking models.
#  Approximates expected event times using median covariate values and
#  compares them against the Weibull-implied baseline.
# ==============================================================================


#' Evaluate landmark risk using median covariate profiles
#'
#' A pseudo-evaluation metric that approximates the predicted event time for a
#' "median borrower" at each landmark period. The procedure:
#' \enumerate{
#'   \item Computes the median of each covariate matrix column (one value per
#'     landmark period).
#'   \item Assembles a (4 × m) percentile matrix and computes the linear
#'     predictor \eqn{\eta_k = \beta_k^{\top} x_{50\%}} for each period.
#'   \item Draws a Weibull random variate scaled by \eqn{\exp(\eta_k)},
#'     yielding an approximate expected event time per period.
#' }
#'
#' @param matrices A list of four numeric covariate matrices (x1, x2, x3, x4),
#'   each of dimension (n × m).
#' @param betas A numeric matrix of regression coefficients (m × 4), one row
#'   per landmark period, one column per covariate.
#' @param thetas A numeric matrix of Weibull parameters (m × 2): column 1 is
#'   the shape and column 2 is the scale.
#' @param probs Numeric. Quantile used to summarise each covariate column.
#'   Default \code{0.5} (median).
#' @param seed Integer. Random seed for the Weibull draws. Default \code{NULL}
#'   (no seed set).
#'
#' @return A list with:
#' \describe{
#'   \item{percentiles_matrix}{4 × m matrix of covariate percentiles.}
#'   \item{dot_products}{m × 1 matrix of linear predictors \eqn{\eta_k}.}
#'   \item{dot_products_exp}{m × 1 matrix of \eqn{\exp(\eta_k)}.}
#'   \item{weibull_vals}{m × 1 matrix of Weibull baseline draws.}
#'   \item{result}{m × 1 matrix of scaled event-time approximations
#'     (\eqn{\text{weibull\_vals} \times \exp(\eta_k)}).}
#' }
#'
#' @details
#' The conclusion from the original study is that event risk increases across
#' landmark periods, and that Landmarking provides limited accuracy while AFT
#' performs even worse on this simulated dataset.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_mortgage_data(n = 500, seed = 42)
#'
#' betas <- matrix(c(
#'    50, 0.9,  0.8,  -0.2,
#'    51, 0.92, 0.81, -0.21,
#'    49, 0.88, 0.79, -0.19,
#'    52, 0.95, 0.82, -0.22,
#'    48, 0.89, 0.78, -0.18,
#'     5, 0.9,  0.8,  -0.2,
#'    51, 0.91, 0.81, -0.21
#' ), nrow = 7, byrow = TRUE)
#'
#' thetas <- matrix(c(
#'   0.5, 22.25000,
#'   0.5, 22.49000,
#'   0.5, 22.49756,
#'   0.5, 22.73315,
#'   0.5, 22.96397,
#'   0.5, 23.19023,
#'   0.5, 23.41217
#' ), nrow = 7, byrow = TRUE)
#'
#' ev <- evaluate_landmark_risk(
#'   matrices = list(sim$x1_matrix, sim$x2_matrix,
#'                   sim$x3_matrix, sim$x4_matrix),
#'   betas    = betas,
#'   thetas   = thetas,
#'   seed     = 42
#' )
#' ev$result
#' }
#'
#' @export
evaluate_landmark_risk <- function(matrices, betas, thetas,
                                    probs = 0.5, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Step 1: percentile per covariate per landmark period
  percentiles_list   <- lapply(matrices, get_percentiles, probs = probs)
  percentiles_matrix <- do.call(rbind, percentiles_list)  # 4 x m

  m <- ncol(percentiles_matrix)
  stopifnot(nrow(betas) == m, nrow(thetas) == m)

  # Step 2: linear predictor for each period
  dot_products <- sapply(seq_len(m), function(i) {
    sum(betas[i, ] * percentiles_matrix[, i])
  })
  dot_products     <- matrix(dot_products, ncol = 1)
  dot_products_exp <- exp(dot_products)

  # Step 3: Weibull baseline draw per period
  weibull_vals <- mapply(function(shape, scale) {
    stats::rweibull(1, shape = shape, scale = scale)
  }, shape = thetas[, 1], scale = thetas[, 2])
  weibull_vals <- matrix(weibull_vals, ncol = 1)

  # Step 4: scaled event-time approximation
  result <- weibull_vals * dot_products_exp

  list(
    percentiles_matrix = percentiles_matrix,
    dot_products       = dot_products,
    dot_products_exp   = dot_products_exp,
    weibull_vals       = weibull_vals,
    result             = result
  )
}
