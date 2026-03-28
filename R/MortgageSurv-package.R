#' MortgageSurv: Survival Analysis Models for Mortgage Default Prediction
#'
#' @description
#' Adapts classical survival analysis models — originally developed in
#' biostatistics — to the context of the UK mortgage market. Mortgage
#' origination is treated as the "birth" event and default as "death", enabling
#' a rigorous framework for studying the dynamics of loan performance in the
#' period preceding default.
#'
#' ## Main workflow
#'
#' 1. **Simulate** a mortgage portfolio with \code{\link{simulate_mortgage_data}}.
#' 2. **Fit** a dynamic Landmarking model with \code{\link{fit_landmark_model}}.
#' 3. **Predict** default probabilities with \code{\link{predict_risk_aft}} /
#'    \code{\link{predict_risk_ph}}.
#' 4. **Evaluate** model performance with \code{\link{evaluate_landmark_risk}}.
#' 5. **Visualise** results with the \code{plot_*()} family of functions.
#'
#' ## Key references
#'
#' - Van Houwelingen, H. C. (2007). Dynamic prediction by landmarking in event
#'   history analysis. *Scandinavian Journal of Statistics*, 34(1), 70–85.
#' - Rizopoulos, D. (2012). *Joint Models for Longitudinal and Time-to-Event
#'   Data*. CRC Press.
#'
#' @keywords internal
"_PACKAGE"
