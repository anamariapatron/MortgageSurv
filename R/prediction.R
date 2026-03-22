# ==============================================================================
#  prediction.R
#  Comparison of predicted default probabilities from AFT and Landmarking models.
# ==============================================================================


#' Extract predicted default probabilities from a PH model
#'
#' For each landmark time-point in \code{landmarks}, computes the probability
#' of default \eqn{P(T \leq t^*)} under a proportional-hazards model fitted with
#' \code{eha::phreg()} (or an equivalent \code{summary()} method accepting
#' \code{type = "survival"} and \code{t}).
#'
#' @param fit A fitted PH model object (e.g., from \code{eha::phreg()}).
#' @param newdata A \code{data.frame} of covariates at which to evaluate the
#'   survival function.
#' @param landmarks Numeric vector of time-points (months).
#'
#' @return A named list of numeric vectors (one per landmark), each containing
#'   the predicted default probability (as a percentage, 0–100) for every row
#'   in \code{newdata}.
#'
#' @examples
#' \dontrun{
#' risk_ph <- predict_risk_ph(fit_weibull_ph, newdata = df_aft,
#'                            landmarks = c(54, 56, 58, 60, 62, 64))
#' }
#'
#' @export
predict_risk_ph <- function(fit, newdata, landmarks) {
  risk_list <- lapply(landmarks, function(t_star) {
    surv_pred <- summary(fit, type = "survival", t = t_star, newdata = newdata)
    sapply(surv_pred, function(x) (1 - x$est) * 100)
  })
  names(risk_list) <- as.character(landmarks)
  risk_list
}


#' Extract predicted default probabilities from an AFT model
#'
#' Analogous to \code{\link{predict_risk_ph}} but for an Accelerated Failure
#' Time model fitted with \code{eha::aftreg()} (or any model with a compatible
#' \code{summary()} method).
#'
#' @param fit A fitted AFT model object.
#' @param newdata A \code{data.frame} of covariates.
#' @param landmarks Numeric vector of evaluation time-points (months).
#'
#' @return A named list of numeric vectors (one per landmark).
#'
#' @examples
#' \dontrun{
#' risk_aft <- predict_risk_aft(fit_weibull_aft, newdata = df_aft,
#'                              landmarks = c(54, 56, 58, 60, 62, 64))
#' }
#'
#' @export
predict_risk_aft <- function(fit, newdata, landmarks) {
  risk_list <- lapply(landmarks, function(t_star) {
    surv_pred <- summary(fit, type = "survival", t = t_star, newdata = newdata)
    sapply(surv_pred, function(x) (1 - x$est) * 100)
  })
  names(risk_list) <- as.character(landmarks)
  risk_list
}


#' Plot predicted-risk density curves from the Landmarking model
#'
#' Produces a 2×3 panel of kernel-density plots, one per landmark, showing the
#' distribution of individual-level default probabilities produced by the
#' dynamic landmarking model.
#'
#' @param landmark_model A fitted landmark model list (from
#'   \code{\link{fit_landmark_model}}).
#' @param landmarks Numeric vector of landmark identifiers matching the names of
#'   \code{landmark_model}. Default \code{c(54, 56, 58, 60, 62, 64)}.
#'
#' @return Invisible \code{NULL}. Draws to the active graphics device.
#'
#' @importFrom graphics par plot polygon lines density
#' @importFrom grDevices adjustcolor
#' @export
plot_landmark_densities <- function(landmark_model,
                                     landmarks = c(54, 56, 58, 60, 62, 64)) {
  old_par <- graphics::par(mfrow = c(2, 3))
  on.exit(graphics::par(old_par))

  for (a in landmarks) {
    preds <- 100 * landmark_model[[as.character(a)]]$data$event_prediction
    d     <- stats::density(preds, na.rm = TRUE)

    graphics::plot(
      d,
      xlab = "Predicted risk of default (%)",
      main = paste("Landmark month", a),
      lwd  = 2,
      col  = "deeppink3",
      ylim = c(0, max(d$y) * 1.2)
    )
    graphics::polygon(d, col = grDevices::adjustcolor("pink", alpha.f = 0.5),
                      border = NA)
    graphics::lines(d, col = "deeppink3", lwd = 2)
  }
  invisible(NULL)
}


#' Plot predicted-risk density curves from the AFT model
#'
#' Produces a 2×3 panel of kernel-density plots, one per landmark, for default
#' probabilities from an AFT model.
#'
#' @param risk_list_aft Named list of predicted risk vectors as returned by
#'   \code{\link{predict_risk_aft}}.
#' @param landmarks Numeric vector of landmark identifiers. Default
#'   \code{c(54, 56, 58, 60, 62, 64)}.
#'
#' @return Invisible \code{NULL}. Draws to the active graphics device.
#'
#' @importFrom graphics par plot lines density
#' @export
plot_aft_densities <- function(risk_list_aft,
                                landmarks = c(54, 56, 58, 60, 62, 64)) {
  old_par <- graphics::par(mfrow = c(2, 3))
  on.exit(graphics::par(old_par))

  for (a in landmarks) {
    d <- stats::density(risk_list_aft[[as.character(a)]], na.rm = TRUE)
    graphics::plot(
      d,
      xlab = "Predicted risk of default (%)",
      main = paste("Landmark month", a),
      lwd  = 2,
      col  = "green3",
      ylim = c(0, max(d$y) * 1.2)
    )
    graphics::lines(d, col = "green3", lwd = 2)
  }
  invisible(NULL)
}


#' Plot overlaid AFT vs. Landmarking density curves
#'
#' Produces a 2×3 panel comparing the kernel-density estimates of predicted
#' default probabilities from the AFT (green) and Landmarking (pink) models,
#' restricted to risks below 4\%.
#'
#' @param landmark_model A fitted landmark model list.
#' @param risk_list_aft Named list of AFT predicted risk vectors.
#' @param landmarks Numeric vector of landmark identifiers. Default
#'   \code{c(54, 56, 58, 60, 62, 64)}.
#' @param output_pdf Character or \code{NULL}. If a file path is provided, the
#'   plot is saved as a PDF at that path. Default \code{NULL} (plots to the
#'   active device).
#'
#' @return Invisible \code{NULL}.
#'
#' @importFrom graphics par plot polygon lines mtext legend
#' @importFrom grDevices adjustcolor pdf dev.off
#' @export
plot_risk_comparison <- function(landmark_model,
                                  risk_list_aft,
                                  landmarks   = c(54, 56, 58, 60, 62, 64),
                                  output_pdf  = NULL) {

  draw_panels <- function() {
    old_par <- graphics::par(
      mfrow   = c(2, 3),
      mar     = c(4, 4, 2, 2),
      oma     = c(6, 0, 0, 0),
      xpd     = TRUE,
      cex.main = 1.5,
      cex.lab  = 1.3,
      cex.axis = 1.2
    )
    on.exit(graphics::par(old_par))

    for (a in landmarks) {
      preds <- 100 * landmark_model[[as.character(a)]]$data$event_prediction
      preds <- preds[preds <= 4]

      d1 <- stats::density(risk_list_aft[[as.character(a)]], na.rm = TRUE,
                           from = 0, to = 4)
      d2 <- stats::density(preds, na.rm = TRUE, from = 0, to = 4)
      ymax <- max(d1$y, d2$y, na.rm = TRUE) * 1.2

      graphics::plot(0, 0, type = "n",
                     xlab = "", ylab = "Density",
                     main = paste("Landmark month", a),
                     xlim = c(0, 4), ylim = c(0, ymax))

      graphics::polygon(c(d1$x, rev(d1$x)), c(rep(0, length(d1$y)), rev(d1$y)),
                        col = grDevices::adjustcolor("green", alpha.f = 0.3),
                        border = "green3")
      graphics::lines(d1, col = "green3", lwd = 2)

      graphics::polygon(c(d2$x, rev(d2$x)), c(rep(0, length(d2$y)), rev(d2$y)),
                        col = grDevices::adjustcolor("pink", alpha.f = 0.3),
                        border = "deeppink3")
      graphics::lines(d2, col = "deeppink3", lwd = 2)
    }

    graphics::mtext("Predicted risk of default (%)", side = 1, outer = TRUE,
                    line = 4, cex = 1.3)
    graphics::par(xpd = NA)
    graphics::legend(
      x      = "bottom",
      legend = c("AFT", "Landmarking"),
      col    = c("green3", "deeppink3"),
      lwd    = 2,
      pt.bg  = c(grDevices::adjustcolor("green", alpha.f = 0.3),
                 grDevices::adjustcolor("pink",  alpha.f = 0.3)),
      pch    = 15,
      bty    = "n",
      cex    = 1.3,
      horiz  = TRUE
    )
  }

  if (!is.null(output_pdf)) {
    grDevices::pdf(output_pdf, width = 12, height = 8)
    draw_panels()
    grDevices::dev.off()
    message("Plot saved to: ", output_pdf)
  } else {
    draw_panels()
  }
  invisible(NULL)
}


#' Compute column-wise percentiles of a covariate matrix
#'
#' A utility function used in model verification: applies a quantile function to
#' each column of a matrix, returning a vector of summary values.
#'
#' @param mat Numeric matrix (rows = individuals, columns = landmark periods).
#' @param probs Numeric vector of probabilities for \code{quantile()}.
#'   Default \code{0.5} (median).
#'
#' @return A numeric matrix with \code{length(probs)} rows and
#'   \code{ncol(mat)} columns.
#'
#' @examples
#' m <- matrix(rnorm(70), nrow = 10)
#' get_percentiles(m, probs = c(0.25, 0.5, 0.75))
#'
#' @export
get_percentiles <- function(mat, probs = 0.5) {
  apply(mat, 2, function(col) stats::quantile(col, probs = probs))
}
