#' Compute mean/median/SD/SE/CI for a numeric vector
#'
#' @keywords internal
.calc_mean_ci <- function(x,
                          ci_method = c("normal", "bootstrap"),
                          n_boot = 1000,
                          seed = NULL,
                          ...) {
  ci_method <- match.arg(ci_method)
  x <- x[is.finite(x)]
  n <- length(x)

  if (!n) {
    return(c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_))
  }

  m   <- mean(x)
  med <- stats::median(x)
  sd  <- stats::sd(x)
  se  <- sd / sqrt(max(1L, n))

  if (ci_method == "normal") {
    ci <- m + c(-1.96, 1.96) * se
  } else {
    if (!is.null(seed)) set.seed(seed)
    means <- replicate(n_boot, mean(sample(x, size = n, replace = TRUE)))
    ci <- stats::quantile(means, probs = c(0.025, 0.975),
                          names = FALSE, na.rm = TRUE)
  }

  c(m, med, sd, se, ci[1], ci[2])
}
