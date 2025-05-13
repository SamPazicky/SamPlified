#' fit_cos
#'
#' Fits cosine function into data.
#'
#' @param data A data frame with columns for x and y. The data must be z-scaled!
#' @param x_col Unquoted name of the x column.
#' @param y_col Unquoted name of the y column.
#' @param period Numeric value of period, e.g. 42 for 42h-long IDC of P. falciparum.
#' @import minpack.lm
#' @import tidyverse
#' @return A list containing:
#' \item{fit}{Fit object}
#' \item{coef}{Fitted coefficients.}
#' \item{R2}{R2 goodness of fit value.}
#' @examples
#' fit_cos(periodic_data,condition,z.log2FPKM,42)
#' @export
fit_cos <- function(data, x_col, y_col, period = 5) {

  require(minpack.lm)
  # Convert column names from arguments to actual data
  x <- data[[deparse(substitute(x_col))]]
  y <- data[[deparse(substitute(y_col))]]
  curfitdata <- data.frame(x = x, y = y)

  # First guess: position of the max y
  C_guess <- curfitdata %>% slice_max(y, n = 1) %>% pull(x)
  C_guess <- ifelse(C_guess == 0, 0, 2 * pi * (C_guess / period))

  # Try first fit
  curfit <- try(nlsLM(y ~ A * cos((2 * pi / period) * x - C),
                      data = curfitdata,
                      start = list(A = 1, C = C_guess),
                      lower = c(A = 0, C = 0),
                      upper = c(A = Inf, C = 2 * pi),
                      control = nls.lm.control(maxiter = 500)),
                silent = TRUE)

  # Retry with second highest point if first fails
  if (class(curfit) != "nls") {
    C_guess2 <- curfitdata %>% arrange(desc(y)) %>% slice(2) %>% pull(x)
    C_guess2 <- 2 * pi * (C_guess2 / period)

    curfit <- try(nlsLM(y ~ A * cos((2 * pi / period) * x - C),
                        data = curfitdata,
                        start = list(A = 1, C = C_guess2),
                        lower = c(A = 0, C = 0),
                        upper = c(A = Inf, C = 2 * pi),
                        control = nls.lm.control(maxiter = 500)),
                  silent = TRUE)
  }

  # Output results
  if (class(curfit) == "nls") {
    observed <- curfitdata$y
    predicted <- predict(curfit, newdata = curfitdata)
    RSS <- sum((observed - predicted)^2)
    TSS <- sum((observed - mean(observed))^2)
    R2 <- 1 - (RSS / TSS)

    return(list(fit = curfit,
                coef = coef(curfit),
                R2 = R2))
  } else {
    return(list(fit = NULL,
                coef = c(A = NA, C = NA),
                R2 = NA))
  }
}
