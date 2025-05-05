#' calculateIC50
#'
#' Calculates IC50 from dose-response data
#'
#' @param data A data frame with columns "x" and "y".
#' @import tidyverse
#' @import drc
#'
#' @return A list containing:
#' \item{IC50}{Calculated IC50.}
#' \item{fit}{Fit object.}
#' \item{coefs}{Fitted coefficients.}
#' \item{plot}{Plot with the fitted data.}
#'
#' @examples
#' calcIC50(data)
#' @export

calcIC50 <- function(
  data=NULL
) {

  data.prep <- prepare.data(data,cols=c("x","y"),no.cols=TRUE)
  if(data.prep$message!="ok") {
    cat("Problem with pathways parameter:\n")
    cat(data.prep$message)
    if(data.prep$status=="error") {
      stop()
    }
  }
  data<- data.prep$data
  rm(data.prep)


  fit <- drm(y ~ x, data=data, fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50"),
                                                              fixed = c(NA,NA,NA,NA)))

  coefs <- setNames(
    fit$coefficients,
    c("hill", "min_value", "max_value","ec_50")
  )


  ic_50 <- with(
    as.list(coefs),
    exp(
      log(ec_50) + (1 / hill) * log(max_value / (max_value - 2 * min_value))
    )
  )

  # predictions and confidence intervals.

  mynthroot = function(x,n) {
    (abs(x)^(1/n))*sign(x)
  }
  scientific_10 <- function(x) {
    parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
  }

  topconc <- max(data$x)
  divfactor <- mynthroot(topconc/min(data$x[data$x!=0]),(200-1))
  fakedata <- data.frame("x"=topconc/divfactor^(0:(200-1)) ) %>%
    add_row(x=0) %>%
    mutate(y=predict(fit,.))

  plot <- ggplot() +
    geom_point(data=data, aes(x=x+min(x[x!=0])/10,y=y)) +
    geom_line(data=fakedata, aes(x=x+min(x[x!=0])/10,y=y)) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_bw()

  output <- list(
    IC50=ic_50,
    fit=fit,
    coefs=coefs,
    plot=plot
  )
  return(output)

}
