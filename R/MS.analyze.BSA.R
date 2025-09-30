#' Analyze BSA MS File and Generate PDF Report
#'
#' Reads a raw mass spectrometry file, computes the Total Ion Chromatogram (TIC),
#' Extracted Ion Chromatograms (XICs) for specified m/z values, calculates half peak widths (FWHM),
#' generates a PDF report with plots and a table, and returns the processed data and plots.
#'
#' @param file Character. Path to the raw MS file. Required.
#' @param mz Numeric vector. m/z values for XIC extraction. Default is c(653.36, 722.32).
#' @param outpdf Character. Name of the output PDF file. Default is "MS_report.pdf".
#'   The file will be written to the current working directory.
#'
#' @return A list containing:
#' \item{tic}{Tibble with TIC data (time, intensity).}
#' \item{tic.plot}{ggplot object of the Total Ion Chromatogram.}
#' \item{halflifes}{Named vector of FWHM values for each XIC.}
#' \item{xic.plots}{List of ggplot objects for each extracted ion chromatogram.}
#'
#' @examples
#' \dontrun{
#' raw_file <- "path/to/BSA.raw"
#' MS.analyze.BSA(raw_file)
#' }
#' @import rawrr
#' @import gsignal
#' @import rmarkdown
#' @import knitr
#' @import tidyverse
#' @export
MS.analyze.BSA <- function(
    file = NULL,
    mz = c(653.36, 722.32),
    outpdf = "MS_report.pdf"
) {

  # --- Load libraries ---
  library(rawrr)
  library(tidyverse)
  library(gsignal)
  library(rmarkdown)
  library(knitr)

  if (is.null(file)) stop("Please provide a file path to a raw MS file.")

  # --- TIC ---
  tic <- readChromatogram(file, type = "tic")
  tictable <- tibble(time = tic$times, intensity = tic$intensities) %>%
    mutate(across(everything(), as.numeric))
  ticmax <- tictable %>% slice_max(intensity, n = 1)
  tic_time <- round(ticmax$time, 3)
  tic_intensity <- format(ticmax$intensity, scientific = TRUE, digits = 3)
  ticplot <- ggplot(tictable, aes(time, intensity)) +
    geom_line() +
    labs(title="Total Ion Chromatogram", x="Retention time (min)", y="Intensity") +
    theme_bw()

  # --- XICs ---
  XICs <- readChromatogram(file, mass = mz, tol = 10, type = "xic")
  XICs <- lapply(XICs, function(x) tibble(time = x$times, intensity = x$intensities))
  names(XICs) <- paste0("mz ", mz)

  halflifes <- sapply(XICs, function(df) gsignal::fwhm(df$time, df$intensity), simplify = TRUE)
  halflifes_df <- data.frame(
    mz = names(halflifes)%>%str_remove("mz"),
    FWHM = paste0(round(as.numeric(halflifes)*60,1)," s")
  )

  plots <- list()
  for (nm in names(XICs)) {
    plots[[nm]] <- ggplot(XICs[[nm]], aes(time, intensity)) +
      geom_line() +
      ggtitle(nm) +
      labs(x="Retention time (min)", y="Intensity") +
      theme_bw()
  }


  # --- Create temporary Rmd file ---
  rmd_file <- tempfile(fileext = ".Rmd")

  cat(file = rmd_file, '
---
title: "BSA analysis report"
output: pdf_document
---

## File Information

File: `r basename(file)`
Date: `r format(Sys.time(), "%Y-%m-%d %H:%M:%S")`

## Total Ion Chromatogram

```{r tic-plot, echo=FALSE, fig.width=7, fig.height=4}
ticplot
```

### Most Intense Peak

Retention time: `r tic_time`
Intensity: `r tic_intensity`

## Extracted Ion Chromatograms (XICs)

```{r xic-plots, echo=FALSE, fig.width=7, fig.height=4, warning=FALSE, message=FALSE}
for (p in plots) {
  print(p)
}
```

## Half Peak Widths (FWHM)

```{r table, echo=FALSE}
kable(halflifes_df, caption="Half Peak Widths (FWHM)")
```
')

  # --- Render PDF ---
  outpdf <- normalizePath(outpdf, mustWork = FALSE)
  rmarkdown::render(rmd_file, output_file = outpdf, quiet = TRUE)
  message("Report written to: ", outpdf)

  output <- list(
    tic <- tictable,
    tic.plot <- ticplot,
    halflifes <- halflifes,
    xic.plots <- plots
  )
  return(output)
}
