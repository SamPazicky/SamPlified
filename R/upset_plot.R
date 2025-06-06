#' upset_plot
#'
#' Plots upset and Euler plots from gene sets.
#'
#' @param sets Data frame with columns 'id' and 'condition'.
#' @param ... Any arguments of UpSetR::upset().
#'
#' @import tidyverse
#' @import UpSetR
#' @return A list with two elements:
#' \item{data}{Enrichment results in table form.}
#' \item{plot}{Summarizing enrichment plot.}
#' $data with the enrichment result in a table form,
#' $plot for the enrichment plot.
#' @examples
#' upset_plot(upsetdata)
#' @export

upset_plot <- function(
    sets=data.frame(),
    ... # any
) {
  require(UpSetR)

  sets.prep <- prepare.data(sets,cols=c("id","condition"),no.cols=TRUE)
  if(sets.prep$message!="ok") {
    cat("Problem with sets parameter:\n")
    cat(sets.prep$message)
    if(sets.prep$status=="error") {
      stop()
    }
  }
  sets <- sets.prep$data
  rm(sets.prep)

  setlist <- sets %>%
    dplyr::select(id,condition) %>%
    unstack()

  upsets <- upset(fromList(setlist), order.by = "freq", nsets=length(setlist), ...)

  return(upsets)

}
