#' prepare.data
#'
#' Checks whether the columns expected from the data are there and filters on them.
#'
#' @param data Data frame to check and prepare.
#' @param cols Column names expected in the data frame
#' @param no.cols Logical: if not all column names are found, then at least
#' if the number of columns matches, they will just be renamed.
#'
#' @import tidyverse
#' @return A list with three elements:
#' $data is a data frame with the prepared data
#' $message is the message to show to the user
#' $error says whether an error should be returned

prepare.data <- function(
  data=NULL,
  cols=NULL,
  no.cols=FALSE
) {

  df.names <- data %>% as.data.frame() %>% names()
  shared.names <- intersect(df.names,cols)
  error.message <- paste0("The data must contain ",length(cols)," columns named: ", paste(cols,collapse=", "),"\n")

  output.data <- data
  output.message <- "ok"
  output.status <- "error"


  if(length(shared.names)!=length(cols)) {
    if(no.cols) {

      if(length(df.names)!=length(cols)) {
        output.message <- error.message
      } else {
        output.message <- paste0("Unexpected column names, but the right number of columns. Assigning new column names: ", paste(cols,collapse=", "),"\n")
        output.data <- data %>%
          setNames(cols)
        output.status <- "ok"
      }
    } else {
      output.message <- error.message
    }
  } else {
    output.data <- data %>%
      dplyr::select(all_of(cols))
    output.status <- "ok"
  }

  return(list(
    data = output.data,
    message = output.message,
    status = output.status
  ))
}
