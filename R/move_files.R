#' move_files
#'
#' Move files into a new destination, with exceptions. Never moves folders.
#'
#' @param from Relative or absolute path from which to copy.
#' @param to Relative or absolute path to which to copy, can be a new folder as well.
#' @param except Vector of names of files that should not be copied.
#' @return Nothing. The function just moves files.
#' @examples
#' move_files(".","newfolder",except="mycode.R")
#' @export

move_files <- function(from, to, except = character()) {
  # Expand paths
  from <- normalizePath(from, mustWork = TRUE)
  to <- normalizePath(to, mustWork = FALSE)

  # Create the destination folder if it doesn't exist
  if (!dir.exists(to)) {
    dir.create(to, recursive = TRUE)
  }

  # List all items in the source directory
  all_items <- list.files(from, full.names = TRUE)

  # Filter: include only files, exclude folders and specified exceptions
  file_info <- file.info(all_items)
  files_to_move <- all_items[
    !file_info$isdir &
      !(basename(all_items) %in% except)
  ]

  # Move the files
  moved <- file.rename(files_to_move, file.path(to, basename(files_to_move)))

  # Return a named logical vector showing which files were moved
  names(moved) <- basename(files_to_move)
  return(moved)
}
