#' Convert Thermo .raw files to .mzML using msconvert
#'
#' This function uses `msconvert` (from ProteoWizard) to convert and merge Thermo `.raw` files
#' grouped by sample name into `.mzML` files. Filenames must follow the pattern
#' `DIA4x[1-4]_SampleName.raw`.
#'
#' @param raw_dir Path to the directory containing `.raw` files.
#' @param samples Optional character vector of sample names to convert. If NULL (default),
#'   all detected samples will be processed.
#' @param verbose Logical; if TRUE, prints each `msconvert` command being run.
#'
#' @return Invisibly returns a named list of the output mzML file paths by sample.
#' @export
#'
#' @examples
#' \dontrun{
#' convert_raw_to_mzML("D:/User MS/Sam/pulseDIA_PISAtrial5_CtrlvsDSM")
#' }
#'

run.msconvert <- function(raw_dir, samples = NULL, verbose = TRUE) {
  # Ensure dependencies
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required.")

  # Check msconvert availability
  if (Sys.which("msconvert") == "") stop("msconvert not found in PATH. Please install ProteoWizard.")

  # List .raw files
  raw_files <- list.files(raw_dir, pattern = "\\.raw$", full.names = TRUE)
  if (length(raw_files) == 0) {
    stop("No .raw files found in the directory: ", raw_dir)
  }

  # Extract sample names
  sample_names <- raw_files %>%
    basename() %>%
    stringr::str_extract("DIA4x[1-4]_(.+)\\.raw") %>%
    stringr::str_remove("^DIA4x[1-4]_") %>%
    stringr::str_remove("\\.raw$")

  # Group files by sample
  file_groups <- split(raw_files, sample_names)

  # Filter samples if requested
  if (!is.null(samples)) {
    file_groups <- file_groups[names(file_groups) %in% samples]
    if (length(file_groups) == 0) stop("No matching samples found.")
  }

  output_files <- list()

  # Convert files
  for (sample in names(file_groups)) {
    files <- file_groups[[sample]]
    outname <- paste0(sample, ".mzML")

    cmd <- paste(
      "msconvert",
      paste(shQuote(files), collapse = " "),
      "--merge",
      "--filter", "\"peakPicking vendor msLevel=1-\"",
      "--filter", "\"zeroSamples removeExtra 1-\"",
      "--mzML",
      "-v",
      "--outdir", shQuote(raw_dir),
      "--outfile", shQuote(outname)
    )

    if (verbose) {
      cat("Running:", cmd, "\n")
    }

    system(cmd)
    output_files[[sample]] <- file.path(raw_dir, outname)
  }

  invisible(output_files)
}
