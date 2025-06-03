#' Run DIA-NN from R
#'
#' Executes DIA-NN on one or more `.mzML` files using system calls.
#'
#' @param files Character vector of full paths to `.mzML` files.
#' @param output.folder Directory where DIA-NN output (e.g., `report.parquet`) will be saved.
#' @param reannotate Logical. If `TRUE`, enables reannotation of the library.
#' @param unrelated.runs Logical. If `TRUE`, treats runs as unrelated for analysis.
#' @param MBR Logical. If `TRUE`, enables match-between-runs (MBR).
#' @param peptide.length.range Integer vector of length 2 specifying min and max peptide length.
#' @param precursor.charge.range Integer vector of length 2 specifying allowed precursor charge range.
#' @param precursor.mz.range Numeric vector of length 2 specifying m/z range for precursor ions.
#' @param fragment.ion.mz.range Numeric vector of length 2 specifying m/z range for fragment ions.
#' @param precursor.FDR Numeric. FDR threshold for precursor quantification (default 0.01).
#' @param diann.path Path to the `diann.exe` binary (default `"diann.exe"` assumes it's in the system PATH).
#' @param lib Path to the spectral library `.tsv` file.
#' @param fasta Path to the FASTA file used for identification and inference.
#' @param threads Integer. Number of CPU threads to use (default is 8).
#'
#' @return Returns the system exit code from the DIA-NN command.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' run.diann(
#'   files = list.files("DIA_data", pattern = "\\.mzML$", full.names = TRUE),
#'   output.folder = "DIA_results",
#'   lib = "library.tsv",
#'   fasta = "fastafile.fasta",
#'   threads = 8
#' )
#' }

run.diann <- function(
    files,
    output.folder,
    reannotate = TRUE,
    unrelated.runs = FALSE,
    MBR = TRUE,
    peptide.length.range = c(5, 30),
    precursor.charge.range = c(1, 4),
    precursor.mz.range = c(300, 1800),
    fragment.ion.mz.range = c(200, 1800),
    precursor.FDR = 0.01,
    diann.path = "diann.exe",
    lib = "D:/DIA_helpfiles/library.tsv",
    fasta = "D:/DIA_helpfiles/Pf3D7_v63_HsUniprot241009_iRT_contam_decoy_formatted.fasta",
    threads = 8
) {
  # Ensure files are in a character vector
  if (!is.character(files)) stop("`files` must be a character vector of mzML paths.")

  # Escape and quote file paths
  file.args <- paste(paste0("--f ", shQuote(files)), collapse = " ")

  # Output file path
  dir.create(output.folder,showWarnings=FALSE)
  out.path <- file.path(output.folder, "report.parquet")

  # Basic arguments
  args <- c(
    file.args,
    paste("--lib", shQuote(lib)),
    paste("--threads", threads),
    "--verbose 1",
    paste("--out", shQuote(out.path)),
    paste("--qvalue", precursor.FDR),
    "--matrices",
    "--xic",
    "--peptidoforms",
    "--met-excision",
    "--missed-cleavages 2",
    "--unimod4",
    "--var-mods 2",
    "--var-mod UniMod:35,15.994915,M",
    "--var-mod UniMod:1,42.010565,*n",
    "--individual-mass-acc",
    "--individual-windows",
    "--rt-profiling",
    paste("--fasta", shQuote(fasta)),
    paste("--min-pep-len", peptide.length.range[1]),
    paste("--max-pep-len", peptide.length.range[2]),
    paste("--min-pr-charge", precursor.charge.range[1]),
    paste("--max-pr-charge", precursor.charge.range[2]),
    paste("--min-pr-mz", precursor.mz.range[1]),
    paste("--max-pr-mz", precursor.mz.range[2]),
    paste("--min-fr-mz", fragment.ion.mz.range[1]),
    paste("--max-fr-mz", fragment.ion.mz.range[2])
  )

  # Optional flags
  if (reannotate) args <- c(args, "--reannotate")
  if (unrelated.runs) args <- c(args, "--unrelated-runs")
  if (MBR) args <- c(args, "--reanalyse")

  # Construct command
  cmd <- paste(shQuote(diann.path), paste(args, collapse = " "))

  # Run command
  cat("Running:\n", cmd, "\n")
  system(cmd)
}
