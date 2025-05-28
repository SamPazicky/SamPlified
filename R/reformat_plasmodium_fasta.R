#' Reformat Plasmodium FASTA Headers
#'
#' This function reads a FASTA file with Plasmodium protein headers and reformats
#' the headers into a standardized UniProt-like format. It handles both normal and
#' reversed entries (with ">rev_" prefix).
#'
#' @param input Character. Path to the input FASTA file to be reformatted.
#' @param output Character. Path where the reformatted FASTA file will be saved.
#'
#' @importFrom stringr str_extract
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Reformat a Plasmodium FASTA file
#' reformat_plasmodium_fasta("raw_plasmodium.fasta", "formatted_plasmodium.fasta")
#' }
reformat_plasmodium_fasta <- function(input, output) {
  lines <- readLines(input)
  output_lines <- character(length(lines))

  pb <- txtProgressBar(min = 0, max = length(lines), style = 3)

  for (i in seq_along(lines)) {
    line <- lines[i]

    if (startsWith(line, ">PF3D7") || startsWith(line, ">rev_PF3D7")) {
      is_rev <- startsWith(line, ">rev_")

      # Remove '>' and optionally 'rev_' from the full header
      raw_line <- sub("^>", "", line)
      clean_line <- if (is_rev) sub("^rev_", "", raw_line) else raw_line

      # Extract ID (without rev_)
      id <- stringr::str_extract(clean_line, "^[^ ]+")

      # Extract gene name
      gene <- stringr::str_extract(clean_line, "gene=PF3D7[^ |]+")
      gene <- sub("gene=", "", gene)

      # Extract description
      description <- stringr::str_extract(clean_line, "gene_product=[^|]+")
      description <- sub("gene_product=", "", description)

      # Fallbacks
      if (is.na(gene)) gene <- id
      if (is.na(description)) description <- "Plasmodium protein"

      # Remove rev_ from ID field in rev_ entries
      id_clean <- sub("^rev_", "", id)

      # Create new header
      prefix <- if (is_rev) "rev_sp" else "sp"
      new_header <- sprintf(">%s|%s|%s %s OS=Plasmodium falciparum OX=5833 GN=%s PE=1 SV=1",
                            prefix, id_clean, id_clean, description, gene)

      output_lines[i] <- new_header
    } else {
      output_lines[i] <- line
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)
  writeLines(output_lines, output)
  message("✅ Reformatted FASTA written to: ", output)
}
