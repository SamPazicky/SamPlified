#' Perform PCA on a protein count matrix using DESeq2 and return labeled plot
#'
#' This function takes a matrix of protein counts (proteins in rows, samples in columns),
#' automatically extracts sample condition from sample names, performs variance-stabilizing transformation (VST),
#' runs PCA using DESeq2, and returns a `ggplot2` PCA plot with sample labels.
#'
#' @param protein_matrix A numeric matrix of protein abundances (rows = proteins, columns = samples).
#'                       Column names must be sample names containing the condition as the first token, separated by `_`.
#' @param ntop Integer. Number of top genes (by variance) to use for PCA. Default is 10.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{data}{A data frame with PCA coordinates and metadata.}
#'   \item{plot}{A `ggplot` object of the PCA with sample labels.}
#' }
#'
#' @import DESeq2
#' @import tidyverse
#' @importFrom ggrepel geom_text_repel
#'
#' @export
#'
#'
DESeq2.pca.analysis <- function(
    protein_matrix = NULL,
    ntop=10
) {

  # Load required package
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Please install the 'DESeq2' package first.")
  }
  library(DESeq2)

  # Match transformation method
  transformation <- match.arg(transformation)

  # Check row/column names
  if (is.null(rownames(protein_matrix)) || is.null(colnames(protein_matrix))) {
    stop("Protein matrix must have row names (proteins) and column names (sample names).")
  } else {
    protein_matrix[is.na(protein_matrix)] <-  min(protein_matrix[protein_matrix > 0], na.rm = TRUE) / 2
    protein_matrix <- apply(protein_matrix, 2, as.numeric)
  }
  overtop <- ceiling(max(protein_matrix)/.Machine$integer.max)
  if(overtop>1) {
    protein_matrix <- protein_matrix/overtop
  }

  # Reorder sample metadata to match column order in matrix
  sample_metadata <- data.frame(
    samples = colnames(protein_matrix)
  ) %>% mutate(condition=str_split_i(samples,"_",1)) %>%
    column_to_rownames("samples")

  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = round(protein_matrix),
    colData = sample_metadata,
    design = as.formula("~ condition")
  )

  dds_trans <- vst(dds, blind = TRUE)

  # PCA plot
  pcadata <- DESeq2::plotPCA(dds_trans, intgroup = "condition",ntop=ntop, returnData = TRUE)
  pcaplot <- pcadata  %>%
    as.data.frame() %>%
    ggplot(aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3) +
    labs(title = "PCA with sample labels",
         x = paste0("PC1: ", round(attr(pcadata, "percentVar")[1] * 100), "% variance"),
         y = paste0("PC2: ", round(attr(pcadata, "percentVar")[2] * 100), "% variance")) +
    theme_bw() +
    theme(legend.position="none")

  output <- list(
    data=pcadata,
    plot=pcaplot
  )
  return(output)
}
