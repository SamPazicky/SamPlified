#' align.sequences
#'
#' Align sequences and export a nice alignment figure
#' 
#' @param file Path to the fasta file.
#' @param substitutionMatrix Argument as in Biostrings::pairwiseAlignment
#' @import msa
#' @import tools
#' @return 
#' @examples 
#' align.sequences("toAlign.fasta")
#' @export

align.sequences <- function(
    file=NULL,
    sequences=NULL,
    substitutionMatrix="BLOSUM50",
    pid.type="PID1",
    paperWidth=7.5,
    shadingMode="similar",
    consensusColors="ColdHot",
    shadingColors="reds",
    showConsensus="none",
    showLegend=TRUE
) {
  
  
  
  if(is.null(file)&is.null(sequences)) {
    stop("Please give path to the fasta file or sequences.")
  } else if(is.null(file)) {
    if(class(sequences)!="AAStringSet") {
      stop("Please give sequences as AAStringSet.")
    }
    su <- paste(names(sequences),collapse="_")
  } else if(is.null(sequences)) {
    if(!str_ends(file,".fasta")) {
      stop("The file must be a fasta file.")
    } else {
      su <- str_remove(file,".fasta$")
    }
    sequences <- readAAStringSet(file)
  } else {
    stop("Please give either sequences or path, not both.")
  }
  
  PAL <- pairwiseAlignment(pattern = sequences[[1]], subject = sequences[[2]], substitutionMatrix=substitutionMatrix)
  pId <- pid(PAL, type=pid.type)
  wPAL <- c(
    paste0(">",names(sequences[1])),
    PAL@pattern %>% as.character(),
    paste0(">",names(sequences[2])),
    PAL@subject %>% as.character()
  )
  # wPALname <- str_split_1(fastafiles,"/") %>% .[length(.)] %>% str_remove(".fasta")
  writeLines(wPAL,paste0(su,"_pairwisealignment.txt"))
  
  
  alignment <- msa(sequences)
  msaPrettyPrint(alignment, output="tex", showNames="none", paperWidth=paperWidth,
                 showLogo="none", askForOverwrite=FALSE, verbose=FALSE,
                 shadingMode=shadingMode,consensusColors=consensusColors, shadingColors=shadingColors,
                 showConsensus=showConsensus,showLegend=showLegend)
  texfile <- readLines("alignment.tex")
  texfile <- sapply(texfile, function(x) 
    str_replace(x,"SAMUEL~1.PAZ","samuel.pazicky")
  ) %>% unname()
  writeLines(texfile,"alignment.tex")
  tools::texi2pdf("alignment.tex", clean=TRUE)
  file.rename("alignment.pdf",paste0(su,"_alignment.pdf"))
  
  saveWidth <- getOption("width")
  options(width=100)
  sink("myAlignment.txt")
  print(alignment, show="complete", halfNrow=-1)
  sink()
  options(width=saveWidth)
  file.rename("myAlignment.txt",paste0(su,"_alignment.txt"))
  
  output <- list(
    alignment=alignment,
    identity=pId
  )
  
  return(output)
  
}
