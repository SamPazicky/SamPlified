#' GSEA
#'
#' Gene set enrichment analysis
#'
#' @param pathways Data frame with columns 'id','pathway.type' and 'pathway.name', or at least a
#' data frame with three columns and this content.
#' @param labels.col ranks Data frame with columns 'id' and 'rank', or at least a data frame with
#' two columns and this content
#' @param types Vector: what types of enrichment to do? Choose from pathway.type of pathways or simply "all".
#' @param minsize minimal size of gene set to test.
#' @param maxsize maximal size of gene set to test.
#' @param p.cutoff Numeric: p-value cutoff. Default is 0.05.
#' @param scoreType Same as in fgseaSimple.
#' @import tidyverse
#' @import fgsea
#' @return A list with four elements:
#' \item{enrichment.table}{Data with the full GSEA results.}
#' \item{enriched.table}{Data with the enriched pathways}
#' \item{plot}{GSEA distribution plots for each enriched term.}
#' \item{enrichment.plot}{Summarizing enrichment plot.}
#'
#' @examples
#' gsea.result <- GSEA(pathways,ranks)
#' @export

GSEA <- function (
  pathways=data.frame(), # data frame with three columns in this order: ids, pathway type, pathway name
  ranks=data.frame(), # data frame with names "id" and "rank"
  types="all",
  minsize=5,
  maxsize=200,
  p.cutoff=0.05,
  scoreType="std"
){

  pathway.prep <- prepare.data(pathways,cols=c("id","pathway.type","pathway.name"),no.cols=TRUE)
  if(pathway.prep$message!="ok") {
    cat("Problem with pathways parameter:\n")
    cat(pathway.prep$message)
    if(pathway.prep$status=="error") {
      stop()
    }
  }
  pathways <- pathway.prep$data
  rm(pathway.prep)

  ranks.prep <- prepare.data(ranks,cols=c("id","rank"), no.cols=TRUE)
  if(ranks.prep$message!="ok") {
    cat("Problem with pathways parameter:\n")
    cat(ranks.prep$message)
    if(ranks.prep$status=="error") {
      stop()
    }
  }
  ranks <- ranks.prep$data
  rm(ranks.prep)

  # convert to required format (named vector)
  ranking <- setNames(ranks$rank, ranks$id)

  # prepare pathway list from annotation table
  pathwaylist <- pathways %>%
    split(.$pathway.type) %>%
    lapply(function(x) x %>% dplyr::select(!pathway.type) %>% unstack())

  # calculate GSEA
  pathtouse <- if(types[1]=="all") names(pathwaylist) else types

  GSEA_result <- list()
  for(p in pathtouse) {
    cat(paste0("Calculating gene set enrichment for pathways in ",p,"..."))
    GSEA_result[[p]] <- fgsea(pathwaylist[[p]], ranking, minSize=minsize, maxSize=maxsize, gseaParam = 1, nproc=1, scoreType=scoreType)
    cat(paste0("\rCalculating gene set enrichment for pathways in ",p,"...done.\n"))
  }
  # add nice pathway names (remove the numbers etc...)
  GSEA_result <- lapply(GSEA_result, function(x)
    x %>% mutate(pathway.name=str_remove(pathway,"GO-[[:digit:]]+-[[:alpha:]]+-|KEGG-pfa[[:digit:]]+-|Sach-|Till-|MPM-MPMP[[:digit:]]+-"))
  )
  # filter by adjusted p.value and put together
  GSEA_enriched <- lapply(GSEA_result, function(x) filter(x,padj<=p.cutoff))
  GSEA_enriched_table <- reduce(GSEA_enriched,bind_rows) %>% arrange(padj)

  if(nrow(GSEA_enriched_table)==0) {
    return(list(enrichment.table=NULL,enriched.table=NULL,plot=NULL,enrichment.plot=NULL))
  }
  # collapse
  collapsedPathways <- collapsePathways(GSEA_enriched_table, reduce(pathwaylist,append), ranking)

  # Convert mainPathways to a data frame and add a column with daughter pathways,
  # then join the entire GSEA table with all statistics
  mainPathways <- data.frame(collapsedPathways$parentPathways) %>%
    setNames("pathway") %>% rownames_to_column("daughter") %>%
    mutate(pathway=ifelse(is.na(pathway),daughter,pathway)) %>%
    group_by(pathway) %>% summarise(corresponding_pathways=paste(daughter,collapse=";;")) %>% ungroup() %>%
    mutate(corresponding_pathways=paste0(pathway,";;",corresponding_pathways)) %>%
    left_join(GSEA_enriched_table)

  # plots
  named_pathwaylist <- reduce(pathwaylist,append)[mainPathways$pathway] %>% setNames(mainPathways$pathway.name)
  plotPathways <- mainPathways %>% dplyr::select(!pathway) %>% mutate(pathway=pathway.name)
  plot <- plotGseaTable(named_pathwaylist, ranking, plotPathways,
                gseaParam=0.5)

  enrichment.plot <- mainPathways %>%
    mutate(Enrichment=ifelse(NES>0,"Over","Under")) %>%
    ggplot(aes(reorder(pathway.name, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score") +
    theme_bw()


  result <- list(enrichment.table=GSEA_result,enriched.table=mainPathways,plot=plot,enrichment.plot=enrichment.plot)

  return(result)
}
