#' GEA
#'
#' Gene  enrichment analysis
#'
#' @param pathways Data frame with columns 'id','pathway.type' and 'pathway.name', or at least a
#' data frame with three columns and this content.
#' @param genes Vector with gene IDs, must have overlap with IDs in the pathways data frame.
#' @param minsize minimal size of gene set to test.
#' @param p.cutoff Numeric: p-value cutoff. Default is 0.05.
#' @param plot.padj.cutoff Numeric: Adjusted p-value cutoff for plotting. Default is 0.05.
#' @param plot.stat Character string: which p-value should be used for plotting? "HP" for hypergeometric
#' p-value and "BP" for binomial.
#'
#' @import tidyverse
#' @return A list with two elements:
#' \item{data}{Data with the enrichment result in a table form}
#' \item{plot}{Enrichment plot}
#'
#' @examples
#' gea.result <- GEA(pathways,genes)
#' @export

GEA <- function(
  pathways=data.frame(), # data frame with three columns in this order: ids, pathway type, pathway name
  genes=NULL,
  minsize=5,
  p.cutoff=0.05,
  p.adj.method="BH",
  plot.padj.cutoff=0.05,
  plot.stat="HP") #HP or BP
{

  if(!plot.stat %in% c("HP","BP")) {
    stop("plot.stat can only be 'HP' or 'BP'.")
  }
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

  if(is.null(genes)) {
    stop("Please include genes.")
  }
  pint <- length(intersect(genes,pathways$id))
  if(pint<minsize) {
    stop(paste0("Insufficient number of genes (",pint,")."))
  }

  r <- r.adj <- data.frame(matrix(ncol = 10, nrow = 0)) %>%
    setNames(c("P","Db","T","D","W","S","UD","HP","BP","GL")) %>% rbind.data.frame()
  # r=data.frame(P="PathwayName",Db="Database", T="Annotated",D="Drawn",W="Anno_inPathway",S="InputPathway",
  #              UD="Direction of deviation from expected", HP="Hypergeometric_P",BP="Binomial_P",GL="GeneList")

  dbs <- pathways$pathway.type %>% unique()

  for(db in dbs)  {
    dbg = pathways %>% dplyr::filter(pathway.type==db) %>% pull(id) %>% unique()
    dbGN <- length(dbg)
    db_inputGN <- length(intersect(dbg,genes))
    groups <- pathways %>% dplyr::filter(pathway.type==db) %>% pull(pathway.name) %>% unique()

    for(gr in groups) {
      gr_genes = pathways %>% dplyr::filter(pathway.name==gr) %>% pull(id) %>% unique()
      grGN=length(gr_genes)
      gr_inputGN = length(intersect(gr_genes,genes))
      #print(gr_inputGN)
      if(gr_inputGN >= minsize) {

        hp=1-phyper(gr_inputGN,grGN,(dbGN-grGN),db_inputGN);
        bp=binom.test(gr_inputGN,db_inputGN,p=grGN/dbGN,alternative="greater")$p.value;

        if (hp<p.cutoff) {
          rr=data.frame(P=gr,Db=db, T=as.character(dbGN),D=as.character(db_inputGN),W=as.character(grGN),S=as.character(gr_inputGN),
                        UD="over", HP=as.character(hp),BP=as.character(bp),GL=paste(intersect(gr_genes,genes),collapse=","))
          r=rbind(r,rr)
        }
        if (hp>(1-p.cutoff)) {
          rr=data.frame(P=gr,Db=db,T=as.character(dbGN),D=as.character(db_inputGN),W=as.character(grGN),S=as.character(gr_inputGN),
                        UD="under", HP=as.character(hp),BP=as.character(bp),GL=paste(intersect(gr_genes,genes),collapse=","))
          r=rbind(r,rr)
        }
      }
    }
  }

  # p.adjustment
  r.adj=data.frame(G="Group",P="PathwayName",Db="Database", T="Annotated",D="Drawn",W="Anno_inPathway",S="InputPathway",
                   UD="Direction of deviation from expected", HP="Hypergeometric_P",BP="Binomial_P",
                   GL="GeneList", HP.adj="Hypergeometric_P.adj", BP.adj="Binomial_P.adj")

  r.adj <- r %>%
    mutate(across(c(HP,BP), as.numeric)) %>%
    group_by(Db) %>%
    mutate(HP.adj=p.adjust(HP,p.adj.method,n())) %>%
    mutate(BP.adj=p.adjust(BP,p.adj.method,n())) %>%
    ungroup() %>%
    setNames(c("pathway.name","pathway.type","annotated","drawn","anno_inPathway","inputPathway",
               "direction","HP","BP","GL","HP.adj","BP.adj"))

  r.plotdata <- r.adj %>%
    dplyr::select(pathway.name,inputPathway,direction,starts_with(plot.stat)) %>%
    setNames(c("Pathway","Size","Enrichment","p","p.adj")) %>%
    mutate(Size=as.numeric(Size)) %>%
    dplyr::filter(p.adj<=plot.padj.cutoff) %>%
    mutate(p.adj=ifelse(p.adj==0,min(p.adj)*0.01,p.adj)) %>%
    mutate(Pathway=str_remove(Pathway,"GO-[[:digit:]]+-[[:alpha:]]+-|KEGG-pfa[[:digit:]]+-|Sach-|Till-|MPM-MPMP[[:digit:]]+-")) %>%
    ggplot(aes(reorder(Pathway,p.adj), -log10(p.adj))) +
    geom_point(aes(fill = Enrichment, size = Size), shape=21) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = -log10(0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Adjusted p-value") +
    theme_bw()

  output <- list(
    data=r.adj,
    plot=r.plotdata
  )

  return(output)
}
