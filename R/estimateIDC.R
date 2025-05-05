#' estimateIDC
#'
#' Estimate IDC of a sample based on a reference
#'
#' @param data A data frame with columns "Accession","Value" and "Sample".
#' @param reference The reference timepoint data frame.
#' @param data.hpi How long is the cycle of the Plasmodium strain in samples?
#' @param reference.hpi How long is the cycle of the Plasmodium strain in the reference?
#' @param scale Logical: Should the data be scaled against reference.hpi?
#' @param log2transform Logical: Should the values be log2 transformed?
#' @param cor.method Corelation method as in argument 'method' of the function cor().
#' @import tidyverse
#' @import gtools
#' @import viridis
#'
#' @return A list containing:
#' \item{data}{Calculated data.}
#' \item{plot}{Correlation plot.}
#'
#' @examples
#' estimateIDC(data)
#' @export


estimateIDC <- function(
  data=NULL,
  reference=samplified::Pf_ref_IDC,
  data.hpi=42,
  reference.hpi=48,
  scale=TRUE,
  log2transform=FALSE,
  cor.method="spearman"
){

  if(is.null(data)) {
    stop("Please include data.")
  } else {
    data.prep <- prepare.data(data,cols=c("Accession","Value","Sample"),no.cols=TRUE)
    if(data.prep$message!="ok") {
      cat("Problem with data parameter:\n")
      cat(data.prep$message)
      if(data.prep$status=="error") {
        stop()
      }
    }
    rm(data.prep)
  }

  ref <- reference %>%
    pivot_longer(cols=starts_with("tp"), values_to="Value",names_to="Ref.point") %>%
    mutate(Ref.point=factor(Ref.point,levels=gtools::mixedsort(unique(Ref.point))))

  if(log2transform) {
    data <- data %>% mutate(Value=log2(Value))
    ref <- ref %>% mutate(Value=log2(Value))
  }

  samples <- gtools::mixedsort(unique(data$Sample))
  reference.points <- levels(ref$Ref.point)

  combs <- expand_grid(samples,reference.points)

  corel<-c()
  for(i in 1:nrow(combs)) {
    rcombs <- combs %>% slice(i)
    rdata <- left_join(
      data %>% dplyr::filter(Sample==rcombs$samples),
      ref %>% dplyr::filter(Ref.point==rcombs$reference.points),
      by="Accession"
    ) %>%
      dplyr::select(starts_with("Value")) %>%
      mutate(across(everything(), ~ ifelse(.x==0,NA,.x)))
    corel[i] <- cor(rdata,use="complete.obs", method=cor.method) %>%
      as.data.frame() %>% pull(1) %>% .[2]
  }
  combs <- combs %>%
    mutate(cor=corel) %>%
    mutate(reference.points=as.numeric(str_remove(reference.points,"tp")))

  reference.levels <- as.numeric(str_remove(reference.points,"tp"))
  if(scale) {
    combs$reference.points <- combs$reference.points*data.hpi/reference.hpi
    reference.levels <- reference.levels*data.hpi/reference.hpi

  }

  corel_plot <- combs %>%
    setNames(c("Sample","Reference","Cor")) %>%
    mutate(Reference=factor(Reference,levels=gtools::mixedsort(unique(reference.levels)))) %>%
    mutate(Sample=factor(Sample,levels=rev(gtools::mixedsort(unique(samples))))) %>%
    ggplot() +
    geom_tile(mapping=aes(x=Reference,y=Sample,fill=Cor)) +
    scale_fill_viridis(name="Correlation", option="turbo") + coord_equal() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.3)) +
    scale_x_discrete(name="Reference timepoint (hpi)",expand=c(0,0)) +
    scale_y_discrete(name="Sample",expand=c(0,0))

  output<-list(
    data=combs,
    plot=corel_plot
  )

  return(output)
}




