#' B.statistics
#'
#' Calculate Bayesian statistics.
#'
#' @param data A data frame with names id, condition, replicate and value.
#' @param control Name of the control condition (as in the data).
#' @param fdr_hit_threshold FDR threshold for hits.
#' @param fdr_candidate_threshold FDR threshold for candidates.
#' @param fc_hit_threshold FC threshold for hits.
#' @param fc_candidate_threshold FC threshold for candidates.
#' @import tidyverse
#' @import limma
#' @import Biobase
#' @return A list containing:
#' \item{data}{Resulting data table.}
#' \item{hitcounts}{Counts of hits.}
#' \item{limma_vs_FDR}{Comparison of limma vs FDR p-value estimations.}
#' \item{volcano}{Volcano plot.}
#' \item{volcano.separate}{Separate volcano plots for each comparison.}
#' \item{MA}{MA plot.}
#' @examples
#' B.statistics(data,"control")
#' @export

B.statistics <- function(
  data=NULL, # data frame with names id, condition, replicate and value.
  control=NA,
  fdr_hit_threshold=0.05,
  fdr_candidate_threshold=0.25,
  fc_hit_threshold=1,
  fc_candidate_threshold=0.5
) {

  if(is.null(data)) {
    stop("Please include data.")
  } else {
    data.prep <- prepare.data(data,cols=c("id","condition","replicate","value"),no.cols=TRUE)
    if(data.prep$message!="ok") {
      cat("Problem with pathways parameter:\n")
      cat(data.prep$message)
      if(data.prep$status=="error") {
        stop()
      }
    }
    rm(data.prep)
  }

  ndata <- data %>%
    mutate(condition=make.names(condition))

  edata <- ndata %>%
    unite(col="sample",condition,replicate,sep=";;") %>%
    pivot_wider(id_cols=id,names_from=sample,values_from=value) %>%
    column_to_rownames("id") %>%
    as.matrix()

  pdata <- ndata %>%
    dplyr::select(!c(id,value)) %>%
    distinct(.keep_all=TRUE) %>%
    unite(col="sample",condition,replicate,sep=";;",remove=FALSE) %>%
    setNames(c("ID","sample","rep"))
  rownames(pdata) <- pdata$ID

  fdata <- ndata %>%
    unite(col="sample",condition,replicate,sep=";;",remove=FALSE) %>%
    pivot_wider(id_cols=id,names_from=sample,values_from=value) %>%
    as.data.frame()
  rownames(fdata) <- fdata$id

  ES_widedata <- ExpressionSet(assayData = edata,
                               phenoData = AnnotatedDataFrame(pdata),
                               featureData = AnnotatedDataFrame(fdata))

  if(!validObject(ES_widedata)) {
    stop("Data is invalid")
  }

  # make contrasts
  Xcontrol <- make.names(control)
  samples <- unique(ndata$condition)
  samples <- samples[-which(samples==Xcontrol)]
  comparisons <- c()
  for(sample in samples) {
    comparisons[sample] <-paste0(sample,"-",Xcontrol)
  }

  limma.results <- NULL

  # comparison <- paste0(sample,"-",control)
  limma.cond <- factor(pdata$sample)
  limma.rep <- factor(pdata$rep)

  if(length(levels(limma.rep))>1) {
    contrast.matrix <- model.matrix( ~ 0 + limma.cond + limma.rep)
  } else {
    contrast.matrix <- model.matrix( ~  0 + limma.cond)

  }
  colnames(contrast.matrix) <- gsub("limma.cond", "", colnames(contrast.matrix))

  # calculate statistics
  limma.object <- limma::eBayes(
    contrasts.fit(
      lmFit(ES_widedata, design = contrast.matrix, weights = NULL), #limma_weights),
      makeContrasts(contrasts = comparisons, levels = contrast.matrix)
    )
  )

  limma_results_all <- list()
  for(comp in seq_along(comparisons)) {
    limma_results_i <- topTable(limma.object, sort.by="t",number=Inf, coef=comparisons[comp])
    names(limma_results_i)[grep("P.Value", names(limma_results_i))] <- "pvalue.limma"
    names(limma_results_i)[grep("adj.P.Val", names(limma_results_i))] <- "fdr.limma"
    limma_results_i <- subset(limma_results_i, !is.na(logFC)) %>%
      mutate(condition=names(comparisons[comp]))
    limma_results_all[[names(comparisons[comp])]] <- limma_results_i
  }
  limma_results_all <- reduce(limma_results_all,bind_rows)

  # names(limma_results_all)[grep("P.Value", names(limma_results_i))] <- "pvalue.limma"
  # names(limma_results_all)[grep("adj.P.Val", names(limma_results_i))] <- "fdr.limma"
  limma_results_all <- subset(limma_results_all, !is.na(logFC))
  fdr_res <- NULL

  try(fdr_res <- fdrtool::fdrtool(limma_results_all$t, plot = FALSE, verbose = FALSE))
  if (!is.null(fdr_res))  {
    limma_results_all$pvalue.fdrtool <- fdr_res$pval
    limma_results_all$qval.fdrtool <- fdr_res$qval
    limma_results_all$lfdr.fdrtool <- fdr_res$lfdr
  }else  {
    limma_results_all$pvalue.fdrtool <- 1
    limma_results_all$qval.fdrtool <- 1
    limma_results_all$lfdr.fdrtool <- 1
  }
  limma_results_all$comparison <- limma_results_all$condition
  limma_results <- limma_results_all
  # limma_results <- bind_rows(limma.results, limma_results_all)
  rm(limma_results_all, fdr_res)


  # LimmaFDRcomparison
  limmaFDR_comp <- ggplot(data = limma_results, aes(abs(t))) +
    geom_line(aes(y = fdr.limma, colour = "limma - adj.P.Val")) +
    geom_line(aes(y = qval.fdrtool, colour = "fdrtool - qval")) +
    geom_point(aes(y = fdr.limma, colour = "limma - adj.P.Val")) +
    geom_point(aes(y = qval.fdrtool, colour = "fdrtool - qval")) +
    facet_wrap( ~ comparison, scale = "free_x") +
    ylab("fdr") +
    theme_bw()

  limmaFDR_comp2 <- ggplot(data = limma_results) +
    geom_histogram(aes(pvalue.limma, alpha = 0.5, fill = "limma"), bins = 40) +
    geom_histogram(aes(pvalue.fdrtool, alpha = 0.5, fill = "fdrtool"), bins = 40) +
    guides(alpha = FALSE) +
    xlab("p-value") +
    facet_wrap( ~  comparison, scale = "free_y") +
    coord_cartesian(xlim = c(0, 1)) +
    theme_bw()

  ## ----LimmaHitAnnotation-------------------------------------------------------------------------------------------------------------------------------------
  limma_results$hit_annotation_method <- NA

  limma_results$pvalue <- NA
  limma_results$fdr <- NA
  limma_results$comparison.old <- limma_results$comparison
  limma_results$comparison <- with(limma_results, comparison)
  for (comparison in unique(limma_results$comparison))
  {
    limma_hits <-
      nrow(limma_results[limma_results$comparison == comparison &
                           limma_results$fdr.limma <= fdr_hit_threshold, ])
    fdrtool_hits <-
      nrow(limma_results[limma_results$comparison == comparison &
                           limma_results$qval.fdrtool <= fdr_hit_threshold, ])
    if (limma_hits >= fdrtool_hits)
    {
      limma_results$hit_annotation_method[limma_results$comparison == comparison] <-
        "limma"
    }
    if (fdrtool_hits > limma_hits)
    {
      limma_results$hit_annotation_method[limma_results$comparison == comparison] <-
        "fdrtool"
    }
    rm(limma_hits, fdrtool_hits)
  }
  # table(limma_results$hit_annotation_method)
  limma_results$hit_annotation_method <- "limma"

  limma_results$pvalue[limma_results$hit_annotation_method == "limma"] <-
    limma_results$pvalue.limma[limma_results$hit_annotation_method == "limma"]
  limma_results$fdr[limma_results$hit_annotation_method == "limma"] <-
    limma_results$fdr.limma[limma_results$hit_annotation_method == "limma"]

  limma_results$pvalue[limma_results$hit_annotation_method == "fdrtool"] <-
    limma_results$pvalue.fdrtool[limma_results$hit_annotation_method == "fdrtool"]
  limma_results$fdr[limma_results$hit_annotation_method == "fdrtool"] <-
    limma_results$qval.fdrtool[limma_results$hit_annotation_method == "fdrtool"]

  limma_results$hit <-
    with(limma_results, ifelse(fdr <= fdr_hit_threshold & abs(logFC) >=
                                 fc_hit_threshold, TRUE, FALSE))
  limma_results$hit_annotation <-
    with(limma_results,
         ifelse(fdr <= fdr_hit_threshold & logFC >= fc_hit_threshold, "hit",
                ifelse(fdr <= fdr_candidate_threshold & logFC >= fc_candidate_threshold, "candidate",
                       ifelse(fdr <= fdr_hit_threshold & logFC <= -fc_hit_threshold, "neghit",
                              ifelse(fdr <= fdr_candidate_threshold & logFC <= -fc_candidate_threshold, "negcandidate",
                                     "nohit")))))
  limma_results$hit_annotation <-
    factor(limma_results$hit_annotation,
           ordered = TRUE,
           levels = c("hit", "candidate", "nohit","neghit","negcandidate"))

  # Keep only enrichedhits
  # ctrl.comparisons <- "Lysate-Incell"
  # limma_results$hit[
  #   with(limma_results, comparison %in% ctrl.comparisons &
  #          hit & logFC < 0)
  # ] <- FALSE
  # limma_results$hit_annotation[
  #   with(limma_results, comparison %in% ctrl.comparisons &
  #          hit_annotation %in% c("hit", "candidate") &
  #          logFC < 0)
  # ] <- "nohit"

  limma_results$hit_annotation <- as.character(limma_results$hit_annotation)
  limma_results$hit_annotation.orig <- limma_results$hit_annotation


  hitcounts <- with(limma_results, table(comparison, hit_annotation))
  # limma_results$comparison <- limma_results$comparison.old
  # limma_results$comparison.old <- NULL

  volcano_plots<- ggplot(data = limma_results, aes(logFC, -log10(fdr), colour = hit_annotation)) +
    geom_vline(aes(xintercept = 0)) +
    geom_point() +
    geom_text(aes(label = gsub("[|].+", "", id)),
              data = subset(limma_results, hit_annotation != "nohit"),
              vjust = 0, nudge_y = 0.1, size = 2, check_overlap = TRUE, show.legend=FALSE) +
    facet_wrap( ~ comparison + hit_annotation_method) +
    xlab("log2(fold change)") +
    theme_bw() +
    scale_color_manual(values=c(candidate="green4",hit="green",negcandidate="pink2",neghit="red",nohit="black"))

  volcano_plots.separate <- list()
  for(comp in unique(limma_results$comparison)) {
    limma_resuts.separate <- limma_results%>%filter(comparison==comp)
    volcano_plots.separate[[comp]]<- ggplot(data = limma_results%>%filter(comparison==comp), aes(logFC, -log10(fdr), colour = hit_annotation)) +
      geom_vline(aes(xintercept = 0)) +
      geom_point() +
      geom_text(aes(label = gsub("[|].+", "", id)),
                data = subset(limma_resuts.separate, hit_annotation != "nohit"),
                vjust = 0, nudge_y = 0.1, size = 2, check_overlap = TRUE, show.legend=FALSE) +
      xlab("log2(fold change)") +
      theme_bw() +
      scale_color_manual(values=c(candidate="green4",hit="green",negcandidate="pink2",neghit="red",nohit="black"))
  }


  ## ----LimmaMA_plot, fig.height = 7---------------------------------------------------------------------------------------------------------------------------
  MA_plot <- ggplot(data = limma_results, aes(AveExpr, logFC, colour = hit_annotation)) +
    geom_hline(aes(yintercept = 0)) +
    geom_point() +
    geom_text(aes(label = gsub("[|].+", "", id)),
              data = subset(limma_results, hit_annotation != "nohit"),
              vjust = 0, nudge_y = 0.1, size = 2, check_overlap = TRUE) +
    facet_wrap( ~ comparison + hit_annotation_method) +
    scale_color_manual(values=c(candidate="green4",hit="green",negcandidate="pink2",neghit="red",nohit="black")) +
    xlab("average log2(signal_sum)") +
    ylab("log2(fold change)") +
    theme_bw()

  output <- list(
    data=limma_results,
    hitcounts=hitcounts,
    limma_vs_FDR=list(limmaFDR_comp,limmaFDR_comp2),
    volcano=volcano_plots,
    volcano.separate=volcano_plots.separate,
    MA=MA_plot
  )
  return(output)
}
