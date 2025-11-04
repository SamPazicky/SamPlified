#' Differential Stability Analysis from DIA-NN Output
#'
#' This function performs differential stability analysis on protein abundance data from a DIA-NN `.parquet` file.
#' It uses the `limma` package to compute moderated statistics, generates volcano plots, and MD plots.
#'
#' The implementation is heavily inspired by code from the limma-for-proteomics
#' GitHub repository: \url{https://github.com/41ison/limma-for-proteomics}
#'
#' @param file Path to the DIA-NN `.parquet` report file for DIA data or to PD-exported txt file for DDA.
#' @param extract.after Regex pattern to extract group names from sample names (after a specific pattern).
#' @param extract.before Regex pattern to extract group names from sample names (before a specific pattern).
#' @param ctrl.name Name of the control group (e.g., "Ctrl").
#' @param pos.ctrl.name Regex pattern to extract the positive control sample.
#' @param pos.ctrl.id ID of the protein that is the hit in the positive sample.
#' @param exclude Vector of patterns, which will be matched to samples and those will be excluded from the analysis.
#' @param replicates Numeric. Number of replicates. Default is 3.
#' @param FC.cutoff Fold change cutoff for determining significance for hits (first value) and candidates (second value) (default is c(1.2,1.2)).
#' @param p.cutoff Adjusted p-value cutoff for determining significance for hits (first value) and candidates (second value) (default is c(0.01,0.05)).
#' @param p.adjust p-value adjustment method in limma::topTable
#' @param compare.to.all If TRUE, one treatment will be compared not only to the control sample but to all other treatment groups. Default is FALSE.
#' @param pulses In how many gas fractions were the samples measured?
#' @param pulse.quant Quantification method for GPF. "pept" is based on maximum peptide intensity, "prot.max" maximum protein intensity and "prot.mean" mean protein intensity.
#' @param report.quant.col Column in the report that contains quantification values. Either Genes.MaxLFQ.Unique or PG.MaxLFQ
#' @param export.xic Logical. Should the XIC of the hit proteins be plotted? This takes additional time (5-10s per hit).
#' @param xic.folder Folder in which the `xic.parquet` files are saved.
#' @param acq Acqusition mode. DIA or DDA?
#' @param DIA.filtering Numeric vector. Threshold for Lib.PG.Q.Value, Lib.Q.Value and PG.Q.Value, in this order.
#' @param TMT.labels String vector of TMT-labelled samples names in the same order as in the data file.
#' @param TMT String: either TMT10 or TMT16.
#'
#' @importFrom diann diann_matrix
#' @importFrom arrow read_parquet
#' @importFrom dplyr distinct filter select mutate %>% sym
#' @importFrom stringr str_split
#' @importFrom tidyr separate
#' @importFrom purrr map
#' @import limma
#' @import tidyverse
#' @importFrom ggrepel geom_text_repel
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{all_results}{A data frame of differential analysis results.}
#'   \item{volcano_facet}{A `ggplot` object showing faceted volcano plots.}
#'   \item{volcano_list}{A named list of individual volcano plots per contrast.}
#'   \item{MD_plot_faceted}{A `ggplot` object showing faceted MD plots.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- analyze_stability_from_parquet(
#'   file = "data/report.parquet",
#'   extract.after = "DIA_",
#'   extract.before = "_Rep",
#'   ctrl.name = "Ctrl",
#'   FC.cutoff = 1.5,
#'   p.cutoff = 0.01
#' )
#' result$all_results
#' result$volcano_facet
#' }

analyze.DIAPISA <- function(file,
                            extract.after = "^",
                            extract.before = "_",
                            ctrl.name = "Ctrl",
                            pos.ctrl.name=NA,
                            pos.ctrl.id="PF3D7_0417200.1-p1",
                            exclude=pos.ctrl.name,
                            replicates=3,
                            FC.cutoff = c(1.2,1.2),
                            p.cutoff = c(0.01,0.05),
                            p.adjust = "BH",
                            compare.to.all=FALSE,
                            pulses=1,
                            pulse.quant="prot.max",
                            report.quant.col="Genes.MaxLFQ.Unique",
                            export.xic=FALSE,
                            xic.folder="report_xic",
                            acq="DIA",
                            DIA.filtering=c(0.01,0.01,0.01),
                            TMT.labels=NA,
                            TMT="TMT16"
) {

  if(acq=="DDA") {
    pulses=1
    export.xic=FALSE
    if(is.na(TMT.labels[1])) {
      stop("TMT.labels must not be NA if acq is DDA.")
    }
  } else if(acq!="DIA") {
    stop("acq can only be either DDA or DIA!")
  }

  exclude <- c(pos.ctrl.name,exclude) %>% unique()

  if(acq=="DIA") {
    diann_report <- arrow::read_parquet(file)
    samples <- diann_report %>% pull(Run) %>% unique()
    diann_report <- diann_report %>%
      dplyr::filter(Lib.PG.Q.Value <= DIA.filtering[1] & Lib.Q.Value <= DIA.filtering[2] & PG.Q.Value <= DIA.filtering[3]) %>%
      dplyr::mutate(File.Name = Run) %>%
      dplyr::filter(str_detect(Protein.Ids, "cRAP|Biognosys", negate = TRUE))
    samples_out <- setdiff(samples,unique(diann_report$Run))
    if(length(samples_out)!=0) {
      stopmessage <- paste0("Some samples fell out of analysis after DIA data filtering! Samples: ",paste(samples_out,collapse=", "))
      stop(stopmessage)
    }

    cnt <- max(diann_report$Run.Index)+1
    if(pulses==1) {
      cat("Analyzing", cnt/pulses, "samples run in DIA mode.\n")
    } else {
      cat("Analyzing", cnt/pulses, "samples run in ",pulses," pulses...\n")

    }

    prot_mtx <- diann::diann_matrix(diann_report,
                                    id.header = "Protein.Ids",
                                    quantity.header = report.quant.col,
                                    proteotypic.only = TRUE,
                                    pg.q = DIA.filtering[3])
  } else if(acq=="DDA") {
    if(TMT=="TMT10") {
      prot_mtx <- fread(file) %>%
        dplyr::select(Accession, contains("Abundances")) %>%
        dplyr::select(!contains("CV",ignore.case=FALSE)) %>%
        column_to_rownames("Accession")
    } else if(TMT=="TMT16") {
      prot_mtx <- fread(file) %>%
        dplyr::select(Accession, contains("Abundance:")) %>%
        column_to_rownames("Accession")
    } else {
      stop("TMT must be TMT10 or TMT16!")
    }

    if(ncol(prot_mtx)==length(TMT.labels)) {
      prot_mtx <- prot_mtx %>% setNames(TMT.labels)
    } else {
      stop("There are more samples than names given in TMT.labels argument!")
    }
  }


  no.pro <- nrow(prot_mtx)
  cat(paste0("Detected ",no.pro," proteins.\n"))
  if(!is.null(exclude)) {
    prot_mtx_filtered <- prot_mtx[,!str_detect(colnames(prot_mtx),paste(exclude,collapse="|"))]
  } else {
    prot_mtx_filtered <- prot_mtx
  }

  excluded <- ncol(prot_mtx)-ncol(prot_mtx_filtered)
  if(excluded==1) {
    cat(paste0("Excluded ",excluded," sample from the analysis."))
  } else if(excluded>1) {
    cat(paste0("Excluded ",excluded," samples from the analysis."))
  }

  if(!is.na(pos.ctrl.name)) {
    sel.score <- prot_mtx[,str_detect(colnames(prot_mtx),paste(c(pos.ctrl.name,ctrl.name),collapse="|"))] %>%
      as.data.frame() %>%
      dplyr::select(contains(pos.ctrl.name),contains(ctrl.name)) %>%
      setNames(c(pos.ctrl.name,paste0("Ctrl",1:replicates))) %>%
      rownames_to_column("id") %>%
      pivot_longer(cols=!id,names_to="condition",values_to="abundance") %>%
      mutate(group=ifelse(str_detect(condition,"Ctrl"),"Ctrl","Pos.Ctrl")) %>%
      group_by(id,group) %>%
      dplyr::summarise(mu=mean(abundance),sd=sd(abundance), .groups="drop") %>%
      ungroup() %>%
      pivot_wider(id_cols=id, names_from=group, names_sep="_",values_from=c(mu,sd)) %>%
      mutate(z=(mu_Pos.Ctrl-mu_Ctrl)/sd_Ctrl) %>%
      filter(id==pos.ctrl.id) %>%
      pull(z) %>%
      round(1)
  } else {
    sel.score <- NA
  }


  cat(paste0("Selection score is ", sel.score,"!\n"))
  cat("Selection score must be above 2. The higher, the better!\n")

  #quantifying pulseDIA
  if(pulses>1) {
    if(pulse.quant=="prot.max") {
      prot_mtx_filtered_p <- prot_mtx_filtered %>%
        as.data.frame() %>%
        rownames_to_column("id") %>%
        # setNames(str_extract(names(.), paste0("(?<=",extract.after,").*"))) %>%
        pivot_longer(
          cols = !id,
          names_to = c("pulse", "sample"),
          names_pattern = ".*_pulseDIA\\d+x(\\d+)_(.+)",
          values_to = "quant"
        ) %>%
        na.omit() %>%
        group_by(id,sample) %>%
        dplyr::summarise(quant=max(quant,na.rm=TRUE),.groups="keep") %>%
        ungroup() %>%
        pivot_wider(id_cols=id,names_from=sample,values_from=quant) %>%
        column_to_rownames("id")
    } else if(pulse.quant=="prot.mean") {
      prot_mtx_filtered_p <- prot_mtx_filtered %>%
        as.data.frame() %>%
        rownames_to_column("id") %>%
        # setNames(str_extract(names(.), paste0("(?<=",extract.after,").*"))) %>%
        pivot_longer(
          cols = !id,
          names_to = c("pulse", "sample"),
          names_pattern = ".*_pulseDIA\\d+x(\\d+)_(.+)",
          values_to = "quant"
        ) %>%
        na.omit() %>%
        group_by(id,sample) %>%
        dplyr::summarise(quant=mean(quant,na.rm=TRUE),.groups="keep") %>%
        ungroup() %>%
        pivot_wider(id_cols=id,names_from=sample,values_from=quant) %>%
        column_to_rownames("id")
    } else if(pulse.quant=="pept") {
      prot_mtx_filtered_p <- diann_report %>%
        filter(Protein.Group != "") %>%
        select(Run, Protein.Group, Precursor.Quantity) %>%
        mutate(Sample = str_remove(Run, "[[:digit:]]x[[:digit:]]")) %>%
        mutate(Sample = str_extract(Sample, paste0("(?<=", extract.after, ").*"))) %>%
        group_by(Sample, Protein.Group) %>%
        summarise(Summed.Quantity = sum(Precursor.Quantity, na.rm = TRUE), .groups = "drop") %>%
        ungroup() %>%
        { if (!is.null(exclude)) filter(., !str_detect(Sample, exclude)) else . } %>%
        pivot_wider(id_cols=Protein.Group,names_from=Sample,values_from=Summed.Quantity) %>%
        column_to_rownames("Protein.Group")
    }

  } else {
    prot_mtx_filtered_p <- prot_mtx_filtered %>%
      as.data.frame() %>%
      setNames(str_extract(names(.), paste0("(?<=",extract.after,").*")))
  }

  # raw protein table
  prot_raw <- prot_mtx %>% as.data.frame() %>%
    rownames_to_column("id")


  groups_for_design <- str_extract(colnames(prot_mtx_filtered_p), paste0(".*(?=", extract.before, ")"))
  design_matrix <- model.matrix(~ 0 + groups_for_design)
  colnames(design_matrix) <- levels(factor(groups_for_design))

  prot_mtx_filtered_p_log <- prot_mtx_filtered_p %>% as.matrix() %>% vsn::justvsn() %>% as.data.frame()

  limma_model <- lmFit(prot_mtx_filtered_p_log, design = design_matrix, method = "ls")

  all_groups <- colnames(design_matrix)
  comparisons <- setdiff(all_groups, ctrl.name)

  contrast_formulae <- setNames(
    lapply(comparisons, function(grp) paste0(grp, "-", ctrl.name)),
    comparisons
  )

  contrast_matrix <- makeContrasts(
    contrasts = unlist(contrast_formulae),
    levels = design_matrix
  )

  if(compare.to.all) {
    contrast_matrix[contrast_matrix == 0|contrast_matrix==-1] <- -1/ncol(contrast_matrix)
  }

  estimated_coef <- contrasts.fit(limma_model, contrast_matrix)
  empirical_Bayes_fit <- eBayes(estimated_coef)

  contrast_names <- colnames(contrast_matrix)

  all_results <- purrr::map_dfr(contrast_names, function(contrast) {
    topTable(empirical_Bayes_fit,
             coef = contrast,
             number = Inf,
             adjust.method = "BH") %>%
      tibble::rownames_to_column("Protein.Ids") %>%
      dplyr::mutate(Contrast = contrast)
  }) %>%
    dplyr::mutate(status = case_when(
      logFC > log2(max(FC.cutoff)) & adj.P.Val <= max(p.cutoff) ~ "Stabilized",
      logFC < -log2(max(FC.cutoff)) & adj.P.Val <= max(p.cutoff) ~ "Destabilized",
      TRUE ~ "Not significant"),
      status = factor(status, levels = c("Destabilized", "Not significant", "Stabilized"))
    ) %>%
    dplyr::mutate(hit = case_when(
      abs(logFC) > log2(FC.cutoff[1]) & adj.P.Val <= p.cutoff[1] ~ "Hit",
      abs(logFC) > log2(FC.cutoff[2]) & adj.P.Val <= p.cutoff[2] ~ "Candidate",
            TRUE ~ ""),
      hit = factor(hit, levels = c("Hit", "Candidate", ""))
    ) %>%
    dplyr::mutate(label = case_when(
      hit != "" ~ paste(hit, status, sep = " - "),
      TRUE ~ status),
      label=as.character(label)
      )

  max_logfc <- ceiling(max(abs(all_results$logFC), na.rm = TRUE))

  volcano_facet <- ggplot(all_results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(alpha = label, color = label), show.legend = c(alpha = FALSE)) +
    scale_alpha_manual(values = c("Destabilized" = 1, "Stabilized" = 1, "Not significant" = 0.3), drop=TRUE) +
    geom_hline(yintercept = -log10(p.cutoff), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-log2(FC.cutoff), log2(FC.cutoff)), linetype = "dashed", color = "black") +
    facet_wrap(~ Contrast) +
    scale_color_manual(values = c("Candidate - Destabilized" = "steelblue2",
                                  "Not significant" = "grey",
                                  "Hit - Destabilized" = "steelblue4",
                                  "Hit - Stabilized" = "firebrick",
                                  "Candidate - Stabilized" = "firebrick2"),
                       drop=TRUE) +
    labs(title = "Volcano plots of limma contrasts",
         x = "log2(Stability Fold Change)",
         y = "-log10(p-value)") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    ) +
    xlim(-max_logfc, max_logfc)

  volcano_list <- purrr::map(unique(all_results$Contrast), function(contrast_name) {

    df <- all_results %>% filter(Contrast == contrast_name)

    count_stabilized <- df %>% filter(status == "Stabilized") %>% nrow()
    count_destabilized <- df %>% filter(status == "Destabilized") %>% nrow()

    max_logfc <- ceiling(max(abs(df$logFC), na.rm = TRUE))
    ymax <- max(-log10(df$adj.P.Val), na.rm = TRUE) + 0.5
    y_offset <- ymax * 0.05
    text_y_pos <- ymax + y_offset

    ggplot(df, aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point(aes(alpha = status, color = status)) +
      scale_alpha_manual(values = c("Destabilized" = 1, "Stabilized" = 1, "Not significant" = 0.3)) +
      geom_hline(yintercept = -log10(p.cutoff), linetype = "dashed", color = "black") +
      geom_vline(xintercept = c(-log2(FC.cutoff), log2(FC.cutoff)), linetype = "dashed", color = "black") +
      annotate("text", x = -max_logfc, y = text_y_pos,
               label = paste0("Destabilized: ", count_destabilized),
               hjust = 0, fontface = "bold", size = 4) +
      annotate("text", x = max_logfc, y = text_y_pos,
               label = paste0("Stabilized: ", count_stabilized),
               hjust = 1, fontface = "bold", size = 4) +
      ggrepel::geom_text_repel(
        data = df %>% dplyr::filter(status != "Not significant"),
        aes(label = Protein.Ids),
        size = 3,
        max.overlaps = 15,
        box.padding = 0.5,
        segment.color = "black"
      ) +
      scale_color_manual(values = c("Destabilized" = "steelblue",
                                    "Not significant" = "grey",
                                    "Stabilized" = "firebrick")) +
      labs(title = contrast_name,
           x = "log2(Stability fold Change)",
           y = "-log10(p-value)") +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold")
      ) +
      xlim(-max_logfc, max_logfc)
  }) %>%
    rlang::set_names(unique(all_results$Contrast))

  label_data <- all_results %>%
    group_by(Contrast) %>%
    summarize(
      max_expr = max(AveExpr, na.rm = TRUE),
      stabilized = sum(status == "Stabilized", na.rm = TRUE),
      destabilized = sum(status == "Destabilized", na.rm = TRUE),
      max_fc = max(abs(logFC), na.rm = TRUE),
      .groups = "drop"
    )

  MD_plot_faceted <- ggplot(all_results, aes(x = AveExpr, y = logFC, color = status)) +
    geom_point(aes(alpha = status)) +
    scale_alpha_manual(values = c("Destabilized" = 1, "Stabilized" = 1, "Not significant" = 0.3)) +
    facet_wrap(~ Contrast) +
    scale_color_manual(values = c("Destabilized" = "steelblue",
                                  "Not significant" = "grey",
                                  "Stabilized" = "firebrick")) +
    labs(title = "Mean-Difference (MD) Plots",
         x = "Average Expression",
         y = "log2(Stability Fold Change)") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    )

  allhits <- all_results %>% filter(hit %in% c("Hit","Candidate")) %>% pull(Protein.Ids)
  rawhits <- prot_mtx_filtered_p_log %>%
    rownames_to_column("id") %>%
    filter(id %in% allhits) %>%
    pivot_longer(cols=!id,names_to="samplerep",values_to="Abundance") %>%
    filter(!samplerep %in% exclude) %>%
    separate_wider_delim(samplerep,extract.before,names=c("Sample","Replicate")) %>%
    group_by(id,Sample) %>%
    dplyr::mutate(meanAbundance=mean(Abundance,na.omit=TRUE)) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(
      relAbundance = Abundance / mean(meanAbundance[Sample == ctrl.name],na.rm=TRUE)
    ) %>%
    ungroup()
  barplots <- list()

  barplots <- rawhits %>%
    # filter(Sample!=ctrl.name) %>%
    group_by(id,Sample) %>%
    mutate(meanrelAbundance=mean(relAbundance,na.rm=TRUE)) %>%
    ungroup() %>%
    split(.$id) %>%
    lapply(function(df) {
      df$Sample <- fct_relevel(df$Sample, ctrl.name)
      ggplot(df) +
        geom_bar(data=.%>%distinct(Sample,meanrelAbundance,.keep_all=TRUE),
                 aes(x=Sample,y=meanrelAbundance,na.rm=TRUE),
                 width=0.8, stat="identity",color="black",fill="gray90") +
        geom_jitter(aes(x=Sample,y=relAbundance,color=Replicate),width=0.15,size=2) +
        scale_y_continuous(expand=c(0,0),
                           limits=c(0,max(rawhits$relAbundance, na.rm=TRUE)+0.1)
        ) +
        ggtitle(unique(df$id)) +
        theme_bw() +
        theme(legend.position="none",
              plot.title=element_text(hjust=0.5))
    })

  #XIC charts
  #xic.folder = "report_xic"
  if(export.xic) {
    hitdata <- list()
    for(x in list.files(xic.folder)) {
      xicreport <- read_parquet(paste0(xic.folder,"/",x))
      samp <- str_remove(x,"\\.xic\\.parquet")

      hitdata[[samp]] <- diann_report %>%
        as.data.frame() %>%
        dplyr::select(Protein.Ids,Precursor.Id) %>%
        filter(Protein.Ids %in% allhits) %>%
        setNames(c("id","pr")) %>%
        distinct(.keep_all=TRUE) %>%
        right_join(xicreport,by="pr") %>%
        na.omit()
    }
    hitdata <- hitdata %>%
      rbindlist(idcol="sample") %>%
      filter(feature!="index") %>%
      group_by(id, sample, pr, feature) %>%
      arrange(rt, .by_group = TRUE) %>%
      mutate(
        # index for row-based window
        idx = row_number(),
        # mark rows with value > 1e3
        hit = value > 1e3,
        # row index of "max value" within the group
        idx_max = which.max(value),
        rt_max  = rt[idx_max],
        # define row-based window (±3 rows)
        keep_rows = idx %in% unlist(lapply(which(hit), function(i) (i-3):(i+3))),
        # define rt-based window (±2.5 around rt_max)
        keep_rt   = abs(rt - rt_max) <= 2.5,
        # combine conditions
        keep = keep_rows | keep_rt
      ) %>%
      filter(keep) %>%
      select(-hit, -idx, -idx_max, -rt_max, -keep_rows, -keep_rt, -keep) %>%
      ungroup() %>%
      mutate(feature=str_remove(feature,"-unknown")) %>%
      mutate(sample=str_extract(sample,paste0("(?<=",extract.after,").*")))


    samples <- unique(hitdata$sample)
    top_strip <- hitdata %>%
      distinct(sample) %>%                 # get unique samples
      ggplot(aes(x = sample, y = 1, label = sample)) +
      geom_text(size = 4) +                # facet strip labels
      theme_void() +                       # remove axes and background
      theme(
        plot.margin = margin(0, 0, 0, 0),
        panel.spacing = unit(0, "cm")
      )

    features <- unique(hitdata$feature)
    fcolors <- rainbow(length(features)) %>% setNames(features)

    dir.create("XICs",showWarnings = FALSE)
    for(hid in unique(hitdata$id)) {
      c.hitdata <- hitdata %>%
        filter(id==hid) %>%
        complete(pr, sample, feature, fill = list(value = NA, rt = NA))
      prs <- unique(c.hitdata$pr) %>% length()
      hitplots <- c.hitdata %>%
        split(.$pr) %>%
        lapply(function(x) {
          pr_label <- unique(x$pr)
          ggplot(x, aes(x = rt, y = value, color = feature)) +
            geom_point(size=1) +
            geom_line(aes(group = feature), linewidth=0.1) +
            facet_wrap(~ sample, nrow = 1, drop=FALSE) +
            scale_color_manual(values = fcolors) +
            labs(title = pr_label) +
            theme_bw() +
            theme(
              legend.position = "none",
              strip.text = element_blank(),        # remove the text
              strip.background = element_blank(),   # remove the background
              plot.title = element_text(size = 8),
              panel.grid.major=element_line(linewidth=0.2),
              panel.grid.minor=element_line(linewidth=0.1),
              axis.text.y=element_text(size=6),
              axis.text.x=element_text(size=6,angle=90,hjust=0,vjust=0.5),
              axis.title=element_blank()
            )
        })
      hitplot <- patchwork::wrap_plots(hitplots,ncol=1)
      hitplot <- top_strip/hitplot + plot_layout(heights=c(0.3,prs))

      ggsave(paste0("XICs/",hid,".pdf"),hitplot,width=2*length(samples), height=4*prs,units="cm")
    }
  }


  return(list(
    raw_data = prot_raw,
    all_results = all_results,
    volcano_facet = volcano_facet,
    volcano_list = volcano_list,
    MD_plot_faceted = MD_plot_faceted,
    sel.score = sel.score,
    raw.bar.plots = barplots
  ))
}
