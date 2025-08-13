analyze.DIAPELSA <- function(file,
                            extract.after = "^",
                            extract.before = "_",
                            ctrl.name = "Ctrl",
                            FC.cutoff = c(1.2,1.2),
                            p.cutoff = c(0.01,0.05),
                            p.adjust = "BH"
) {

  diann_report <- arrow::read_parquet(file) %>%
    dplyr::filter(Lib.PG.Q.Value <= 0.01 & Lib.Q.Value <= 0.01) %>%
    dplyr::mutate(File.Name = Run) %>%
    dplyr::filter(str_detect(Protein.Ids, "cRAP|Biognosys", negate = TRUE))


  pep.mtx <- diann_report %>%
    dplyr::select(Run,Protein.Ids,Precursor.Id,Precursor.Normalised) %>%
    mutate(Run=str_extract(Run,paste0("(?<=",extract.after,").*"))) %>%
    unite("id",Protein.Ids,Precursor.Id,sep=";_;") %>%
    pivot_wider(id_cols=id,names_from=Run,values_from=Precursor.Normalised) %>%
    column_to_rownames("id")

  pep.raw <- pep.mtx %>% as.data.frame() %>% rownames_to_column("id") %>%
    separate_wider_delim(id,";_;",names=c("protein","peptide"))

  groups_for_design <- str_extract(colnames(pep.mtx), paste0(".*(?=", extract.before, ")"))
  design_matrix <- model.matrix(~ 0 + groups_for_design)
  colnames(design_matrix) <- levels(factor(groups_for_design))

  limma_model <- lmFit(log2(pep.mtx), design = design_matrix, method = "ls")

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
    ) %>%
    separate_wider_delim(Protein.Ids,";_;",names=c("protein.id","peptide"))

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

  return(list(
    raw_data = pep.raw,
    all_results = all_results,
    volcano_facet = volcano_facet))

}
