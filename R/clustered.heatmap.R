clustered.heatmap <- function(
    data=NULL,
    cl.group="all",
    plot.reverse=TRUE,
    plot.ddg.space=0.05
) {
  
  require(ggdendro)
  require(patchwork)
  
  get_legend<-function(myggplot) { # to extract legend from a single ggplot
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  data.prep <- prepare.data(data,cols=c("id","group","condition","value"),no.cols=TRUE)
  if(data.prep$message!="ok") {
    cat("Problem with pathways parameter:\n")
    cat(data.prep$message)
    if(data.prep$status=="error") {
      stop()
    }
  }
  data<- data.prep$data
  rm(data.prep)
  
  if(!cl.group %in% c(unique(data$group),"all")) {
    stop("cl.group must be all or one of the categories in the data frame.")
  }
  
  if(cl.group!="all") {
    widedata <- data %>%
      filter(group==cl.group)
  } else {
    widedata <- data
  }
  widedata <- widedata %>%
    pivot_wider(id_cols=id,names_from=c(group,condition),values_from=value) %>%
    column_to_rownames("id")
  
  hc <- hclust(dist(widedata))
  ord <- rownames(widedata)[hc$order]
  if(plot.reverse) {
    ord <- rev(ord)
  }
  ddg <- ggdendrogram(hc, rotate=FALSE, labels=FALSE, theme_dendro=TRUE, size=1) +
    theme_void() + coord_flip()
  ddg$layers[[2]]$aes_params$size <- 0.1
  
  plot <- data %>%
    mutate(id=factor(id,levels=ord)) %>%
    ggplot(aes(x=condition,y=id)) +
    geom_tile(aes(fill=value)) +
    scale_fill_gradient2(high="red3",mid="cornsilk",low="blue3") +
    facet_wrap(~group) +
    scale_x_discrete() + scale_y_discrete() +
    theme_bw(base_size = 12) +
    theme(axis.text.y=element_blank(),
          legend.position="bottom",
          panel.border=element_blank(),
          plot.margin = ggplot2::margin(0, 0, 0, 0),
          panel.spacing = unit(0, "lines"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()
    )
  legend <- get_legend(plot)
  plot <- plot +
    theme(legend.position="none")
  comp.plot <- plot + ddg + plot_layout(ncol=2, widths=c(1,plot.ddg.space))
  output <- list(
    plot=comp.plot,
    clustering=hc,
    dendrogram=ddg,
    heatmap=plot,
    legend=legend
  )
}
