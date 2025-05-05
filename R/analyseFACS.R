analyzeFACS <- function() {
  
  require(flowCore)
  require(flowGate)
  require(flowWorkspace)
  require(ncdfFlow)
  require(tidyverse)
  require(data.table)
  require(toOrdinal)
  
  FCSfolder <- choose.dir(default=".", caption="Select a folder with fcs files.")
  fcsfiles <- list.files(FCSfolder,pattern=".fcs$",ignore.case=TRUE,full.names=TRUE)
  l.fcs <- length(fcsfiles)
  if(l.fcs==0) {
    stop("There are no FCS files in the specified FCSfolder.")
  } else {
    cat(paste0("Found ", l.fcs, " FCS files.\n"))
  }
  
  gatingsample <- choose.files(default=FCSfolder, "Choose the FCS file that will be used for gating",
                               multi=FALSE,filters="fcs$|FCS$") %>%
    str_split_1("\\\\") %>% .[length(.)]
  
  
  saveRAM <- askYesNo("Do you want to save RAM when handling the data? Recommended for large data sets.")
  
 
  if(saveRAM) {
    fs <- read.ncdfFlowSet(
      files=fcsfiles
    )
  } else {
    fs <- read.flowSet(path = FCSfolder,
                       pattern = ".fcs$|.FCS$",
                       full.names = TRUE) 
  }
  gs <- GatingSet(fs)
  
  spill.mat <- spillover(fs[[1]]) %>% .[!sapply(., is.null)]
  if(length(spill.mat)==1) {
    spill.mat <- spill.mat[[1]]
  }
  
  comp <- compensation(spill.mat)
  gsc <- compensate(gs, comp)
  
  cat("The following channels are available:\n")
  cat("Channel names from the data: ", paste(colnames(gsc),collapse=", "),"\n")
  cat("Channel names explainer:\n")
  print(markernames(gsc))
  cat("Please choose the channels to use for gating in the order of gating.\n")
  cat("(follow the format Gate name, Parent gate: channel-on-x-axis, channel-on-y-axis etc)\n")
  
  channel_redoing <- TRUE
  while(channel_redoing) {
    channels <- c()
    channeling <- TRUE; i<-1
    while(channeling) {
      cat(paste0("Define the ",toOrdinal(i)," gate or insert empty line when done."))
      channels[i] <- readline()
      if(channels[i]=="") {
        channeling=FALSE
      } else {
        i=i+1
      }
    }
    
    ctable <- data.frame(ch=channels) %>%
      slice_head(n=-1) %>%
      separate_wider_delim(cols=ch,delim=":",names=c("Gates","Channels")) %>%
      separate_wider_delim(cols=Gates,delim=",",names=c("Gate","Parent"),too_few="align_start") %>%
      separate_wider_delim(cols=Channels,delim=",",names=c("ChannelX","ChannelY"),too_few="align_start") %>%
      mutate(across(everything(),function(x) str_replace(str_trim(x)," ","."))) %>%
      mutate(Parent=ifelse(is.na(Parent),"root",Parent)) %>%
      mutate(ChannelY=ifelse(is.na(Parent),NULL,ChannelY))
    
    
    if(any(duplicated(ctable$Gate))) {
      dup.gate <- ctable$Gate[duplicated(ctable$Gate)]
      message("The following gates have duplicate name. Please make sure that the gates have the same names.")
      cat(paste0("Duplicated gate names: ", paste(unique(dup.gate),collapse=", ")),"\n")
      cat("Please redefine the gates with unique gate names.\n")
      channels <- c(); i=1
    } else {
      channel.names <- c(ctable$ChannelX,ctable$ChannelY) %>% unique()
      if(any(!channel.names %in% colnames(gsc))) {
        unknown.channels <- channel.names[!channel.names %in% colnames(gsc)] %>% unique() %>%paste(collapse=", ")
        message("You chose unknown channels for gating.")
        cat(paste0("Unknown channels: ", unknown.channels,"\n"))
        cat("Channel names from the data: ", paste(colnames(gsc),collapse=", "),"\n")
        cat("Channel names explainer:\n")
        print(markernames(gsc))
        cat("Please redo the gates using channel names from the data.\n")
        channels <- c(); i=1
      } else {
        cat("Following channels will be used for gating, in this order.")
        print(ctable)
        isok <- menu(choices=c("Yes","No"),title="Is this correct?")
        if(isok==1) {
          channel_redoing=FALSE
        } else {
          message("No problem, let's redo the gates.")
          channels <- c(); i=1
        }
      }
    }
  }
  
  # gating
  gates <- list()
  gsample <- which(sampleNames(gsc)==gatingsample)
  for(i in 1:nrow(ctable)) {
    gates[[i]] <- gs_gate_interactive(gsc,
                                      filterId = ctable$Gate[i],
                                      sample=gsample,
                                      subset=ctable$Parent[i],
                                      dims = list(ctable$ChannelX[i],ctable$ChannelY[i]))
  }
  
  regate <- askYesNo("Are you happy with the gates? If you want to redo them, click No.")
  while(!regate) {
    gates <- list()
    for(i in 1:nrow(ctable)) {
      gates[[i]] <- gs_gate_interactive(gsc,
                                        filterId = ctable$Gate[i],
                                        sample=gsample,
                                        subset=ctable$Parent[i],
                                        dims = list(ctable$ChannelX[i],ctable$ChannelY[i]),
                                        regate=TRUE)
    }
    regate <- askYesNo("Are you happy with the gates? If you want to redo them, click No.")
  }
  
  # plotting
  
  gatepaths <- gs_get_pop_paths(gsc, path = "full") %>% .[.!="root"]
  ctable$Path = gatepaths
  

  cplots <- list()
  for(i in 1:nrow(ctable)) {
    
    cplots[[ctable$Path[i]]] <-
      ggcyto(gsc, aes(x =!!sym(ctable$ChannelX[i]), y = !!sym(ctable$ChannelY[i])), subset = ctable$Parent[i]) +
      geom_hex(bins = 200) +
      geom_gate(ctable$Path[i],colour="black",lwd=1) +
      scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)),
                    name=ctable$ChannelX[i]) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)),
                    name=ctable$ChannelX[i]) +
      theme_bw() +
      ggtitle(paste0("Gate ", ctable$Path[i])) +
      theme(panel.grid=element_blank())
  }
  
  stats <- gs_pop_get_count_fast(gsc,"count",format="wide") %>%
    t() %>% as.data.frame() %>%
    rownames_to_column("Sample")
  
  saveoutput <- askYesNo("Do you want to save the table and plots to your hard drive?")
  if(saveoutput) {
    save.dir <- choose.dir(default=".")
    for(plot in names(cplots)) {
      ggsave(paste0(save.dir,"\\",str_remove(str_replace_all(plot,"/","_"),"^_"),".png"), cplots[[plot]], dpi=600)
    }
    fwrite(stats,paste0(save.dir,"\\counts.csv"))
  }
  
  output <- list(
    counts=stats,
    plots=cplots
  )
  
  return(output)
}
