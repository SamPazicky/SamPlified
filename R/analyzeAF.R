analyzeAF <- function(
  folder=NA
) {

  require(jsonlite)
  require(tidyverse)

  if(is.na(folder)) {
    stop("Please include a folder")
  } else if(!file.exists(folder)) {
    stop(paste0("The folder ", folder, " does not exist."))
  } else {
    if(is.na(str_extract(folder,"/$"))) {
      folder <- paste0(folder,"/")
    }
  }

  fd.files <- list.files(folder,pattern="full_data.*json$", full.names=TRUE) %>%
    gtools::mixedsort()
  jr.file <- list.files(folder,pattern="job_request.json$", full.names=TRUE)

  # get sequences etc.
  jrfile <- fromJSON(jr.file)
  sequences <- rep(
    jrfile$sequences[[1]][["proteinChain"]][["sequence"]],
    jrfile$sequences[[1]][["proteinChain"]][["count"]]
  )
  sequences.pasted <- paste(sequences,collapse="")
  sequences.lengths <- sapply(sequences, str_length) %>% unname()
  sequences.res <- data.frame(AA.seq=str_split_1(sequences.pasted,"")) %>%
    rownames_to_column("AA.no")
  sequences.ends <- cumsum(sequences.lengths)
  sequences.starts <- sequences.ends-sequences.lengths + 1
  sequences.rects <- lapply(seq_along(sequences.starts),
                                      function(x) list(
                                        xmin=sequences.starts[x],
                                        xmax=sequences.ends[x],
                                        ymin=sequences.starts[x],
                                        ymax=sequences.ends[x]))

  for(fdf in fd.files) {
    jfile <- fromJSON(fdf)

    cor.data <- jfile$pae %>% as.data.frame() %>%
      pivot_longer(cols=everything(),names_to="AA.no",values_to="PAE") %>%
      mutate(AA.no=str_remove(AA.no,"V")) %>%
      left_join(sequences.res,by="AA.no") %>%
      rename_with(~paste0("aligned",.x), starts_with("AA")) %>%
      rownames_to_column("AA.no") %>%
      left_join(sequences.res,by="AA.no") %>%
      rename_with(~paste0("scored",.x), starts_with("AA")) %>%
      dplyr::select(starts_with("scored"),starts_with("aligned"),PAE)

    cor.plot <- ggplot(cor.data) +
      geom_tile(aes(x=scoredAA.no,y=alignedAA.no, fill=PAE)) +
      scale_fill_gradient(low="darkgreen",high="white")
    for(i in seq_along(sequences.rects)) {
      cor.plot <- cor.plot +
        geom_rect(aes(!!!sequences.rects[[i]]))
    }

    ggsave("corplot.png",cor.plot, width=10,height=10,units="cm")


  }
}
