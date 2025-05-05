setwd("../210930_SLO_parasite_enrichment_trial3")
library(flowViz)
library(flowCore)
library(flowStats)
library(flowUtils)
library(geneplotter)
library(ggcyto)
library(gridExtra)
library(scales)
library(ggpointdensity)
library(viridis)
library(drc)
library(dplyr)
library(tidyr)


######
# USER IMPUT
######



#define general plot design and variables
customPlot <- list(
  theme_bw(base_size = 12),
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
)
scientific_10 <- function(x) {   parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))) }


#sample information
samples<-c("Control",0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0)
repeats<-c("noAbMCM","MCM")
HU<-12.33
save_flowcharts <- "REDUCED" 
#if YES, the flow cytometry charts will be saved but this is computationally-intensive process and it takes 5 min per chart
#if REDUCED, the flow cytometetry chart will be plotted only from first N points, hoping that they are representative enough.
reduced_n <- 5000 #this sets the first n charts number
#if NO (or anything else), the charts will not be plotted

#########
# GATES #
#########

# uRBC gate
gate_RBCs <- as.matrix(log10(read.csv("gate_RBCs.csv", header=FALSE)))
colnames(gate_RBCs) <- c("FSC.H","SSC.H")
RBCGate <- polygonGate(filterId="RBCs", .gate= gate_RBCs)

# single cells gate
gate_SCs <- as.matrix(log10(read.csv("gate_SCs.csv", header=FALSE)))
colnames(gate_SCs) <- c("FSC.A","FSC.H")
SCGate <- polygonGate(filterId="Single_cells", .gate= gate_SCs)

#infected RBCs
gate_iRBCs <- as.matrix(log10(read.csv("gate_iRBCs.csv", header=FALSE)))
colnames(gate_iRBCs) <- c("FSC.A","FL1.A")
iRBCGate <- polygonGate(filterId="iRBCs", .gate= gate_iRBCs) 

#uninfected RBCs
gate_uRBCs <- as.matrix(log10(read.csv("gate_uRBCs.csv", header=FALSE)))
colnames(gate_uRBCs) <- c("FSC.A","FL1.A")
uRBCGate <- polygonGate(filterId="uRBCs", .gate= gate_iRBCs) 

used_gates<-c(RBCGate,SCGate,iRBCGate)
subset_gates<-c(RBCGate,SCGate)
calc_gate<-setdiff(used_gates,subset_gates)

sbgates<-NULL
for (i in subset_gates) {
  if (is.null(sbgates)==TRUE) {
    sbgates<-i
  } else{
    sbgates<-sbgates&i
  }
}
rm(i)

clcgate<-NULL
for (i in calc_gate) {
  if (is.null(clcgate)==TRUE) {
    clcgate<-i
  } else{
    clcgate<-clcgate&i
  }
}
rm(i)

#initiating variables
data_counts <- data.frame(samples)

#read the files and construct flowSet
pathstart<-"20210930"
pathcore<-""
dirs<-list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
dirs <- dirs[grepl(paste0(pathstart,".*",pathcore,".*"), dirs) ]

for (j in seq_along(dirs)) {
  filepath=dirs[[j]]
  fclist <- list.files(path=filepath, pattern="*.fcs", ignore.case=TRUE)
  fs <-read.flowSet(fclist, path=filepath, transformation = FALSE, alter.names=TRUE) #, column.pattern=".H"
  
  pData(fs)$SLO <-as.character(samples)
  
  #compensate with spillover matrices
  
  comp <- fsApply(fs, function(x) spillover(x)[[3]], simplify=FALSE)
  fs_comp <- compensate(fs, comp)
  
  
  # transformation
  
  aTrans <- transformList(as.vector(colnames(summary(fs_comp)[[1]])),truncateTransform(a=1))
  lTrans <- transformList(as.vector(colnames(summary(fs_comp)[[1]])),log10)
  fs_trans <- transform(fs_comp, aTrans)
  fs_trans2 <- transform(fs_trans, lTrans)
  
  fs_chart <- fs_trans2
  #plot transformed data with gating
  gates<-NULL
  for (i in used_gates) {
    x=parameters(i)[[1]]
    y=parameters(i)[[2]]
    gate<-i@filterId
    
    if (save_flowcharts=="YES") {
      
      FSCHdata_gcyto <- ggcyto(fs_chart, aes_string(x=x, y=y))
      FSCHchart <- FSCHdata_gcyto +
        geom_pointdensity(fs_chart, mapping=aes_string(x=x, y=y), adjust=.1, size=.02, alpha=0.8) +
        facet_wrap(~SLO) +
        geom_gate(data=i, colour="springgreen3",size=0.75) +
        customPlot +
        scale_color_viridis(option="magma") +
        labs(colour = "Counts") +
        coord_cartesian(xlim=c(3,7),ylim=c(3,7)) +
        xlab(paste0("log ",x)) +
        ylab(paste0("log ",y))
      ggsave(paste0("FSCHcharts_","gate_",gate,"_",j,".png"),plot=FSCHchart,height=20,width=30,units="cm")
      
    } else if (save_flowcharts=="REDUCED") {
      
      sample_table <- data.frame(frame=sort(names(fs_chart@frames)), SLO=samples)
      chartdata <- data.frame()
      for (frame in sample_table$frame) {  
        framedata <- as.data.frame(fs_chart@frames[[frame]]@exprs)
        framedata$frame <- frame
        if(nrow(framedata)>reduced_n) {
          framedata <- slice_head(framedata,n=reduced_n)
        }
        framedata <- left_join(framedata,sample_table)
        chartdata <- bind_rows(chartdata, framedata)
      }
      
      FSCHchart <- ggplot(data=chartdata) +
        geom_pointdensity(mapping=aes_string(x=x, y=y), adjust=.1, size=.02, alpha=0.8) +
        facet_wrap(~SLO) +
        geom_polygon(data=as.data.frame(i@boundaries), mapping=aes_string(x=x,y=y), colour="springgreen3",size=0.75, fill=NA) +
        customPlot +
        scale_color_viridis(option="magma") +
        labs(colour = "Counts") +
        coord_cartesian(xlim=c(3,7),ylim=c(3,7)) +
        xlab(paste0("log ",x)) +
        ylab(paste0("log ",y))
      ggsave(paste0("FSCHcharts_","gate_",gate,"_",repeats[j],".png"),plot=FSCHchart,height=20,width=30,units="cm")
    } 
    
    fs_chart<-flowCore::Subset(fs_chart,i)
    if (is.null(gates)==TRUE) {
      gates<-i
    } else{
      gates<-gates&i
    }
  }
  rm(i)
  #applying gates
  
  RBCs_SCs <- flowCore::Subset(fs_trans2,sbgates)
  iRBC_counts = flowCore::filter(RBCs_SCs,clcgate)
  
  countall=vector()
  totalall=vector()
  for (i in seq_along(iRBC_counts)) {
    count <- summary(iRBC_counts[[i]])$true
    total <- summary(iRBC_counts[[i]])$n
    countall<-c(countall,count)
    totalall<-c(totalall,total)
  }
  rm(i)
  # Constituting a data frame
  newdata <- data.frame(countall,totalall)
  colnames(newdata) <- c(paste0("Count",j),paste0("Total",j))
  data_counts <- cbind(data_counts,newdata)
}
data_counts

# Calculating lysed percentages

data_counts <- data_counts %>% slice(-1)

for (i in seq_along(dirs)) {
  colname=paste0("Count",i,"_uRBCs")
  assign(colname,eval(parse(text=paste0("data_counts$Total",i,"-data_counts$Count",i))))
  data_counts <-cbind(data_counts, eval(parse(text=colname)))
  names(data_counts)[ncol(data_counts)] <- colname
  
  
  eval(parse(text=paste0(
    "data_counts <- data_counts %>% dplyr::mutate(Lysed",i,"iRBCs=(first(Count",i,")-(Count",i,"))/first(Count",i,"),
                                               Lysed",i,"Total=(first(Total",i,")-(Total",i,"))/first(Total",i,"),
                                               Lysed",i,"uRBCs=(first(Count",i,"_uRBCs)-(Count",i,"_uRBCs))/first(Count",i,"_uRBCs))"
  )))
  
}
data_counts <- data_counts %>% dplyr::mutate(Units=as.numeric(samples)*HU)


data_i <- data_counts %>% dplyr::select(starts_with("Lysed"))
mydata <- cbind(data_counts$Units, data_i)
colnames(mydata)[[1]] <- "Units"

# mydata<-mydata+min(mydata)
for (j in 1:length(dirs)) {
  ilist <- names(mydata%>%dplyr::select(setdiff(contains(as.character(j)),contains("Total"))))
  
  for (i in seq_along(ilist)) {
    fit <- drm(mydata[[ilist[i]]]~Units,data=mydata, fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50"),
                                                                fixed = c(NA,NA,NA,NA)))
    fitplot <- plot(fit, type="all")
    
    coefs <- setNames(
      fit$coefficients,
      c("hill", "min_value", "max_value","ec_50")
    )
    coefs 
    
    ic_50 <- with(
      as.list(coefs),
      exp(
        log(ec_50) + (1 / hill) * log(max_value / (max_value - 2 * min_value))
      )
    )
    ic_50 
    
    data.LL.4 <- drm(data = mydata,mydata[[ilist[i]]]~Units,fct=LL.4(fixed = c(NA,NA,NA,NA)),na.action = na.omit) # run model.
    # predictions and confidence intervals.
    data.fits <- expand.grid(conc=exp(seq(log(56.00e0), log(1.00e-01), length=100)))
    # new data with predictions
    pm <- predict(data.LL.4, newdata=data.fits, interval="confidence") 
    data.fits$p <- pm[,1]
    data.fits$pmin <- pm[,2]
    data.fits$pmax <- pm[,3]
    assign(paste0("data.fits",ilist[i]),data.fits)
  } 
  
  eval(parse(text=paste0(
    "iRBC_fits<-data.fits",ilist[i-1]
  )))
  eval(parse(text=paste0(
    "uRBC_fits <- data.fits",ilist[i] 
  )))
  EC50plot <- ggplot(mydata) +
    geom_point(mapping=aes(x=Units,y=get(ilist[1])*100,color=as.factor("iRBCs"))) +
    geom_point(mapping=aes(x=Units,y=get(ilist[2])*100,color=as.factor("uRBCs"))) +
    # geom_ribbon(data=iRBC_fits, aes(x=conc, y=p*100, ymin=pmin*100, ymax=pmax*100,fill=as.factor("iRBCs")), alpha=0.1) +
    geom_line(data=iRBC_fits, aes(x=conc, y=p*100, color=as.factor("iRBCs"))) +
    # geom_ribbon(data=uRBC_fits, aes(x=conc, y=p*100, ymin=pmin*100, ymax=pmax*100,fill=as.factor("uRBCs")), alpha=0.1) +
    geom_line(data=uRBC_fits, aes(x=conc, y=p*100, color=as.factor("uRBCs"))) +
    scale_color_manual(values=c("brown","red")) +
    scale_fill_manual(values=c("brown","red")) +
    # scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
    customPlot +
    xlab("Streptolysin O units (HU)") +
    ylab("Lysed RBCs (%)") +
    labs(color="Population", fill="Population") +
    xlim(c(0,max(mydata$Units*1.1))) +
    theme(legend.position=c(0.8,0.4))
  
  eval(parse(text=paste0(
    "mydata <- mydata %>% mutate(Parasitemia=data_counts$Count",j,"/data_counts$Total",j,"*100)"
  )))
  paraPlot <- ggplot(mydata) +
    geom_point(mapping=aes(x=Units,y=Parasitemia)) +
    geom_line(mapping=aes(x=Units,y=Parasitemia)) +
    customPlot +
    xlab(label="Streptolysin O units (HU)")+ylab("Parasitemia (%)") +
    ylim(c(0,100)) + xlim(c(0,max(mydata$Units*1.1)))
  # ggsave("parasitemia.png",plot=paraPlot, height=12, width=16, units="cm")
  
  combinedPlot<- grid.arrange(EC50plot,paraPlot,ncol=1)
  ggsave(paste0("combinedPlot",repeats[j],".png"),plot=combinedPlot, height=12, width=16, units="cm")
}
