library(reshape2)
library(ez)




###################
# Set up
###################

# parDir <- "../Analyses/"
parDir <- "/run/user/1000/gvfs/afp-volume:host=10.5.72.28,user=exp,volume=Yorick/STT_reml/Analyses/"
writeGraph <- 1


### ROI variables
roiDir <- paste0(parDir,"roiAnalysis/")
nsDir <- paste0(roiDir,"ns_betas_RR1/")
ns_outDir <- paste0(roiDir,"ns_stats_RR1/")




###################
# Functions
###################

BehNames.Function <- function(dataString){
  if(dataString=="SpT1"){out<-c("RH1","RF1","RC1","RM1"); return(out)}
  else if(dataString=="SpT1pT2"){out<-c("RF1H2","RF1F2"); return(out)}
  else if(dataString=="T1"){out<-c("H1","F1","C1","M1"); return(out)}
  else if(dataString=="T1pT2"){out<-c("F1H2","F1F2"); return(out)}
  else if(dataString=="T2"){out<-c("H2","F2"); return(out)}
  else if(dataString=="T2fT1"){out<-c("F1H2","F1F2"); return(out)}
}

SE.Function <- function(x,plot_data){
  SD <- sd(plot_data[,x])/sqrt(length(plot_data[,x]))
  return(SD)
}

LWC.Function <- function(w,x,y,z){
  
  list.A <- w
  list.B <- x
  N <- as.numeric(y)
  data.raw <- z
  
  n.A <- as.numeric(length(list.A))
  n.B <- as.numeric(length(list.B))
  range1<-n.B*N 
  
  df.long <- as.data.frame(matrix(NA,nrow=N*n.A*n.B,ncol=5))
  names(df.long) <- c("Mask","Hemisphere","Behavior","Subj","Value")
  
  # col 1 = Mask
  a<-1; aa<-range1
  for(b in list.A){
    df.long[a:aa,1] <- b
    a<-aa+1; aa<-aa+range1
  }
  
  # col 2 = Hemisphere
  if(grepl("_",list.A[1])==TRUE){
    a<-1; aa<-range1
    for(b in list.A){
      df.long[a:aa,2] <- gsub("_.*","",b)
      a<-aa+1; aa<-aa+range1
    }
  }else{
    df.long[,2] <- "NA"
  }
  
  # col 3 = Behavior
  a<-1; aa<-N
  while(aa <= dim(df.long)[1]){
    for(b in 1:length(list.B)){
      df.long[a:aa,3] <- rep(list.B[b],N)
      a<-aa+1; aa<-aa+N
    }
  }

  # col 4 = subj
  a<-1; aa<-N
  while(aa <= dim(df.long)[1]){
    df.long[a:aa,4] <- 1:N
    a<-aa+1; aa<-aa+N
  }
  
  # col 5 = data
  data.hold <- melt(data.raw)
  df.long[,5] <- data.hold[,3]
  
  #write
  return(df.long)
}

Mdata.Function <- function(x){
  
  #masks
  DF <- x
  ind.mask <- grep("Mask", DF[,1])
  assign("num.mask", length(ind.mask), envir = .GlobalEnv)
  # num.mask <- length(ind.mask)
  
  #subjects
  ind.subj <- grep("File", DF[,1])
  len.subj <- length(ind.subj)
  assign("num.subj", len.subj/num.mask, envir = .GlobalEnv)
  # num.subj <- len.subj/num.mask
  
  #betas
  ind.betas <- grep("+tlrc", DF[,1])
  len.betas <- length(ind.betas)
  assign("num.betas", (len.betas/num.mask)/num.subj, envir = .GlobalEnv)
  # num.betas <- (len.betas/num.mask)/num.subj
  
  # organize data
  ind.data <- matrix(as.numeric(as.character(DF[grep("+tlrc",DF[,1]),3])),ncol=num.betas,byrow=T)
  df.hold <- matrix(0,nrow=num.subj, ncol=num.mask*num.betas)
  for(i in 1:num.mask){
    df.hold[,(num.betas*i-(num.betas-1)):(num.betas*i)] <- ind.data[(num.subj*i-(num.subj-1)):(num.subj*i),1:num.betas]
  }
  colnames(df.hold) <- c(as.character(rep(DF[ind.mask,1],each=num.betas)))
  h.mask.names <- sub("Mask ","",DF[ind.mask,1])
  assign("mask.names",h.mask.names, envir = .GlobalEnv)
  
  return(df.hold)
}

PostHoc.Function <- function(a,b,x,y){
  
  beh1 <- a
  beh2 <- b
  funcDF <- x
  outString <- y
  funcOut <- NA
  
  colA <- 1; for(col in 1:num.mask){
    
    colB <- colA+1
    tOut <- t.test(funcDF[,colA],funcDF[,colB],paired = T)
    funcOut <- c(funcOut,colnames(funcDF)[colA])
    funcOut <- c(funcOut,capture.output(print(tOut)))
    
    tsig <- tOut$p.value
    if(tsig < 0.05){
      
      holdDF <- matrix(0,ncol=2,nrow=dim(funcDF)[1])
      colnames(holdDF) <- c(beh1,beh2)
      holdDF[,1:2] <- cbind(funcDF[,colA],funcDF[,colB])
      holdMask <- gsub("Mask ","",colnames(funcDF)[colA])
      holdString <- paste0(comp,"_",beh1,"-",beh2,"_",holdMask)
      write.table(holdDF,file=paste0(ns_outDir,"Table_",holdString,".txt"),sep="\t",row.names=F,col.names=T,quote=F)
    }
    colA <- colB+1
  }
  writeLines(funcOut,paste0(ns_outDir,outString))
}




###################
# NeuroSynth Analysis
###################
# j <- "Betas_T1_NS_data.txt"

NSmaster_list <- read.table(paste0(nsDir,"Master_list_NS.txt"))
fdrOut <- NA

for(j in t(NSmaster_list)){
  
  # get info
  hold <- read.delim(paste0(nsDir,j),header=F)
  Mdata <- Mdata.Function(hold)
  
  hold <- gsub("^.*?_","",j)
  comp <- gsub("_.*$","",hold)
  beh <- BehNames.Function(comp)
  
  ### RR1 update - hardcode HvM, CvF for SpT1, T1
  if(grepl("_T1_",j)==T || grepl("_SpT1_",j)==T){
    
    # Anova on H/M
    dfHM <- matrix(0,ncol=2*num.mask,nrow=num.subj)
    if(comp == "T1"){
      dfHM <- Mdata[,c(1,4,5,8,9,12,13,16,17,20,21,24,25,28,29,32,33,36,37,40)]
    }else{
      dfHM <- Mdata[,c(1,4,5,8,9,12,13,16)]
    }
      
    dfHM_long <- LWC.Function(mask.names,c("Hit","Miss"),num.subj,dfHM)
    statsHM <- ezANOVA(dfHM_long,dv=Value,wid=Subj,within=c(Mask,Behavior),type='III')
    outputHM <- capture.output(print(statsHM))
    writeLines(outputHM,paste0(ns_outDir,"Stats_AN_",comp,"_MxHM.txt"))

    # Post-hoc
    sigHM <- statsHM$`Sphericity Corrections`$`p[GG]`[2]
    if(sigHM < 0.05){
      fdrOut <- rbind(fdrOut,c(paste0(comp,"-HM"),sigHM))
      PostHoc.Function("Hit","Miss",dfHM,paste0("Stats_TT_HvM_",comp,".txt"))
    }
    
    # Anova on C/F
    dfCF <- matrix(0,ncol=2*num.mask,nrow=num.subj)
    if(comp == "T1"){
      dfCF <- Mdata[,c(3,2,7,6,11,10,15,14,19,18,23,22,27,26,31,30,35,34,39,38)]
    }else{
      dfCF <- Mdata[,c(3,2,7,6,11,10,15,14)]
    }
    dfCF_long <- LWC.Function(mask.names,c("CR","FA"),num.subj,dfCF)
    statsCF <- ezANOVA(dfCF_long,dv=Value,wid=Subj,within=c(Mask,Behavior),type='III')
    outputCF <- capture.output(print(statsCF))
    writeLines(outputCF,paste0(ns_outDir,"Stats_AN_",comp,"_MxCF.txt"))

    sigCF <- statsCF$`Sphericity Corrections`$`p[GG]`[2]
    if(sigCF < 0.05){
      fdrOut <- rbind(fdrOut,c(paste0(comp,"-CF"),sigCF))
      PostHoc.Function("CR","FA",dfCF,paste0("Stats_TT_CvF_",comp,".txt"))
    }
    
  }else{
    
    ### Convert to Long, run ANOVA 
    Mdata_long <- LWC.Function(mask.names,beh,num.subj,Mdata)
    stats <- ezANOVA(Mdata_long,dv=Value,wid=Subj,within=c(Mask,Behavior),type='III')

    # write, patch for multple T1 tests
    if(grepl("NSEnc",j)==T){
      holdOut<-paste0(comp,"Enc")
    }else{
      holdOut<-comp
    }
    output <- capture.output(stats)
    writeLines(output,paste0(ns_outDir,"Stats_AN_",holdOut,".txt"))
    
    # post hoc
    pval <- stats$`Sphericity Corrections`$`p[GG]`[2]
    if(pval < 0.05){
      fdrOut <- rbind(fdrOut,c(comp,pval))
      PostHoc.Function(beh[1],beh[2],Mdata, paste0("Stats TT_",comp,".txt"))
    }
  }
}


### FDR correction
fdrDF <- as.data.frame(fdrOut[-1,],row.names = T)
fdrDF <- cbind(fdrDF,NA)
colnames(fdrDF) <- c("Comp","Orig","Corr")

p.input <- as.numeric(as.character(fdrDF[,2]))
fdrCor <- p.adjust(p.input, method = "fdr", n=dim(fdrDF)[1])
fdrDF[,3] <- fdrCor

write.table(fdrDF,file=paste0(ns_outDir,"Stats_FDR.txt"),sep="\t",row.names=F,col.names=T,quote=F)




### Make graphs
# Graph out whatever survives t-tests
tableList <- list.files(path = ns_outDir, pattern = "Table")

# j <- tableList[1]
for(j in tableList){
  
  tmp <- gsub("Table_","",j)
  holdStr <- gsub(".txt","",tmp)
  comp <- gsub("_.*$","",holdStr)
  tmp2 <- gsub("^.*?_","",holdStr)
  roi <- gsub("^.*?_","",tmp2)
  tmp3 <- gsub("_.*$","",tmp2)
  behA <- gsub("-.*$","",tmp3)
  behB <- gsub("^.*-","",tmp3)
  
  df <- read.delim(paste0(ns_outDir,j))

  MEANS <- colMeans(df)
  E.BARS<-NA
  for(a in 1:dim(df)[2]){
    E.BARS[a] <- SE.Function(a,df)
  }
  RANGE <- range(c(MEANS,MEANS-E.BARS,MEANS+E.BARS,0))
  

  plotable <- matrix(0,nrow=2,ncol=2)
  plotable[1,] <- MEANS
  plotable[2,] <- E.BARS
  
  if(writeGraph == 1){
    graphOut <- paste0(ns_outDir,"Graph_",holdStr,".png")
    bitmap(graphOut, width = 6.5, units = 'in', res=300)
  }
  barCenters <- barplot(plotable[1,], names.arg = c(behA,behB), main=paste(comp,roi), ylab="Parameter Estimate",ylim=RANGE)
  segments(barCenters, MEANS-E.BARS, barCenters, MEANS+E.BARS)
  arrows(barCenters, MEANS-E.BARS, barCenters, MEANS+E.BARS, lwd = 1, angle = 90, code = 3, length = 0.05)
  set.pos <- rowMeans(plotable); if(set.pos[1]>0){POS<-3}else{POS<-1}
  text(barCenters,0,round(plotable[1,],4),cex=1,pos=POS,font=2)
  if(writeGraph == 1){
    dev.off()
  }
}
