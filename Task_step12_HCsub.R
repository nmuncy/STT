library(reshape2)
library(ez)




###--- Notes
#
# This script will conduct two-factor within-subject ANOVAs
# on HC subfield betas. Factor A = ROI, Factor B = behavior.
# Significant interactions will be investigated.
#
# FDR correction, graphs removed (no sig detected)




###################
# Set up
###################

# parDir <- "../Analyses/"
parDir <- "/run/user/1000/gvfs/afp-volume:host=10.5.72.28,user=exp,volume=Yorick/STT_reml/Analyses/"
roiDir <- paste0(parDir,"roiAnalysis/")
subDir <- paste0(roiDir,"sub_betas_RR1/")
hc_outDir <- paste0(roiDir,"sub_stats_RR1/")




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

PostHoc.Function <- function(x,y){
  
  funcDF <- x
  outString <- y
  funcOut <- NA
  
  colA <- 1; for(col in 1:num.mask){
    colB <- colA+1
    tOut <- t.test(funcDF[,colA],funcDF[,colB],paired = T)
    funcOut <- c(funcOut,colnames(funcDF)[colA])
    funcOut <- c(funcOut,capture.output(print(tOut)))
    colA <- colB+1
  }
  writeLines(funcOut,paste0(hc_outDir,outString))
}




###################
# HC Sub Analysis
###################

HCmaster_list <- read.table(paste0(subDir,"Master_list_Sub.txt"))

# j <- "Betas_SpT1_sub_data.txt"
for(j in t(HCmaster_list)){
  
  # get info
  hold <- read.delim(paste0(subDir,j),header=F)
  Mdata <- Mdata.Function(hold)

  hold <- gsub("^.*?_","",j)
  comp <- gsub("_.*$","",hold)
  beh <- BehNames.Function(comp)

  ### RR1 update
  # SpT1, T1 have 4 betas, while all others have 2
  # so hardcode HvM, CvF
  
  if(grepl("_T1_",j)==T || grepl("_SpT1_",j)==T){
    
    ## Anova on H/M
    # split, organize df
    dfHM <- matrix(0,ncol=2*num.mask,nrow=num.subj)
    dfHM <- Mdata[,c(1,4,5,8,9,12,13,16)]
    
    # convert to long, stats
    dfHM_long <- LWC.Function(mask.names,c("Hit","Miss"),num.subj,dfHM)
    statsHM <- ezANOVA(dfHM_long,dv=Value,wid=Subj,within=c(Mask,Behavior),type='III')
    
    # write
    outputHM <- capture.output(print(statsHM))
    writeLines(outputHM,paste0(hc_outDir,"Stats_AN_",comp,"_MxHM.txt"))

    # Post-hoc, if sig - could be a function
    sigHM <- statsHM$`Sphericity Corrections`$`p[GG]`[2]
    if(sigHM < 0.05){
      PostHoc.Function(dfHM,paste0("Stats_TT_HvM_",comp,".txt"))
    }
    
    # Anova on C/F
    dfCF <- matrix(0,ncol=2*num.mask,nrow=num.subj)
    dfCF <- Mdata[,c(3,2,7,6,11,10,15,14)]
    dfCF_long <- LWC.Function(mask.names,c("CR","FA"),num.subj,dfCF)
    statsCF <- ezANOVA(dfCF_long,dv=Value,wid=Subj,within=c(Mask,Behavior),type='III')
    outputCF <- capture.output(print(statsCF))
    writeLines(outputCF,paste0(hc_outDir,"Stats_AN_",comp,"_MxCF.txt"))
    
    sigCF <- statsCF$`Sphericity Corrections`$`p[GG]`[2]
    if(sigCF < 0.05){
      PostHoc.Function(dfCF,paste0("Stats_TT_CvF_",comp,".txt"))
    }
    
  }else{
    
    # omnibus anova
    Mdata_long <- LWC.Function(mask.names,beh,num.subj,Mdata)
    stats <- ezANOVA(Mdata_long,dv=Value,wid=Subj,within=c(Mask,Behavior),type='III')
    
    output <- capture.output(stats)
    writeLines(output,paste0(hc_outDir,"Stats_AN_",comp,".txt"))

    # post-hoc, if sig
    pgg <- stats$`Sphericity Corrections`$`p[GG]`[2]
    if(pgg < 0.05){
      PostHoc.Function(Mdata,paste0("Stats_TT_",comp,".txt"))
    }
  }
}


