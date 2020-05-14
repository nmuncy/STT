library(reshape2)
library(ez)



### --- Notes
#
# This script will run stats on ETAC betas



###################
# Set up
###################

# parDir <- "../Analyses/"
parDir <- "/run/user/1000/gvfs/afp-volume:host=10.5.72.28,user=exp,volume=Yorick/STT_reml/Analyses/"

doWrite <- 1
doGraphs <- 1

etacDir <- paste0(parDir,"grpAnalysis/")
etacData_location <- paste0(etacDir,"etac_betas/")
etacData_list <- read.table(paste0(etacData_location,"All_List.txt"))
etac_outDir <- etacData_location
etac_statsDir <- paste0(etacDir,"etac_stats_RR1/")




###################
# Functions
###################

GraphEtacNames.Function <- function(dataString){
  if(dataString=="SpT1"){return(list(n1="RpHit", n2="RpFA"))}
  else if(dataString=="SpT1pT2"){return(list(n1="RpFpH", n2="RpFpF"))}
  else if(dataString=="T1"){return(list(n1="Hit", n2="FA"))}
  else if(dataString=="T1pT2"){return(list(n1="FpH", n2="FpF"))}
  else if(dataString=="T2"){return(list(n1="Hit", n2="FA"))}
  else if(dataString=="T2fT1"){return(list(n1="HfH", n2="HfF"))}
}

BehNames.Function <- function(dataString){
  if(dataString=="SpT1"){out<-c("RH1","RF1","RC1","RM1"); return(out)}
  else if(dataString=="SpT1pT2"){out<-c("RF1H2","RF1F2"); return(out)}
  else if(dataString=="T1"){out<-c("H1","F1","C1","M1"); return(out)}
  else if(dataString=="T1pT2"){out<-c("F1H2","F1F2"); return(out)}
  else if(dataString=="T2"){out<-c("H2","F2"); return(out)}
  else if(dataString=="T2fT1"){out<-c("F1H2","F1F2"); return(out)}
}

EtacNames.Function <- function(x,y){
  if(grepl("SpT1",x)==T){
    if(y=="RVS"){
      return("R. Visual Stream")
    }else if(y=="LVS"){
      return("L. Visual Stream")
    }else if(y=="LPPCU"){
      return("L. Posterior Precuneus")
    }else if(y=="LAG"){
      return("L. Angular Gyrus")
    }else if(y=="RIPS"){
      return("R. Intraparietal Sulcus")
    }else if(y=="LDMPFC"){
      return("L. Dorsal Medial PFC")
    }else if(y=="RAG"){
      return("R. Angular Gyrus")
    }
  }else if(grepl("T1_",x)==T){
    if(y=="RMPFC"){
      return("R. Medial PFC")
    }else if(y=="RAG"){
      return("R. Angular Gyrus")
    }else if(y=="RAMTG"){
      return("R. Anterior Middle Temporal Gyrus")
    }else if(y=="RIPS"){
      return("R. Intraparietal Sulcus")
    }else if(y=="RPCU"){
      return("R. Precuneus")
    }else if(y=="LINS"){
      return("L. Insular Cortex")
    }else if(y=="RIFG"){
      return("R. Inferior Frontal Gyrus")
    }else if(y=="RPUT"){
      return("R. Putamen")
    }else if(y=="LDMPFC"){
      return("L. Dorsal Medial PFC")
    }
  }else if(grepl("T1pT2",x)==T){
    if(y=="LINS"){
      return("L. Insular Cortex")
    }else if(y=="RSMG"){
      return("R. Supramarginal Gyrus")
    }
  }else if(x == "T2fT1"){
    if(y=="RPCU"){
      return("R. Precuneus")
    }else if(y=="RAG"){
      return("R. Angular Gyrus")
    }
  }
}

SE.Function <- function(x,plot_data){
  SD <- sd(plot_data[,x])/sqrt(length(plot_data[,x]))
  return(SD)
}

Graph.Function <- function(DF,output_name,maskN,out_place){
  
  TITLE <- maskN
  MEANS <- colMeans(DF)
  E.BARS<-NA
  for(a in 1:dim(DF)[2]){
    E.BARS[a] <- SE.Function(a,DF)
  }

  RANGE <- range(c(MEANS,MEANS-E.BARS,MEANS+E.BARS,0))
  ROI <- EtacNames.Function(output_name,maskN)
  XNAMES <- GraphEtacNames.Function(output_name)

  plotable <- matrix(0,nrow=2,ncol=num.betas)
  plotable[1,] <- MEANS
  plotable[2,] <- E.BARS
  
  if(doWrite == 1){
    graphOut <- paste0(out_place,"Graph_",output_name,"_",TITLE,".png")
    # bitmap(graphOut, width = 6.5, units = 'in', type="tiff24nc", res=1200)
    bitmap(graphOut, width = 6.5, units = 'in', res=300)
  }
  barCenters <- barplot(plotable[1,], names.arg = c(XNAMES), main=ROI, ylab="Beta Coefficient",ylim=RANGE)
  segments(barCenters, MEANS-E.BARS, barCenters, MEANS+E.BARS)
  arrows(barCenters, MEANS-E.BARS, barCenters, MEANS+E.BARS, lwd = 1, angle = 90, code = 3, length = 0.05)
  set.pos <- rowMeans(plotable); if(set.pos[1]>0){POS<-3}else{POS<-1}
  text(barCenters,0,round(plotable[1,],4),cex=1,pos=POS,font=2)
  if(doWrite == 1){
    dev.off()
  }
}

TT.Function <- function(x, y, dataN, maskN, out_place){
  ttest_out <- t.test(x,y,paired=T)
  meanX <- mean(x)
  meanY <- mean(y)
  if(doWrite == 1){
    output <- c(meanX, meanY)
    output <- c(output, capture.output(print(ttest_out)))
    writeLines(output,paste0(out_place,"Stats_TT_",dataN,"_",maskN,".txt"))
    return(ttest_out)
  }else{
    print(ttest_out)
  }
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

MkTable.Function <- function(x,y){
  
  DF <- x
  DF.perm <- y
  hold.post <- matrix(NA,nrow=1,ncol=1+(num.comp*6))
  hold.post[1,1] <- hold.mask
  
  d<-2; for(k in 1:dim(DF.perm)[1]){
    
    colA <- DF.perm[k,1]; colB <- DF.perm[k,2]
    t.hold <- TT.Function(DF[,colA],DF[,colB],comp,paste0(hold.mask,"_",beh[colA],"-",beh[colB]),ns_outDir)

    t.hold.comp <- paste0(beh[colA],"-",beh[colB])
    t.hold.t <- as.numeric(t.hold$statistic)
    t.hold.df <- as.numeric(t.hold$parameter)
    t.hold.p <- as.numeric(t.hold$p.value)
    t.hold.lb <- as.numeric(t.hold$conf.int[1])
    t.hold.ub <- as.numeric(t.hold$conf.int[2])

    t.hold.capture <- c(t.hold.comp,t.hold.t,t.hold.df,t.hold.p,t.hold.lb,t.hold.ub)
    dd <- d+5
    hold.post[1,d:dd]<-t.hold.capture
    d <- dd+1
  }
  return(hold.post)
}




###################
# ETAC
###################
# # For testing
# i <- "All_Betas_SpT1.txt"
# j <- "Betas_SpT1_LAG.txt"

for(i in t(etacData_list)){
  
  beta_list <- read.table(paste0(etacData_location,i))
  
  df.post <- matrix(NA,nrow = dim(beta_list)[1],ncol = 6)
  colnames(df.post) <- c("ROI","T","df","p","LB","RB")
  
  count<-1
  for(j in t(beta_list)){
    
    ### Get, clean data
    raw_data <- read.delim2(paste0(etacData_location,j),header=F)
    
    # num subjects
    ind.subj <- grep("File", raw_data[,1])
    num.subj <- as.numeric(length(ind.subj))
    
    # num betas - will be 2 for ETAC
    ind.beta <- grep("beh", raw_data[,2])
    num.betas <- as.numeric(length(ind.beta)/num.subj)
    
    ind.beh1 <- ind.beta[c(TRUE, FALSE)]
    ind.beh2 <- ind.beta[c(FALSE, TRUE)]
    
    # fill df
    df <- matrix(0,ncol=num.betas,nrow=num.subj)
    df[,1] <- as.numeric(as.character(raw_data[ind.beh1,3]))
    df[,2] <- as.numeric(as.character(raw_data[ind.beh2,3]))
    
    # write out
    hold<-gsub("^.*_","",i); comp<-gsub("\\..*","",hold)
    hold1<-gsub("^.*_","",j); anat<-gsub("\\..*","",hold1)
    
    if(doWrite == 1){
      write.table(df,paste0(etacData_location,"Avg_Betas_",comp,"_",anat,".txt"),col.names=F, row.names=F,sep = "\t")
    }
    
    
    ### Stats
    t.out <- TT.Function(df[,1],df[,2],comp,anat,etac_statsDir)
    
    hold.t <- round(t.out$statistic,digits=2)
    hold.df <- t.out$parameter
    hold.p <- round(t.out$p.value,digits=6)
    hold.ci <- round(t.out$conf.int,digits=3)
    hold.write <- c(anat,hold.t,hold.df,hold.p,hold.ci)
    df.post[count,] <- hold.write

    
    ### Graphs
    if(doGraphs == 1){
      Graph.Function(df,comp,anat,etac_statsDir)
    }
    count <- count+1
  }
  
  if(doWrite==1){
    write.table(df.post,paste0(etac_statsDir,"Table_",comp,".txt"),row.names = F, quote = F, sep = "\t")
  }
}

