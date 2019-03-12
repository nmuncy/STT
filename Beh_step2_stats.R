


### --- Notes
#
# This script will read in raw data from the excel sheet,
# determine behaviors, and perform stats




#--------------------------------#
#   Part Zero: Set up
#--------------------------------#

### load data and packages
library(stringr)
library("openxlsx")
library(tiff)
# library("extrafont")


### for testing
#subjnames <- "s3852"
# j<-"s1875"


### Variables
dataDir <- "/Volumes/Yorick/STT_reml/Analyses/behAnalysis/"
iStats <- paste0(dataDir,"data_indiv/")
gStats <- paste0(dataDir,"data_group/")

subjnames <- read.xlsx(paste0(dataDir,"STT_all_data.xlsx"),"Subjects")
subjnames <- t(subjnames)


### Toggles
part1 <- 1   # make df
part2 <- 1   # calc bins
part3 <- 1   # d', props
part4 <- 1   # stats
part5 <- 1   # figs



### Functions
Outlier.Function <- function(x){
  
  hold.list <- " "
  df.x <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])
  df.x <- x
  
  for(xx in 2:dim(df.x)[2]){
    
    assign(paste0("list.",colnames(df.x)[xx]),0,envir = .GlobalEnv)
    
    h.iqr <- IQR(df.x[,xx])
    h.quant <- quantile(as.numeric(df.x[,xx]),names=F)
    h.max <- (1.5*h.iqr)+h.quant[4]
    h.min <- h.quant[2]-(1.5*h.iqr)
    
    hold.list <- c(hold.list,colnames(df.x)[xx])
    hold.list <- c(hold.list,paste("max = ",h.max))
    hold.list <- c(hold.list,paste("min = ",h.min))
    
    for(yy in 1:dim(df.x)[1]){
      if(as.numeric(df.x[yy,xx]) > h.max || as.numeric(df.x[yy,xx]) < h.min){
        
        assign(paste0("list.",colnames(df.x)[xx]),c(get(paste0("list.",colnames(df.x)[xx])),df.x[yy,1]),envir = .GlobalEnv)
        hold.out <- paste0(df.x[yy,1]," : ",df.x[yy,xx])
        hold.list <- c(hold.list,hold.out)
      }
    }
    assign(paste0("list.",colnames(df.x)[xx]),get(paste0("list.",colnames(df.x)[xx]))[-1],envir = .GlobalEnv)
  }
  return(hold.list)
}

Prop.Function <- function(a,b){
  func.df <- matrix(0,nrow=length(a),ncol=2)
  func.df[,1] <- a
  func.df[,2] <- b
  input <- colMeans(func.df)
  output <- prop.test(input[1],input[2],p=0.5,correct=TRUE)
  return(output)
}




#--------------------------------#
#   Part One: Make data frames
#--------------------------------#
if(part1 == 1){
  
  for(j in subjnames){
    
    data <- read.xlsx(paste0(dataDir,"STT_all_data.xlsx"),as.character(j))
    
    ### 0) extract raw data
    # stimulus
    ind <- NA
    StimRow <- NA
    for(i in 1:dim(data)[1]){
      if(grepl("StimImage",as.character(data[i,]))==T){
        ind <- c(ind,i)
        StimRow <- rbind(StimRow, str_split(data[i,]," ")[[1]][2])
      }
    }
    
    # response time
    ind1 <- NA
    RTRow <- NA
    for(i in 1:dim(data)[1]){
      if(grepl(".RT:",as.character(data[i,]))==T){
        ind1 <- c(ind1,i)
        RTRow <- rbind(RTRow, str_split(data[i,]," ")[[1]][2])
      }
    }
    
    # response
    ind2 <- NA
    RespRaw <- NA
    for(i in 1:dim(data)[1]){
      if(grepl(".RESP:",as.character(data[i,]))==T){
        ind2 <- c(ind2,i)
        RespRaw <- rbind(RespRaw, str_split(data[i,]," ")[[1]][2])
      }
    }
    
    # clean RespCol
    Nth.keep <- function(RespRaw, n)RespRaw[(seq(n,to=nrow(RespRaw),by=n)),]
    RespRow <- Nth.keep(RespRaw,2)
    
    
    # transpose data
    StimCola <- t(StimRow)
    StimColb <- t(StimCola)
    RTCola <- t(RTRow)
    RTColb <- t(RTCola)
    RespCola <- t(RespRow)
    RespCol <- t(RespCola)
    
    # clean NAs from StimCol, RTCol
    StimColc <- StimColb[!is.na(StimColb)]
    RTCol <- RTColb[!is.na(RTColb)]
    
    # fix redundancies (from Test2) in StimCol, clip off front
    StimCold <- StimColc[-seq(434,864,2)]
    StimCol <- substr(StimCold,9,nchar(StimCold))
    
    
    
    ### 1) Organize by Trial
    
    # Study
    StudyStim <- StimCol[1:216]
    for(i in 1:216){
      StudyStim[i] <- str_split(StudyStim[i], '\\\\')[[1]][2]
    }
    StudyRT <- RTCol[1:216]
    StudyResp <- RespCol[1:216]
    Study_1 <- data.frame(StudyStim, StudyResp, StudyRT)
    colnames(Study_1) <- c("Stimulus","KeyPress","RT")
    
    
    # Test1
    Test1Stim <- StimCol[217:432]
    for(i in 1:216){
      Test1Stim[i] <- str_split(Test1Stim[i], '\\\\')[[1]][2]
    }
    Test1RT <- RTCol[217:432]
    Test1Resp <- RespCol[217:432]
    Test1_1 <- data.frame(Test1Stim, Test1Resp, Test1RT)
    colnames(Test1_1) <- c("Stimulus","KeyPress","RT")
    
    
    # Test2
    Test2Stim <- StimCol[433:648]
    Test2RT <- RTCol[433:648]
    Test2Resp <- RespCol[433:648]
    Test2_1 <- data.frame(Test2Stim, Test2Resp, Test2RT)
    colnames(Test2_1) <- c("Stimulus","KeyPress","RT")
    
    
    ### 2) encode data
    
    ## Study
    # 1 = In, 2 = Out
    
    Study_2 <- Study_1
    Study_2$Behavior <- NA
    Study_2$Image <- NA
    
    for(i in 1:216){
      if(Study_2[i,3]==0 || grepl(999, Study_2[i,1])==T){
        Study_2[i,4] <- 0 }
      else{
        if(Study_2[i,2]==1){
          Study_2[i,4] <- "In"}
        else{Study_2[i,4] <- "Out"}
      }
      if(grepl(999, Study_2[i,1])==T){
        Study_2[i,5] <- substr(as.character(Study_2[i,1]),1,3)}
      else{Study_2[i,5] <- substr(as.character(Study_2[i,1]),1,nchar(as.character(Study_2[i,1]))-5)}
    }
    
    # clean jitter stimuli
    Study_2 <- Study_2[!(Study_2$Image=="999"),]
    row.names(Study_2) <- 1:nrow(Study_2)
    Study <- Study_2
    
    
    ## Test1
    # a-1 = Hit; a-2 = Miss
    # b-1 = FA; b-2 = CR
    
    Test1_2 <- Test1_1
    Test1_2[,2] <- as.numeric(as.character(Test1_1[,2]))
    Test1_2$Behavior <- NA
    Test1_2$Image <- NA
    Test1_2$StimType <- NA
    
    for(i in 1:216){
      if(grepl('a',Test1_2[i,1])){
        Test1_2[i,6] <- "Targ"
      }else if(grepl('b',Test1_2[i,1])){
        Test1_2[i,6] <- "Lure"
      }else{
        Test1_2[i,6] <- "Jit"
      }  
      if(as.numeric(as.character(Test1_2[i,3]))==0|Test1_2[i,2]>2){
        Test1_2[i,2] <- "NR"
        Test1_2[i,4] <- "NR"
      }else{
        if(grepl('a',Test1_2[i,1])){
          if(Test1_2[i,2]==1){
            Test1_2[i,4] <- "Hit"}
          else{Test1_2[i,4] <- "Miss"}
        }
        else{
          if(Test1_2[i,2]==1){
            Test1_2[i,4] <- "FA"}
          else{Test1_2[i,4] <- "CR"}
        }
      }
      if(grepl(999, Test1_2[i,1])==T){
        Test1_2[i,5] <- substr(as.character(Test1_2[i,1]),1,3)}
      else{Test1_2[i,5] <- substr(as.character(Test1_2[i,1]),1,nchar(as.character(Test1_2[i,1]))-5)
      }
    }

    Test1_2 <- Test1_2[!(Test1_2$Image=="999"),]
    row.names(Test1_2) <- 1:nrow(Test1_2)
    T1 <- Test1_2
    
    
    ## Test2
    # a-1,2 = Hit; a-3,4 = FA
    # b-1,2 = FA; b-3,4 = Hit
    
    Test2_2 <- Test2_1
    Test2_2[,2] <- as.numeric(as.character(Test2_1[,2]))
    Test2_2$Behavior <- NA
    Test2_2$Conf <- NA
    Test2_2$Image <- NA
    
    for(i in 1:216){
      if(Test2_2[i,3]==0|as.numeric(as.character(Test2_2[i,2]))>5 || grepl(999, Test2_2[i,1])==T){
        Test2_2[i,4] <- Test2_2[i,5] <- Test2_2[i,2] <- "NR"
      }else{
        if(grepl("a",Test2_2[i,1])){
          if((Test2_2[i,2]==1|Test2_2[i,2]==2)){
            Test2_2[i,4] <- "Hit"
            Test2_2[i,5] <- as.numeric(as.character(Test2_2[i,2]))}
          else{
            Test2_2[i,4] <- "FA"
            Test2_2[i,5] <- as.numeric(as.character(Test2_2[i,2]))
          }
        }
        else{
          if(Test2_2[i,2]==3|Test2_2[i,2]==4){
            Test2_2[i,4] <- "Hit"
            Test2_2[i,5] <- abs(as.numeric(as.character(Test2_2[i,2]))-5)}
          else{
            Test2_2[i,4] <- "FA"
            Test2_2[i,5] <- abs(as.numeric(as.character(Test2_2[i,2]))-5)
          }
        }
      }
      if(grepl(999, Test2_2[i,1])==T){
        Test2_2[i,6] <- substr(as.character(Test2_2[i,1]),1,3)}
      else{Test2_2[i,6] <- substr(as.character(Test2_2[i,1]),1,nchar(as.character(Test2_2[i,1]))-5)}
    }
    
    Test2_2 <- Test2_2[!(Test2_2$Image=="999"),]
    row.names(Test2_2) <- 1:nrow(Test2_2)

    
    ## parse T2 CI info
    # 2,3 = Somewhat Conf
    # 1,4 = Very Conf
    
    Test2_CI_2 <- Test2_2
    Test2_CI_2$VerySome <- NA
    
    for(i in 1:200){
      if(Test2_CI_2[i,5]==1|Test2_CI_2[i,5]==4){Test2_CI_2[i,7] <- "V"}
      else if(Test2_CI_2[i,5]==2|Test2_CI_2[i,5]==3){Test2_CI_2[i,7] <- "S"}
      else{Test2_CI_2[i,7] <- "NR"}
    }
    T2 <- Test2_CI_2
    
    
    ## T1pT2 - match up in T1 time
    T1pT2 <- matrix(NA,nrow=200,ncol=6)
    colnames(T1pT2) <- colnames(Test1_2)
    
    for(i in c(1:3,5:6)){
      T1pT2[,i] <- Test1_2[,i]
    }
    for(i in 1:200){
      hold1 <- as.character(Test1_2[i,5],1,nchar(as.character(Test1_2[i,5])))
      position1 <- which(Test2_2[,6]==as.character((hold1),arr.ind=T))
      T1pT2[i,4] <- paste0(Test1_2[i,4], "-", Test2_2[position1,4])
    }
    
    
    ## T1pT2CI - match up in T1 time
    T1pT2CI <- matrix(NA,nrow = 200,ncol=6)
    colnames(T1pT2CI) <- colnames(Test1_2)
    
    for(i in c(1:3,5:6)){
      T1pT2CI[,i] <- Test1_2[,i]
    }
    for(i in 1:200){
      hold1 <- as.character(Test1_2[i,5],1,nchar(as.character(Test1_2[i,5])))
      position1 <- which(Test2_CI_2[,6]==as.character((hold1),arr.ind=T))
      T1pT2CI[i,4] <- paste0(Test1_2[i,4], "-", Test2_CI_2[position1,4], "-", Test2_CI_2[position1,7])
    }
    
    
    ## T2fT1 - match in T2 time, include confidence
    T2fT1 <- matrix(NA,nrow=200,ncol=7)
    colnames(T2fT1) <- colnames(Test2_CI_2)
    
    for(i in c(1:3,5:7)){
      if(i == 1){
        T2fT1[,i] <- as.character(Test2_CI_2[,i])
      }else{
        T2fT1[,i] <- Test2_CI_2[,i]
      }
    }
    
    for(i in 1:200){
      hold1 <- as.character(Test2_2[i,6],1,nchar(as.character(Test2_2[i,6])))
      position1 <- which(Test1_2[,5]==as.character((hold1),arr.ind=T))
      T2fT1[i,4] <- paste0(Test1_2[position1,4], "-", Test2_2[i,4])
      # T2r1_DF_3_a[i,2] <- as.numeric(as.character(Test2_2[i,3]))
    }
    
    
    ## Store data as txt files
    # "Study"
    writeList <- c("T1", "T2","T2fT1", "T1pT2","T1pT2CI")
    for(i in writeList){
      write.table(get(i), paste0(iStats,j, "_DF_", i, ".txt"), row.names = F, quote = FALSE, sep = "\t")
    }
  }
}




#--------------------------------#
#   Part Two: Determine bin counts
#--------------------------------#
if(part2 == 1){
  
  ### Set up
  # make empty data frames
  
  #T1
  df.T1 <- matrix(0,nrow=0,ncol=6)
  colnames(df.T1) <- c("Subj","Hit","FA","CR","Miss","NR")
  
  #T2
  df.T2 <- matrix(0,nrow=0,ncol=4)
  colnames(df.T2) <- c("Subj","Hit","FA","NR")
  
  #T2CI
  df.T2CI <- matrix(0,nrow=0,ncol=10)
  h.list<-"Subj"
  for(k in c("Hit","FA","NR")){
    for(m in c("V","S","NR")){
      h.list <- c(h.list,paste0(k,"-",m))
    }
  }
  colnames(df.T2CI) <- h.list
  
  #T1pT2, T2fT1
  df.T1pT2 <- matrix(0,nrow=0,ncol=16)
  df.T2fT1 <- matrix(0,nrow=0,ncol=16)
  h.list <- "Subj"
  for(k in c("Hit","FA","Miss","CR","NR")){
    for(m in c("Hit","FA","NR")){
      h.list <- c(h.list,paste0(k,"-",m))
    }
  }
  colnames(df.T1pT2) <- h.list
  colnames(df.T2fT1) <- h.list
  
  #T1pT2CI
  df.T1pT2CI <- matrix(0,nrow=0,ncol=46)
  h.list <- "Subj"
  for(k in c("Hit","FA","Miss","CR","NR")){
    for(m in c("Hit","FA","NR")){
      for(n in c("V","S","NR")){
        h.list <- c(h.list,paste0(k,"-",m,"-",n))
      }
    }
  }
  colnames(df.T1pT2CI) <- h.list
  
  
  ### Determine counts for e/participant
  for(j in subjnames){
    
    ## Test1
    # set indices for e/behavior, determine where in DF behavior occurs,
    # find length of behavior index, augment master DF
    
    ind.Hit <- ind.FA <- ind.Miss <- ind.CR <- ind.NR <-0
    df.raw <- read.delim(paste0(iStats,j,"_DF_T1.txt"))
    for(i in 1:200){
      for(k in c("Hit","FA","Miss","CR","NR")){
        if(grepl(k,as.character(df.raw[i,4]))==T){
          assign(paste0("ind.",k),c(get(paste0("ind.",k)),i))
        } 
      }
    }
    for(k in c("Hit","FA","Miss","CR","NR")){
      hold <- as.numeric(length(get(paste0("ind.",k))))-1
      assign(paste0("out.",k),hold)
    }
    df.T1 <- rbind(df.T1,c(j,out.Hit,out.FA,out.CR,out.Miss,out.NR))

    
    ## Test2
    ind.Hit <- ind.FA <- ind.NR <- 0
    df.raw <- read.delim(paste0(iStats,j,"_DF_T2.txt"))
    for(i in 1:200){
      for(k in c("Hit","FA","NR")){
        if(grepl(k,as.character(df.raw[i,4]))==T){
          assign(paste0("ind.",k),c(get(paste0("ind.",k)),i))
        } 
      }
    }
    for(k in c("Hit","FA","NR")){
      hold <- as.numeric(length(get(paste0("ind.",k))))-1
      assign(paste0("out.",k),hold)
    }
    df.T2 <- rbind(df.T2,c(j,out.Hit,out.FA,out.NR))

    
    ## Test2CI
    for(k in c("Hit","FA","NR")){
      for(m in c("V","S","NR")){
        assign(paste0("ind.",k,"-",m),0)
      }
    }
    df.raw <- read.delim(paste0(iStats,j,"_DF_T2.txt"))
    for(i in 1:200){
      for(k in c("Hit","FA","NR")){
        for(m in c("V","S","NR")){
          if(grepl(k,as.character(df.raw[i,4]))==T && grepl(m,as.character(df.raw[i,7]))==T){
           assign(paste0("ind.",k,"-",m),c(get(paste0("ind.",k,"-",m)),i))
          }
        }
      }
    }
    for(k in c("Hit","FA","NR")){
      for(m in c("V","S","NR")){
        hold <- as.numeric(length(get(paste0("ind.",k,"-",m))))-1
        assign(paste0("out.",k,"-",m),hold)
      }
    }
    out <- j
    for(k in c("Hit","FA","NR")){
      for(m in c("V","S","NR")){
        hold <- as.numeric(length(get(paste0("ind.",k,"-",m))))-1
        assign(paste0("ind.",k,"-",m),hold)
        out <- c(out,get(paste0("ind.",k,"-",m)))
      }
    }
    df.T2CI <- rbind(df.T2CI,out)
    
    
    ## T1pT2
    for(k in c("Hit","FA","Miss","CR","NR")){
      for(m in c("Hit","FA","NR")){
        assign(paste0("ind.",k,"-",m),0)
      }
    }
    df.raw <- read.delim(paste0(iStats,j,"_DF_T1pT2.txt"))
    for(i in 1:200){
      for(k in c("Hit","FA","Miss","CR","NR")){
        for(m in c("Hit","FA","NR")){
          if(grepl(paste0(k,"-",m),as.character(df.raw[i,4]))==T){
            assign(paste0("ind.",k,"-",m),c(get(paste0("ind.",k,"-",m)),i))
          }
        }
      }
    }
    out <- j
    for(k in c("Hit","FA","Miss","CR","NR")){
      for(m in c("Hit","FA","NR")){
        hold <- as.numeric(length(get(paste0("ind.",k,"-",m))))-1
        assign(paste0("ind.",k,"-",m),hold)
        out <- c(out,get(paste0("ind.",k,"-",m)))
      }
    }
    df.T1pT2 <- rbind(df.T1pT2,out)

    
    ## T1pT2CI
    for(k in c("Hit","FA","Miss","CR","NR")){
      for(m in c("Hit","FA","NR")){
        for(n in c("V","S","NR")){
          assign(paste0("ind.",k,"-",m,"-",n),0)
        }
      }
    }
    df.raw <- read.delim(paste0(iStats,j,"_DF_T1pT2CI.txt"))
    for(i in 1:200){
      for(k in c("Hit","FA","Miss","CR","NR")){
        for(m in c("Hit","FA","NR")){
          for(n in c("V","S","NR")){
            if(grepl(paste0(k,"-",m,"-",n),as.character(df.raw[i,4]))==T){
              assign(paste0("ind.",k,"-",m,"-",n),c(get(paste0("ind.",k,"-",m,"-",n)),i))
            }
          }
        }
      }
    }
    out <- j
    for(k in c("Hit","FA","Miss","CR","NR")){
      for(m in c("Hit","FA","NR")){
        for(n in c("V","S","NR")){
          hold <- as.numeric(length(get(paste0("ind.",k,"-",m,"-",n))))-1
          assign(paste0("ind.",k,"-",m,"-",n),hold)
          out <- c(out,get(paste0("ind.",k,"-",m,"-",n)))
        }
      }
    }
    df.T1pT2CI <- rbind(df.T1pT2CI,out)
    
    
    ## T2fT1 (not doing CI)
    for(k in c("Hit","FA","Miss","CR","NR")){
      for(m in c("Hit","FA","NR")){
        assign(paste0("ind.",k,"-",m),0)
      }
    }
    df.raw <- read.delim(paste0(iStats,j,"_DF_T2fT1.txt"))
    for(i in 1:200){
      for(k in c("Hit","FA","Miss","CR","NR")){
        for(m in c("Hit","FA","NR")){
          if(grepl(paste0(k,"-",m),as.character(df.raw[i,4]))==T){
            assign(paste0("ind.",k,"-",m),c(get(paste0("ind.",k,"-",m)),i))
          }
        }
      }
    }
    out <- j
    for(k in c("Hit","FA","Miss","CR","NR")){
      for(m in c("Hit","FA","NR")){
        hold <- as.numeric(length(get(paste0("ind.",k,"-",m))))-1
        assign(paste0("ind.",k,"-",m),hold)
        out <- c(out,get(paste0("ind.",k,"-",m)))
      }
    }
    df.T2fT1 <- rbind(df.T2fT1,out)
  }
  
  
  ### write out
  write.table(df.T1,paste0(gStats,"Group_table_respBin_T1.txt"), row.names = F, quote = FALSE, sep = "\t")
  write.table(df.T2,paste0(gStats,"Group_table_respBin_T2.txt"), row.names = F, quote = FALSE, sep = "\t")
  write.table(df.T2CI,paste0(gStats,"Group_table_respBin_T2CI.txt"), row.names = F, quote = FALSE, sep = "\t")
  write.table(df.T1pT2,paste0(gStats,"Group_table_respBin_T1pT2.txt"), row.names = F, quote = FALSE, sep = "\t")
  write.table(df.T1pT2CI,paste0(gStats,"Group_table_respBin_T1pT2CI.txt"), row.names = F, quote = FALSE, sep = "\t")
  write.table(df.T2fT1,paste0(gStats,"Group_table_respBin_T2fT1.txt"), row.names = F, quote = FALSE, sep = "\t")
}




#--------------------------------#
#   Part Three: d's, proportions
#--------------------------------#
if(part3 == 1){
  
  ### T1 d', prop
  df.T1 <- read.delim(paste0(gStats,"Group_table_respBin_T1.txt"))
  
  df.T1d <- matrix(0,nrow=dim(df.T1)[1],ncol=2)
  colnames(df.T1d) <- c("Subj","d'")
  df.T1p <- matrix(0,nrow=dim(df.T1)[1],ncol=5)
  colnames(df.T1p) <- c("Subj","pHit","nTR","pCR","nLR")
  
  for(i in 1:dim(df.T1)[1]){
    
    n.TargR <- df.T1[i,2]+df.T1[i,5]
    n.LureR <- df.T1[i,3]+df.T1[i,4]
    pHit <- df.T1[i,2]/n.TargR
    pFA <- df.T1[i,3]/n.LureR
    zHit <- qnorm(pHit)
    zFA <- qnorm(pFA)
    
    d.T1out <- zHit - zFA
    df.T1d[i,] <- c(as.character(df.T1[i,1]),d.T1out)
    df.T1p[i,] <- c(as.character(df.T1[i,1]),pHit,n.TargR,1-pFA,n.LureR)
  }
  write.table(df.T1d,paste0(gStats,"Group_table_dprime_T1.txt"), row.names = F, quote = FALSE, sep = "\t")
  write.table(df.T1p,paste0(gStats,"Group_table_prop_T1.txt"), row.names = F, quote = FALSE, sep = "\t")

  
  ### T2 d', prop
  df.T2 <- read.delim(paste0(gStats,"Group_table_respBin_T2.txt"))
  
  df.T2d <- matrix(0,nrow=dim(df.T2)[1],ncol=2)
  colnames(df.T2d) <- c("Subj","d'")
  df.T2p <- matrix(0,nrow=dim(df.T2)[1],ncol=3)
  colnames(df.T2p) <- c("Subj","pHit","nSR")
  
  for(i in 1:dim(df.T2)[1]){

    n.StimR <- df.T2[i,2]+df.T2[i,3]
    pHit <- df.T2[i,2]/n.StimR
    pFA <- df.T2[i,3]/n.StimR
    zHit <- qnorm(pHit)
    zFA <- qnorm(pFA)

    d.T2out <- (zHit - zFA)/sqrt(2)
    df.T2d[i,] <- c(as.character(df.T2[i,1]),d.T2out)
    df.T2p[i,] <- c(as.character(df.T2[i,1]),pHit,n.StimR)
  }
  write.table(df.T2d,paste0(gStats,"Group_table_dprime_T2.txt"), row.names = F, quote = FALSE, sep = "\t")
  write.table(df.T2p,paste0(gStats,"Group_table_prop_T2.txt"), row.names = F, quote = FALSE, sep = "\t")
  
  
  ### T2CI prop
  df.T2CI <- read.delim(paste0(gStats,"Group_table_respBin_T2CI.txt"))
  df.T2CIp <- matrix(0,nrow=dim(df.T2CI)[1],ncol=5)
  colnames(df.T2CIp) <- c("Subj","pHV","nHit","pFV","nFA")
  
  for(i in 1:dim(df.T2CI)[1]){
    
    n.Hit <- df.T2CI[i,2]+df.T2CI[i,3]
    n.FA <- df.T2CI[i,5]+df.T2CI[i,6]
    pHV <- df.T2CI[i,2]/n.Hit
    pFV <- df.T2CI[i,5]/n.FA
    df.T2CIp[i,] <- c(as.character(df.T2CI[i,1]),pHV,n.Hit,pFV,n.FA)
  }
  write.table(df.T2CIp,paste0(gStats,"Group_table_prop_T2CI.txt"), row.names = F, quote = FALSE, sep = "\t")
  
  
  ### T1pT2 prop
  df.T1pT2 <- read.delim(paste0(gStats,"Group_table_respBin_T1pT2.txt"))
  df.T1pT2p <- matrix(0,nrow=dim(df.T1pT2)[1],ncol=5)
  colnames(df.T1pT2p) <- c("Subj","pHH","nHpR","pFH","nFpR")
  
  for(i in 1:dim(df.T1pT2)[1]){
    
    n.HpR <- df.T1pT2[i,2]+df.T1pT2[i,3]
    n.FpR <- df.T1pT2[i,5]+df.T1pT2[i,6]
    pHH <- df.T1pT2[i,2]/n.HpR
    pFH <- df.T1pT2[i,5]/n.FpR
    df.T1pT2p[i,] <- c(as.character(df.T1pT2[i,1]),pHH,n.HpR,pFH,n.FpR)
  }
  write.table(df.T1pT2p,paste0(gStats,"Group_table_prop_T1pT2.txt"), row.names = F, quote = FALSE, sep = "\t")
  
  
  ### T1pT2CI prop
  df.T1pT2CI <- read.delim(paste0(gStats,"Group_table_respBin_T1pT2CI.txt"))
  df.T1pT2CIp <- matrix(0,nrow=dim(df.T1pT2CI)[1],ncol=5)
  colnames(df.T1pT2CIp) <- c("Subj","pFHV","nFHC","pFFV","nFFC")
  
  for(i in 1:dim(df.T1pT2CI)[1]){
    
    n.FHC <- df.T1pT2CI[i,11]+df.T1pT2CI[i,12]
    n.FFC <- df.T1pT2CI[i,14]+df.T1pT2CI[i,15]
    pFHV <- df.T1pT2CI[i,11]/n.FHC
    pFFV <- df.T1pT2CI[i,14]/n.FFC
    df.T1pT2CIp[i,] <- c(as.character(df.T1pT2CI[i,1]),pFHV,n.FHC,pFFV,n.FFC)
  }
  write.table(df.T1pT2CIp,paste0(gStats,"Group_table_prop_T1pT2CI.txt"), row.names = F, quote = FALSE, sep = "\t")
}




#--------------------------------#
#   Part Four: Stats
#--------------------------------#
if(part4 == 1){
  
  ### d' ttests
  
  ## T1vT2, T1v0, T2v0
  hold1 <- read.delim(paste0(gStats,"Group_table_dprime_T1.txt"))
  hold2 <- read.delim(paste0(gStats,"Group_table_dprime_T2.txt"))
  
  df.hold <- matrix(0,nrow=35,ncol=2)
  colnames(df.hold) <- c("T1d","T2d")
  df.hold[,1] <- hold1[,2]
  df.hold[,2] <- hold2[,2]
  
  hold.mean <- colMeans(df.hold)
  hold.1v2 <- t.test(df.hold[,1],df.hold[,2],paired=TRUE)
  hold.1v0 <- t.test(df.hold[,1],mu=0,paired=FALSE)
  hold.2v0 <- t.test(df.hold[,2],mu=0,paired=FALSE)
  
  # write it
  output <- "D' ttests"
  output <- c(output,capture.output(print(hold.mean)))
  output <- c(output,"T1vT2")
  output <- c(output, capture.output(print(hold.1v2)))
  output <- c(output,"T1v0")
  output <- c(output, capture.output(print(hold.1v0)))
  output <- c(output,"T2v0")
  output <- c(output, capture.output(print(hold.2v0)))
  writeLines(output,paste0(gStats,"Stats_ttest_dprimes.txt"))
  
  # capture, adjust p-values
  p.1v2 <- hold.1v2$p.value
  p.1v0 <- hold.1v0$p.value
  p.2v0 <- hold.2v0$p.value

  p.input <- c(p.1v2,p.1v0,p.2v0)
  p.output <- p.adjust(p.input, method="fdr", n=length(p.input))

  output <- matrix(NA,nrow=length(p.output),ncol=2)
  colnames(output) <- c("orig","adj")
  rownames(output) <- c("1v2","1v0","2v0")
  output[,1] <- p.input
  output[,2] <- p.output
  write.table(output,file=paste0(gStats,"Stats_ttest_dprimes_corrected.txt"),sep="\t",row.names=T,col.names=T)
  

  
  ## d' outliers
  df.out <- cbind("Subj" = as.character(hold1[,1]),df.hold)
  outL.d <- Outlier.Function(df.out)
  writeLines(outL.d,paste0(gStats,"Stats_outliers_dprimes.txt"))
  
  
  # d' ttest no outliers
  df.xo <- cbind("Subj" = as.character(hold1[,1]),df.hold)
  for(i in 2:dim(df.xo)[2]){
    assign(paste0("v.",colnames(df.xo)[i]),NA)
    h.name <- paste0("list.",colnames(df.xo)[i])
    for(j in 1:dim(df.xo)[1]){
      if(df.xo[j,1] %in% get(h.name) != TRUE){
        assign(paste0("v.",colnames(df.xo)[i]),c(get(paste0("v.",colnames(df.xo)[i])),as.numeric(df.xo[j,i])))
      }
    }
    assign(paste0("v.",colnames(df.xo)[i]),get(paste0("v.",colnames(df.xo)[i]))[-1])
  }
  
  hold.mx1 <- mean(v.T1d)
  hold.mx2 <- mean(v.T2d)
  hold.1v2 <- t.test(v.T1d,v.T2d,paired=TRUE)
  hold.1v0 <- t.test(v.T1d,mu=0,paired=FALSE)
  hold.2v0 <- t.test(v.T2d,mu=0,paired=FALSE)
  
  # write it
  output <- "D' ttests"
  output <- c(output,hold.mx1,hold.mx2)
  output <- c(output,"T1vT2")
  output <- c(output, capture.output(print(hold.1v2)))
  output <- c(output,"T1v0")
  output <- c(output, capture.output(print(hold.1v0)))
  output <- c(output,"T2v0")
  output <- c(output, capture.output(print(hold.2v0)))
  writeLines(output,paste0(gStats,"Stats_ttest_dprimes_xout.txt"))
  
  # capture, adjust p-values
  p.1v2 <- hold.1v2$p.value
  p.1v0 <- hold.1v0$p.value
  p.2v0 <- hold.2v0$p.value
  
  p.input <- c(p.1v2,p.1v0,p.2v0)
  p.output <- p.adjust(p.input, method="fdr", n=length(p.input))
  
  output <- matrix(NA,nrow=length(p.output),ncol=2)
  colnames(output) <- c("orig","adj")
  rownames(output) <- c("1v2","1v0","2v0")
  output[,1] <- p.input
  output[,2] <- p.output
  write.table(output,file=paste0(gStats,"Stats_ttest_dprimes_xout_corrected.txt"),sep="\t",row.names=T,col.names=T)
  
  
  
  ### prop testing - one sample prop testing w/cont corrections

  # T1 - pHit, pCR
  hold1 <- read.delim(paste0(gStats,"Group_table_respBin_T1.txt"))
  t1H <- Prop.Function(hold1[,2],(hold1[,2]+hold1[,5]))
  t1C <- Prop.Function(hold1[,4],(hold1[,3]+hold1[,4]))
  
  # T2 - pHit
  hold2 <- read.delim(paste0(gStats,"Group_table_respBin_T2.txt"))
  t2H <- Prop.Function(hold2[,2],(hold2[,2]+hold2[,3]))
  
  # T2CI - pHV, pFV
  hold3 <- read.delim(paste0(gStats,"Group_table_respBin_T2CI.txt"))
  t2HV <- Prop.Function(hold3[,2],(hold3[,2]+hold3[,3]))
  t2FV <- Prop.Function(hold3[,5],(hold3[,5]+hold3[,6]))
  
  # T1pT2 - pHH, pFH
  hold4 <- read.delim(paste0(gStats,"Group_table_respBin_T1pT2.txt"))
  t1HH <- Prop.Function(hold4[,2],(hold4[,2]+hold4[,3]))
  t1FH <- Prop.Function(hold4[,5],(hold4[,5]+hold4[,6]))
  
  # T1pT2CI - pFHV, pFFV
  hold5 <- read.delim(paste0(gStats,"Group_table_respBin_T1pT2CI.txt"))
  t1FHV <- Prop.Function(hold5[,11],(hold5[,11]+hold5[,12]))
  t1FFV <- Prop.Function(hold5[,14],(hold5[,14]+hold5[,15]))

  # write it
  output <- "Prop testing"
  output <- c(output,"T1 pHit")
  output <- c(output,capture.output(print(t1H)))
  output <- c(output,"T1 pCR")
  output <- c(output,capture.output(print(t1C)))
  output <- c(output,"T2 pHit")
  output <- c(output,capture.output(print(t2H)))
  output <- c(output,"T2CI pHV")
  output <- c(output,capture.output(print(t2HV)))
  output <- c(output,"T2CI pFV")
  output <- c(output,capture.output(print(t2FV)))
  output <- c(output,"T1pT2 pHH")
  output <- c(output,capture.output(print(t1HH)))
  output <- c(output,"T1pT2 pFH")
  output <- c(output,capture.output(print(t1FH)))
  output <- c(output,"T1pT2CI pFHV")
  output <- c(output,capture.output(print(t1FHV)))
  output <- c(output,"T1pT2CI pFFV")
  output <- c(output,capture.output(print(t1FFV)))
  writeLines(output,paste0(gStats,"Stats_prop_testing.txt"))
  
  # capture, adjust p-values
  p.t1H <- t1H$p.value
  p.t1C <- t1C$p.value
  p.t2H <- t2H$p.value
  p.t2HV <- t2HV$p.value
  p.t2FV <- t2FV$p.value
  p.t1HH <- t1HH$p.value
  p.t1FH <- t1FH$p.value
  p.t1FHV <- t1FHV$p.value
  p.t1FFV <- t1FFV$p.value

  p.input <- c(p.t1H,p.t1C,p.t2H,p.t2HV,p.t2FV,p.t1HH,p.t1FH,p.t1FHV,p.t1FFV)
  p.output <- p.adjust(p.input, method="fdr", n=length(p.input))
  
  output <- matrix(NA,nrow=length(p.output),ncol=2)
  colnames(output) <- c("orig","adj")
  rownames(output) <- c("t1H", "t1C", "t2H", "t2HV", "t2FV", "t1HH", "t1FH", "t1FHV", "t1FFV")
  output[,1] <- p.input
  output[,2] <- p.output
  write.table(output,file=paste0(gStats,"Stats_prop_testing_corrected.txt"),sep="\t",row.names=T,col.names=T)

  
  
  ## prop outliers
  hold11 <- read.delim(paste0(gStats,"Group_table_prop_T1.txt"))
  df.out <- cbind("Subj" = as.character(hold11[,1]),"T1Hit" = hold11[,2], "T1CR" = hold11[,4])
  outL.T1 <- Outlier.Function(df.out)
  
  hold21 <- read.delim(paste0(gStats,"Group_table_prop_T2.txt"))
  df.out <- cbind("Subj" = as.character(hold21[,1]),"T2Hit" = hold21[,2])
  outL.T2 <- Outlier.Function(df.out)
  
  hold31 <- read.delim(paste0(gStats,"Group_table_prop_T2CI.txt"))
  df.out <- cbind("Subj" = as.character(hold31[,1]),"T2HV" = hold31[,2], "T2FV" = hold31[,4])
  outL.T2CI <- Outlier.Function(df.out)
  
  hold41 <- read.delim(paste0(gStats,"Group_table_prop_T1pT2.txt"))
  df.out <- cbind("Subj" = as.character(hold41[,1]),"T1HH" = hold41[,2], "T1FH" = hold41[,4])
  outL.T1pT2 <- Outlier.Function(df.out)
  
  hold51 <- read.delim(paste0(gStats,"Group_table_prop_T1pT2CI.txt"))
  df.out <- cbind("Subj" = as.character(hold51[,1]),"T1FHV" = hold51[,2], "T1FFV" = hold51[,4])
  outL.T1pT2CI <- Outlier.Function(df.out)
  
  # write it
  output <- "Prop testing"
  output <- c(output,outL.T1)
  output <- c(output,outL.T2)
  output <- c(output,outL.T2CI)
  output <- c(output,outL.T1pT2)
  output <- c(output,outL.T1pT2CI)
  writeLines(output,paste0(gStats,"Stats_outliers_props.txt"))
  
  
  ### prop test, no out
  # T1 - pHit, pCR
  df.t1xo <- matrix(0,nrow=0,ncol=2)
  for(i in 1:dim(hold1)[1]){
    if(as.character(hold1[i,1]) %in% list.T1Hit != TRUE){
      df.t1xo <- rbind(df.t1xo,c(hold1[i,2],hold1[i,2]+hold1[i,5]))
    }
  }
  t1Hxo <- Prop.Function(df.t1xo[,1],df.t1xo[,2])
  
  df.t1cxo <- matrix(0,nrow=0,ncol=2)
  for(i in 1:dim(hold1)[1]){
    if(as.character(hold1[i,1]) %in% list.T1CR != TRUE){
      df.t1cxo <- rbind(df.t1cxo,c(hold1[i,4],hold1[i,3]+hold1[i,4]))
    }
  }
  t1Cxo <- Prop.Function(df.t1cxo[,1],df.t1cxo[,2])
  
  
  #T2 - pHit
  df.t2hxo <- matrix(0,nrow=0,ncol=2)
  for(i in 1:dim(hold2)[1]){
    if(as.character(hold2[i,1]) %in% list.T2Hit != TRUE){
      df.t2hxo <- rbind(df.t2hxo,c(hold2[i,2],hold2[i,2]+hold2[i,3]))
    }
  }
  t2Hxo <- Prop.Function(df.t2hxo[,1],df.t2hxo[,2])

  
  # T2CI - pHV, pFV
  df.t2hcxo <- matrix(0,nrow=0,ncol=2)
  for(i in 1:dim(hold3)[1]){
    if(as.character(hold3[i,1]) %in% list.T2HV != TRUE){
      df.t2hcxo <- rbind(df.t2hcxo,c(hold3[i,2],hold3[i,2]+hold3[i,3]))
    }
  }
  t2HCxo <- Prop.Function(df.t2hcxo[,1],df.t2hcxo[,2])
  
  df.t2fcxo <- matrix(0,nrow=0,ncol=2)
  for(i in 1:dim(hold3)[1]){
    if(as.character(hold3[i,1]) %in% list.T2FV != TRUE){
      df.t2fcxo <- rbind(df.t2fcxo,c(hold3[i,5],hold3[i,5]+hold3[i,6]))
    }
  }
  t2FCxo <- Prop.Function(df.t2fcxo[,1],df.t2fcxo[,2])


  # T1pT2 - pHH, pFH
  df.t1hhxo <- matrix(0,nrow=0,ncol=2)
  for(i in 1:dim(hold4)[1]){
    if(as.character(hold4[i,1]) %in% list.T1HH != TRUE){
      df.t1hhxo <- rbind(df.t1hhxo,c(hold4[i,2],hold4[i,2]+hold4[i,3]))
    }
  }
  t1HHxo <- Prop.Function(df.t1hhxo[,1],df.t1hhxo[,2])
  
  df.t1fhxo <- matrix(0,nrow=0,ncol=2)
  for(i in 1:dim(hold4)[1]){
    if(as.character(hold4[i,1]) %in% list.T1FH != TRUE){
      df.t1fhxo <- rbind(df.t1fhxo,c(hold4[i,5],hold4[i,5]+hold4[i,6]))
    }
  }
  t1FHxo <- Prop.Function(df.t1fhxo[,1],df.t1fhxo[,2])
  
  
  # T1pT2CI - pFHV, pFFV
  df.t1fhvxo <- matrix(0,nrow=0,ncol=2)
  for(i in 1:dim(hold5)[1]){
    if(as.character(hold5[i,1]) %in% list.T1FHV != TRUE){
      df.t1fhvxo <- rbind(df.t1fhvxo,c(hold5[i,11],hold5[i,11]+hold5[i,12]))
    }
  }
  t1FHVxo <- Prop.Function(df.t1fhvxo[,1],df.t1fhvxo[,2])
  
  df.t1ffvxo <- matrix(0,nrow=0,ncol=2)
  for(i in 1:dim(hold5)[1]){
    if(as.character(hold5[i,1]) %in% list.T1FFV != TRUE){
      df.t1ffvxo <- rbind(df.t1ffvxo,c(hold5[i,14],hold5[i,14]+hold5[i,15]))
    }
  }
  t1FFVxo <- Prop.Function(df.t1ffvxo[,1],df.t1ffvxo[,2])

  
  # write it
  output <- "Prop testing"
  output <- c(output,"T1 pHit")
  output <- c(output,capture.output(print(t1Hxo)))
  output <- c(output,"T1 pCR")
  output <- c(output,capture.output(print(t1Cxo)))
  output <- c(output,"T2 pHit")
  output <- c(output,capture.output(print(t2Hxo)))
  output <- c(output,"T2CI pHV")
  output <- c(output,capture.output(print(t2HCxo)))
  output <- c(output,"T2CI pFV")
  output <- c(output,capture.output(print(t2FCxo)))
  output <- c(output,"T1pT2 pHH")
  output <- c(output,capture.output(print(t1HHxo)))
  output <- c(output,"T1pT2 pFH")
  output <- c(output,capture.output(print(t1FHxo)))
  output <- c(output,"T1pT2CI pFHV")
  output <- c(output,capture.output(print(t1FHVxo)))
  output <- c(output,"T1pT2CI pFFV")
  output <- c(output,capture.output(print(t1FFVxo)))
  writeLines(output,paste0(gStats,"Stats_prop_testing_xout.txt"))
  
  
  # capture, adjust p-values
  p.t1H <- t1Hxo$p.value
  p.t1C <- t1Cxo$p.value
  p.t2H <- t2Hxo$p.value
  p.t2HV <- t2HCxo$p.value
  p.t2FV <- t2FCxo$p.value
  p.t1HH <- t1HHxo$p.value
  p.t1FH <- t1FHxo$p.value
  p.t1FHV <- t1FHVxo$p.value
  p.t1FFV <- t1FFVxo$p.value
  
  p.input <- c(p.t1H,p.t1C,p.t2H,p.t2HV,p.t2FV,p.t1HH,p.t1FH,p.t1FHV,p.t1FFV)
  p.output <- p.adjust(p.input, method="fdr", n=length(p.input))
  
  output <- matrix(NA,nrow=length(p.output),ncol=2)
  colnames(output) <- c("orig","adj")
  rownames(output) <- c("t1H", "t1C", "t2H", "t2HV", "t2FV", "t1HH", "t1FH", "t1FHV", "t1FFV")
  output[,1] <- p.input
  output[,2] <- p.output
  write.table(output,file=paste0(gStats,"Stats_prop_testing_xout_corrected.txt"),sep="\t",row.names=T,col.names=T)
  
  
  
  
  ### Ttest prop bin
  
  # T1pT2CI - pFHV v pFFV
  hold.df <- matrix(NA,nrow=35,ncol=2)
  colnames(hold.df) <- c("FHV","FFV")
  hold.df[,1] <- hold5[,11]
  hold.df[,2] <- hold5[,14]
  hold.FHVvFFV <- t.test(hold.df[,1],hold.df[,2],paired=TRUE)
  hold.means <- colMeans(hold.df)
  
  output <- hold.means
  output <- c(output,capture.output(print(hold.FHVvFFV)))
  writeLines(output,paste0(gStats,"Stats_ttest_T1pT2CI_FHVvsFFV.txt"))
}




#--------------------------------#
#   Part Five: Figures
#--------------------------------#
if(part5 == 1){
  
  ### proportions
  df.hold1 <- read.delim(paste0(gStats,"Group_table_prop_T1.txt"))
  df.hold2 <- read.delim(paste0(gStats,"Group_table_prop_T2.txt"))
  df.hold3 <- read.delim(paste0(gStats,"Group_table_prop_T2CI.txt"))
  df.hold4 <- read.delim(paste0(gStats,"Group_table_prop_T1pT2.txt"))
  
  df.graph <- matrix(NA,nrow=35,ncol=7)
  df.graph[,1] <- df.hold1[,2]  #pT1H
  df.graph[,2] <- df.hold1[,4]  #pT1C 
  df.graph[,3] <- df.hold2[,2]  #pT2H 
  df.graph[,4] <- df.hold3[,2]  #pT2HV 
  df.graph[,5] <- df.hold3[,4]  #pT2FV 
  df.graph[,6] <- df.hold4[,2]  #pHHF 
  df.graph[,7] <- df.hold4[,4]  #pFHF 
  colnames(df.graph) <- c("T1 Hit Rate", "T1 CR Rate", "T2 Hit Rate", "T2 Hit-Very Rate", "T2 FA-Very Rate", "T1 Hit p T2 Hit Rate", "T1 FA p T2 Hit Rate")
  XLAB <- c("T1 Hit Rate", "T1 CR Rate", "T2 Hit Rate", "T2 Hit-Very Rate", "T2 FA-Very Rate", "T1 Hit p T2 Hit Rate", "T1 FA p T2 Hit Rate")

  tiff(paste0(gStats,"Fig_props.tiff"), height = 5.5, width = 7, units = 'in', res=300)
  par(mar=c(8,5,4,2),family="Times New Roman")
  hold.graph <- boxplot(df.graph, at=c(1,2, 4,5,6, 8,9) ,ylim=c(0,1), xaxt="n", xlab="", ylab="Proportion", col="white", cex.lab=1.5)
  axis(1,at=c(1,2, 4,5,6, 8,9),labels=F)
  text(x = c(1,2, 4,5,6, 8,9), y = -0.075, srt = 45, adj = 1, labels = XLAB, xpd = T)
  mtext("Comparison",side=1,line=6.5,cex=1.5)
  abline(h = 0.5)
  title(main=list("Behavioral Responses", cex=2))
  dev.off()
  
  
  # ### d'
  # df.hold1 <- read.delim(paste0(gStats,"Group_table_dprime_T1.txt"))
  # df.hold2 <- read.delim(paste0(gStats,"Group_table_dprime_T2.txt"))
  # 
  # df.graph <- matrix(NA,nrow=35,ncol=2)
  # df.graph[,1] <- df.hold1[,2]
  # df.graph[,2] <- df.hold2[,2]
  # # colnames(df.graph) <- c("T1d'","T2d'")
  # XLAB <- c("T1d'","T2d'")
  # 
  # # tiff(paste0(gStats,"Fig_props.tiff"), height = 5.5, width = 7, units = 'in', res=300)
  # # par(mar=c(8,5,4,2),family="Times New Roman")
  # hold.graph <- boxplot(df.graph, ylim=c(-0.5,1.5), xlab=XLAB, ylab="Proportion", col="white", cex.lab=2)
  # abline(h = 0)
  # title(main=list("Behavioral Responses", cex=2))
  # # dev.off()
}
















