library(tiff)


outDir <- "/Volumes/Yorick/STT_reml/write/figures/"
gStats <- "/Volumes/Yorick/STT_reml/Analyses/behAnalysis/data_group/"
h.colors <- c("darkred","darkblue","darkgray")

SE.Function <- function(plot_data){
  hold.xx <- NA
  for(xx in 1:dim(plot_data)[2]){
    hold.xx <- c(hold.xx, sd(plot_data[,xx])/sqrt(length(plot_data[,xx])))
  }
  hold.xx <- hold.xx[-1]
  return(hold.xx)
}



### Beh Graphs

# prop graph
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

tiff(paste0(outDir,"Fig_beh1.tiff"), height = 5.5, width = 7, units = 'in', res=300)
par(mar=c(8,5,4,2),family="Times New Roman")
hold.graph <- boxplot(df.graph, at=c(1,2, 4,5,6, 8,9) ,ylim=c(0,1), xaxt="n", xlab="", ylab="Proportion", col="white", cex.lab=1.5)
axis(1,at=c(1,2, 4,5,6, 8,9),labels=F)
text(x = c(1,2, 4,5,6, 8,9), y = -0.075, srt = 45, adj = 1, labels = XLAB, xpd = T)
mtext("Comparison",side=1,line=6.5,cex=1.5)
abline(h = 0.5)
title(main=list("Behavioral Responses", cex=2))
dev.off()




### NS graphs
draw_txt1 <- "*"; draw_txt2 <- "**"; draw_txt3 <- "***"

# T1 - read in, reorganize data frame
df.hold <- read.table("/Volumes/Yorick/STT_reml/Analyses/roiAnalysis/ns_stats/Betas_Table_T1.txt",header=T, sep = "\t")
df.hold <- as.matrix(df.hold)

df.reorganize <- matrix(NA,nrow=dim(df.hold)[1],ncol=dim(df.hold)[2])
df.reorganize[,1:3] <- df.hold[,1:3]
df.reorganize[,4:6] <- df.hold[,7:9]
df.reorganize[,7:9] <- df.hold[,10:12]
df.reorganize[,10:12] <- df.hold[,4:6]

# make plotable matrix
df.plotable <- matrix(NA,nrow=2,ncol=dim(df.hold)[2])
df.plotable[1,] <- colMeans(df.reorganize)
df.plotable[2,] <- SE.Function(df.reorganize)

# determine max/min values for plotting range
max_position <- which(df.plotable == max(df.plotable[1,]), arr.ind = T)
hold_max <- df.plotable[max_position] + df.plotable[2,max_position[2]]
graph_max <- hold_max+(0.1*hold_max)
min_position <- which(df.plotable == min(df.plotable[1,]), arr.ind=T)
hold_min <- df.plotable[min_position] - df.plotable[2,min_position[2]]
graph_min <- hold_min+(0.1*hold_min)

# save out as tiff
tiff(paste0(outDir,"Fig_Plot_NS_T1.tiff"), height = 5.5, width = 7, units = 'in', res=300)

# draw graph - ylim max was hardcoded for the legend issue
par(mar=c(5,5,4,2),family="Times New Roman")
hold.graph <- barplot(matrix(df.plotable[1,],nrow=3), beside=T, ylim=c(graph_min,0.6),ylab="Beta Estimate",cex.lab=1.5,legend=c("Hit","FA","CR"),col=h.colors)
abline(h=0)
segments(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,])
arrows(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
mtext(c("LMTL","LVMPFC","LDMPFC","LMFG"), side=1, at=hold.graph[1,]+0.5, cex=1.2)
title(main=list("NeuroSynth: Test1", cex=2))

# build matrix of sig levels, for determining how many stars to use
# row 1 = Hit-FA; row 2 = FA-CR; row 3 = Hit-CR
# 1 = *; 2 = **; 3 = ***
sig_matrix <- matrix(0,nrow=dim(hold.graph)[1],ncol=dim(hold.graph)[2])
sig_matrix[2:3,] <- 3
sig_matrix[3,1] <- 2 #LMTL
sig_matrix[2,1] <- 0 #LMTL
sig_matrix[2,2] <- 2 #LVMPFC

# determine location of bars, build list of every third (for 3 bars on graph)
hold_sequence <- 1:dim(df.plotable)[2]
hold_list <- hold_sequence[seq(1,length(hold_sequence),3)]

# for lists of length 3
count<-1; while(count <= dim(hold.graph)[2]){
  
  # get subset list, will permute through it. determine if any values are positive
  perm_hold <- c(hold_list[count], hold_list[count]+1, hold_list[count]+2)
  pos_value<-F; hold_start<-1
  for(x in perm_hold){
    if(df.plotable[1,x] > 0){
      pos_value<-T
    }
  }
  
  # loop through contrasts (permutations)
  c<-1; for(a in perm_hold){
    
    # set vars
    if(c < length(perm_hold)){
      start<-a; end<-a+1
      begin<-c; finish<-c+1
      h_star<-hold.graph[begin,count]+0.5
      inc<-0.1
    }else{
      start<-a-2; end<-a
      begin<-c-2; finish<-c
      h_star<-hold.graph[begin+1,count]
      inc<-0.25
    }
    
    # if bars need to be drawn
    if(sig_matrix[c,count] != 0){
      
      # determine num of stars
      hold_text <- sig_matrix[c,count]
      draw_star <- get(paste0("draw_txt",hold_text))
      
      # draw differently, depending on pos/neg values
      if(pos_value == T){
        larger <- which(df.plotable == max(df.plotable[1,start:end]), arr.ind=T)
        segment_y <- df.plotable[larger]+df.plotable[2,larger[2]]+inc*graph_max
        segments(hold.graph[begin,count], segment_y,hold.graph[finish,count],segment_y,lwd=2, xpd=T)
        text(x=h_star, y=segment_y+0.05*graph_max,draw_star, xpd=T) 
      }else{
        segment_y <- hold_start*0.1*graph_max
        segments(hold.graph[begin,count], segment_y,hold.graph[finish,count],segment_y,lwd=2)
        text(x=h_star, y=segment_y+0.05*graph_max,draw_star) 
        hold_start<-hold_start+2
      }
    }
    c<-c+1
  }
  count<-count+1
}

dev.off()


# T1pT2 - this one needs some tweaking with the sig bars
df.hold <- read.table("/Volumes/Yorick/STT_reml/Analyses/roiAnalysis/ns_stats/Betas_Table_T1pT2.txt",header=T, sep = "\t")
df.hold <- as.matrix(df.hold)

df.reorganize <- matrix(NA,nrow=dim(df.hold)[1],ncol=dim(df.hold)[2])
df.reorganize[,1:3] <- df.hold[,7:9]
df.reorganize[,4:9] <- df.hold[,1:6]

df.plotable <- matrix(NA,nrow=2,ncol=dim(df.hold)[2])
df.plotable[1,] <- colMeans(df.reorganize)
df.plotable[2,] <- SE.Function(df.reorganize)

max_position <- which(df.plotable == max(df.plotable[1,]), arr.ind = T)
hold_max <- df.plotable[max_position] + df.plotable[2,max_position[2]]
graph_max <- hold_max+(0.1*hold_max)
min_position <- which(df.plotable == min(df.plotable[1,]), arr.ind=T)
hold_min <- df.plotable[min_position] - df.plotable[2,min_position[2]]
graph_min <- hold_min+(0.1*hold_min)

tiff(paste0(outDir,"Fig_Plot_NS_T1pT2.tiff"), height = 7.5, width = 9, units = 'in', res=300)

par(mar=c(5,5,5,2),family="Times New Roman")
hold.graph <- barplot(matrix(df.plotable[1,],nrow=3), beside=T, ylim=c(graph_min,0.45),ylab="Beta Estimate",cex.lab=1.5,legend=c("FpH","FpF","CpH"),col=h.colors)
abline(h=0)
segments(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,])
arrows(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
mtext(c("LDMPFC","LAG","LMFG"), side=1, at=hold.graph[1,]+0.5, cex=1.2)
title(main=list("NeuroSynth: Test1 FAs preceding Test2", cex=2))

# row 1 = FpH-FpF; row 2 = FpF-CpH; row 3 = FpH-CpH
sig_matrix <- matrix(0,nrow=dim(hold.graph)[1],ncol=dim(hold.graph)[2])
sig_matrix[1,1] <- 1 #LDMPFC
sig_matrix[2:3,] <- 3
sig_matrix[2,3] <- 2 #LMFG
sig_matrix[3,3] <- 2 #LMFG

hold_sequence <- 1:dim(df.plotable)[2]
hold_list <- hold_sequence[seq(1,length(hold_sequence),3)]

count<-1; while(count <= dim(hold.graph)[2]){
  
  perm_hold <- c(hold_list[count], hold_list[count]+1, hold_list[count]+2)
  pos_value<-F; hold_start<-1
  for(x in perm_hold){
    if(df.plotable[1,x] > 0){
      pos_value<-T
    }
  }
  
  c<-1; for(a in perm_hold){
    
    if(c < length(perm_hold)){
      start<-a; end<-a+1
      begin<-c; finish<-c+1
      h_star<-hold.graph[begin,count]+0.5
      inc<-0.05
    }else{
      start<-a-2; end<-a
      begin<-c-2; finish<-c
      h_star<-hold.graph[begin+1,count]
      inc<-0.1
    }
    
    if(sig_matrix[c,count] != 0){
      
      hold_text <- sig_matrix[c,count]
      draw_star <- get(paste0("draw_txt",hold_text))
      
      if(pos_value == T){
        larger <- which(df.plotable == max(df.plotable[1,start:end]), arr.ind=T)
        segment_y <- df.plotable[larger]+df.plotable[2,larger[2]]+inc*graph_max
        segments(hold.graph[begin,count], segment_y,hold.graph[finish,count],segment_y,lwd=2, xpd=T)
        text(x=h_star, y=segment_y+0.02*graph_max,draw_star, xpd=T) 
      }else{
        segment_y <- hold_start*0.05*graph_max
        segments(hold.graph[begin,count], segment_y,hold.graph[finish,count],segment_y,lwd=2)
        text(x=h_star, y=segment_y+0.02*graph_max,draw_star) 
        hold_start<-hold_start+2
      }
    }
    c<-c+1
  }
  count<-count+1
}

dev.off()


# T2fT1
df.hold <- read.table("/Volumes/Yorick/STT_reml/Analyses/roiAnalysis/ns_stats/Betas_Table_T2fT1.txt",header=T, sep = "\t")
df.hold <- as.matrix(df.hold)

df.reorganize <- matrix(NA,nrow=dim(df.hold)[1],ncol=dim(df.hold)[2])
df.reorganize[,1:18] <- df.hold[,1:18]

df.plotable <- matrix(NA,nrow=2,ncol=dim(df.hold)[2])
df.plotable[1,] <- colMeans(df.reorganize)
df.plotable[2,] <- SE.Function(df.reorganize)

max_position <- which(df.plotable == max(df.plotable[1,]), arr.ind = T)
hold_max <- df.plotable[max_position] + df.plotable[2,max_position[2]]
graph_max <- hold_max+(0.1*hold_max)
min_position <- which(df.plotable == min(df.plotable[1,]), arr.ind=T)
hold_min <- df.plotable[min_position] - df.plotable[2,min_position[2]]
graph_min <- hold_min+(0.1*hold_min)

tiff(paste0(outDir,"Fig_Plot_NS_T2fT1.tiff"), height = 7.5, width = 9, units = 'in', res=300)

par(mar=c(5,5,4,2),family="Times New Roman")
hold.graph <- barplot(matrix(df.plotable[1,],nrow=3), beside=T, ylim=c(graph_min,graph_max),ylab="Beta Estimate",cex.lab=1.5,legend=c("HfF","FfF","HfC"),col=h.colors)
abline(h=0)
segments(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,])
arrows(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
mtext(c("LAG","LDMPFC","LMFG","LPCU","LTPJ","RPOS"), side=1, at=hold.graph[1,]+0.5, cex=1.2)
title(main=list("NeuroSynth: Test2 following Test1 FAs", cex=2))

# row 1 = HfF-FfF; row 2 = FfF-HfC; row 3 = HfF-HfC
sig_matrix <- matrix(0,nrow=dim(hold.graph)[1],ncol=dim(hold.graph)[2])
sig_matrix[2:3,]<-3
sig_matrix[2,1]<-2 #LAG

hold_sequence <- 1:dim(df.plotable)[2]
hold_list <- hold_sequence[seq(1,length(hold_sequence),3)]

count<-1; while(count <= dim(hold.graph)[2]){
  
  perm_hold <- c(hold_list[count], hold_list[count]+1, hold_list[count]+2)
  pos_value<-F; hold_start<-1
  for(x in perm_hold){
    if(df.plotable[1,x] > 0){
      pos_value<-T
    }
  }
  
  c<-1; for(a in perm_hold){
    
    if(c < length(perm_hold)){
      start<-a; end<-a+1
      begin<-c; finish<-c+1
      h_star<-hold.graph[begin,count]+0.5
      inc<-0.05
    }else{
      start<-a-2; end<-a
      begin<-c-2; finish<-c
      h_star<-hold.graph[begin+1,count]
      inc<-0.15
    }
    
    if(sig_matrix[c,count] != 0){
      
      hold_text <- sig_matrix[c,count]
      draw_star <- get(paste0("draw_txt",hold_text))
      
      if(pos_value == T){
        larger <- which(df.plotable == max(df.plotable[1,start:end]), arr.ind=T)
        segment_y <- df.plotable[larger]+df.plotable[2,larger[2]]+inc*graph_max
        segments(hold.graph[begin,count], segment_y,hold.graph[finish,count],segment_y,lwd=2, xpd=T)
        text(x=h_star, y=segment_y+0.04*graph_max,draw_star, xpd=T) 
      }else{
        segment_y <- hold_start*0.05*graph_max
        segments(hold.graph[begin,count], segment_y,hold.graph[finish,count],segment_y,lwd=2)
        text(x=h_star, y=segment_y+0.04*graph_max,draw_star) 
        hold_start<-hold_start+2
      }
    }
    c<-c+1
  }
  count<-count+1
}

dev.off()


# # Combined T1pT2, T2fT1
# df.hold1 <- read.table("/Volumes/Yorick/STT_reml/Analyses/roiAnalysis/ns_stats/Betas_Table_T1pT2.txt",header=T, sep = "\t")
# df.hold1 <- as.matrix(df.hold1)
# df.reorganize1 <- matrix(NA,nrow=dim(df.hold1)[1],ncol=dim(df.hold1)[2])
# df.reorganize1[,1:3] <- df.hold1[,7:9]
# df.reorganize1[,4:9] <- df.hold1[,1:6]
# df.plotable1 <- matrix(NA,nrow=2,ncol=dim(df.hold1)[2])
# df.plotable1[1,] <- colMeans(df.reorganize1)
# df.plotable1[2,] <- SE.Function(df.reorganize1)
# 
# df.hold2 <- read.table("/Volumes/Yorick/STT_reml/Analyses/roiAnalysis/ns_stats/Betas_Table_T2fT1.txt",header=T, sep = "\t")
# df.hold2 <- as.matrix(df.hold2)
# df.reorganize2 <- matrix(NA,nrow=dim(df.hold2)[1],ncol=dim(df.hold2)[2])
# df.reorganize2[,1:18] <- df.hold2[,1:18]
# df.plotable2 <- matrix(NA,nrow=2,ncol=dim(df.hold2)[2])
# df.plotable2[1,] <- colMeans(df.reorganize2)
# df.plotable2[2,] <- SE.Function(df.reorganize2)
#  
# 
# tiff(paste0(outDir,"Fig_Plot_NS_Combined.tiff"), height = 9, width = 9, units = 'in', res=300)
# par(mfrow=c(2,1),family="Times New Roman")
# 
# par(mar=c(1,5,3,2))
# hold.graph1 <- barplot(matrix(df.plotable1[1,],nrow=3), beside=T, ylim=c(-0.2,0.4),ylab="Beta Estimate",cex.lab=1.5,legend=c("FpH","FpF","CpH"),col=h.colors)
# abline(h=0)
# segments(hold.graph1, df.plotable1[1,]-df.plotable1[2,], hold.graph1, df.plotable1[1,]+df.plotable1[2,])
# arrows(hold.graph1, df.plotable1[1,]-df.plotable1[2,], hold.graph1, df.plotable1[1,]+df.plotable1[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
# text(x=hold.graph1[2,],y=-0.18,labels=c("LDMPFC","LAG","LMFG"),cex=1.2)
# title(main=list("Test1 preceding Test2", cex=2))
# 
# par(mar=c(1,5,2,2))
# hold.graph2 <- barplot(matrix(df.plotable2[1,],nrow=3), beside=T, ylim=c(-0.7,0.3),ylab="Beta Estimate",cex.lab=1.5,legend=c("HfF","FfF","HfC"),col=h.colors)
# abline(h=0)
# segments(hold.graph2, df.plotable2[1,]-df.plotable2[2,], hold.graph2, df.plotable2[1,]+df.plotable2[2,])
# arrows(hold.graph2, df.plotable2[1,]-df.plotable2[2,], hold.graph2, df.plotable2[1,]+df.plotable2[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
# text(x=hold.graph2[2,],y=-0.62,labels=c("LAG","LDMPFC","LMFG","LPCU","LTPJ","RPOS"),cex=1.2)
# title(main=list("Test2 following Test1", cex=2),line=-1)
# dev.off()    




### ETAC graphs
etacDir <- "/Volumes/Yorick/STT_reml/Analyses/grpAnalysis/etac_betas/"

# SpT1 - read in and reorganize
df.hold <- read.table(paste0(etacDir,"Avg_Betas_SpT1_LVS.txt"),header=F, sep = "\t")
df.reorganize <- matrix(NA,nrow=dim(df.hold)[1],ncol=14)
df.reorganize[,1:2] <- as.matrix(df.hold)   # LVS
df.reorganize[,3:4] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_SpT1_RVS.txt"),header=F,sep = "\t"))     #RVS
df.reorganize[,5:6] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_SpT1_RIPS.txt"),header=F,sep = "\t"))   #RIPS
df.reorganize[,7:8] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_SpT1_LAG.txt"),header=F,sep = "\t"))   #LAG
df.reorganize[,9:10] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_SpT1_RAG.txt"),header=F,sep = "\t"))   #RAG
df.reorganize[,11:12] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_SpT1_LPPCU.txt"),header=F,sep = "\t"))  #LPCU
df.reorganize[,13:14] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_SpT1_LDMPFC.txt"),header=F,sep = "\t"))  #LMPFC

# make plotable matrix
df.plotable <- matrix(NA,nrow=2,ncol=dim(df.reorganize)[2])
df.plotable[1,] <- colMeans(df.reorganize)
df.plotable[2,] <- SE.Function(df.reorganize)

# determine graph size
max_position <- which(df.plotable == max(df.plotable[1,]), arr.ind = T)
hold_max <- df.plotable[max_position] + df.plotable[2,max_position[2]]
graph_max <- hold_max+(0.1*hold_max)
min_position <- which(df.plotable == min(df.plotable[1,]), arr.ind=T)
hold_min <- df.plotable[min_position] - df.plotable[2,min_position[2]]
graph_min <- hold_min+(0.1*hold_min)

#draw graph
tiff(paste0(outDir,"Fig_Plot_ETAC_SpT1.tiff"), height = 7.5, width = 12, units = 'in', res=300)

par(mar=c(5,5,4,2),family="Times New Roman")
hold.graph <- barplot(matrix(df.plotable[1,],nrow=2), beside=T, ylim=c(graph_min,graph_max),ylab="Beta Estimate",cex.lab=1.5,legend=c("Hit","FA"),col=h.colors[1:2])
abline(h=0)
segments(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,])
arrows(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
mtext(c("A","B","C","D","E","F","G"), side=1, at=hold.graph[1,]+0.5, cex=1.2)
title(main=list("Exploratory: Study preceding Test1", cex=2))

# add sig bars for every pair of bars
hold_sequence <- 1:dim(df.plotable)[2]
hold_list <- hold_sequence[seq(1,length(hold_sequence),2)]
count<-1; for(x in hold_list){
  if(df.plotable[1,x] > 0 || df.plotable[1,(x+1)] > 0){
    larger <- which(df.plotable == max(df.plotable[1,x:(x+1)]), arr.ind = T)
    segment_y <- df.plotable[larger]+df.plotable[2,larger[2]]+0.05*graph_max
  }else{
    segment_y <- 0.05*graph_max
  }
  segments(hold.graph[1,count], segment_y, hold.graph[2,count], segment_y, lwd=2)
  draw_txt <- "***"
  text(x=hold.graph[1,count]+0.5,y=segment_y+0.02*graph_max,draw_txt)
  count<-count+1
}

dev.off()



# T1
df.hold <- read.table(paste0(etacDir,"Avg_Betas_T1_RPUT.txt"),header=F, sep = "\t")
df.reorganize <- matrix(NA,nrow=dim(df.hold)[1],ncol=18)
df.reorganize[,1:2] <- as.matrix(df.hold)   # RPUT
df.reorganize[,3:4] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_T1_RMPFC.txt"),header=F,sep = "\t"))   #RMPFC
df.reorganize[,5:6] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_T1_RAG.txt"),header=F,sep = "\t"))     #RAG
df.reorganize[,7:8] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_T1_RAMTG.txt"),header=F,sep = "\t"))   #RAMTG
df.reorganize[,9:10] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_T1_RPCU.txt"),header=F,sep = "\t"))   #RPCU
df.reorganize[,11:12] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_T1_LDMPFC.txt"),header=F,sep = "\t"))   #LDMPFC
df.reorganize[,13:14] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_T1_LINS.txt"),header=F,sep = "\t"))  #LINS
df.reorganize[,15:16] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_T1_RIPS.txt"),header=F,sep = "\t"))  #RIPS
df.reorganize[,17:18] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_T1_RIFG.txt"),header=F,sep = "\t"))  #RIFG


df.plotable <- matrix(NA,nrow=2,ncol=dim(df.reorganize)[2])
df.plotable[1,] <- colMeans(df.reorganize)
df.plotable[2,] <- SE.Function(df.reorganize)

max_position <- which(df.plotable == max(df.plotable[1,]), arr.ind = T)
hold_max <- df.plotable[max_position] + df.plotable[2,max_position[2]]
graph_max <- hold_max+(0.1*hold_max)
min_position <- which(df.plotable == min(df.plotable[1,]), arr.ind=T)
hold_min <- df.plotable[min_position] - df.plotable[2,min_position[2]]
graph_min <- hold_min+(0.1*hold_min)

tiff(paste0(outDir,"Fig_Plot_ETAC_T1.tiff"), height = 7.5, width = 12, units = 'in', res=300)

par(mar=c(5,5,4,2),family="Times New Roman")
hold.graph <- barplot(matrix(df.plotable[1,],nrow=2), beside=T, ylim=c(graph_min,graph_max),ylab="Beta Estimate",cex.lab=1.5,legend=c("Hit","FA"),col=h.colors[1:2])
abline(h=0)
segments(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,])
arrows(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
mtext(c("A","B","C","D","E","F","G","H","I"), side=1, at=hold.graph[1,]+0.5, cex=1.2)
title(main=list("Exploratory: Test1", cex=2))

hold_sequence <- 1:dim(df.plotable)[2]
hold_list <- hold_sequence[seq(1,length(hold_sequence),2)]
count<-1; for(x in hold_list){
  if(df.plotable[1,x] > 0 || df.plotable[1,(x+1)] > 0){
    larger <- which(df.plotable == max(df.plotable[1,x:(x+1)]), arr.ind = T)
    segment_y <- df.plotable[larger]+df.plotable[2,larger[2]]+0.05*graph_max
  }else{
    segment_y <- 0.05*graph_max
  }
  segments(hold.graph[1,count], segment_y, hold.graph[2,count], segment_y, lwd=2)
  draw_txt <- "***"
  text(x=hold.graph[1,count]+0.5,y=segment_y+0.02*graph_max,draw_txt)
  count<-count+1
}

dev.off()


# T1pT2
df.hold <- read.table(paste0(etacDir,"Avg_Betas_T1pT2_LINS.txt"),header=F, sep = "\t")
df.reorganize <- matrix(NA,nrow=dim(df.hold)[1],ncol=4)
df.reorganize[,1:2] <- as.matrix(df.hold)   # LINS
df.reorganize[,3:4] <- as.matrix(read.table(paste0(etacDir,"Avg_Betas_T1pT2_RSMG.txt"),header=F,sep = "\t"))   #RSMG

df.plotable <- matrix(NA,nrow=2,ncol=dim(df.reorganize)[2])
df.plotable[1,] <- colMeans(df.reorganize)
df.plotable[2,] <- SE.Function(df.reorganize)

max_position <- which(df.plotable == max(df.plotable[1,]), arr.ind = T)
hold_max <- df.plotable[max_position] + df.plotable[2,max_position[2]]
graph_max <- hold_max+(0.1*hold_max)
min_position <- which(df.plotable == min(df.plotable[1,]), arr.ind=T)
hold_min <- df.plotable[min_position] - df.plotable[2,min_position[2]]
graph_min <- hold_min+(0.1*hold_min)

tiff(paste0(outDir,"Fig_Plot_ETAC_T1pT2.tiff"), height = 7.5, width = 6, units = 'in', res=300)

par(mar=c(5,5,4,5),family="Times New Roman")
hold.graph <- barplot(matrix(df.plotable[1,],nrow=2), beside=T, ylim=c(graph_min,0.6),ylab="Beta Estimate",cex.lab=1.5,legend=c("FpH","FpF"),col=h.colors[1:2])
abline(h=0)
segments(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,])
arrows(hold.graph, df.plotable[1,]-df.plotable[2,], hold.graph, df.plotable[1,]+df.plotable[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
mtext(c("A","B"), side=1, at=hold.graph[1,]+0.5, cex=1.2)
title(main=list("Exploratory: Test1 FAs preceding Test2", cex=2))

hold_sequence <- 1:dim(df.plotable)[2]
hold_list <- hold_sequence[seq(1,length(hold_sequence),2)]
count<-1; for(x in hold_list){
  if(df.plotable[1,x] > 0 || df.plotable[1,(x+1)] > 0){
    larger <- which(df.plotable == max(df.plotable[1,x:(x+1)]), arr.ind = T)
    segment_y <- df.plotable[larger]+df.plotable[2,larger[2]]+0.05*graph_max
  }else{
    segment_y <- 0.05*graph_max
  }
  segments(hold.graph[1,count], segment_y, hold.graph[2,count], segment_y, lwd=2)
  draw_txt <- "***"
  text(x=hold.graph[1,count]+0.5,y=segment_y+0.02*graph_max,draw_txt)
  count<-count+1
}

dev.off()
