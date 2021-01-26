### --- Notes
#
# Post hoc tests to see if subjects are more
# likely to respond one way in Test2 if they 
# saw the relevant image more recently in 
# Test1.


# Set paths
parent_dir = "/Users/nmuncy/Projects/STT/data/behAnalysis/"
data_dir = paste0(parent_dir, "data_indiv/")
out_dir = paste0(parent_dir, "data_group/")

# Get list of subjects
ref_table = read.delim(paste0(out_dir, "Group_table_dprime_T1.txt"), header = T)
subj_list = ref_table$Subj

# start group df
df_group <- as.data.frame(matrix(NA, nrow = dim(ref_table)[1], ncol = 5))
colnames(df_group) <- c("Subj", "Num.CrFa", "Dist.CrFa", "Num.CrHit", "Dist.CrHit")

out_row <- 1
for(subj in subj_list){

  # get data
  df_t1 <- read.delim(paste0(data_dir, subj, "_DF_T1.txt"), header=T)
  df_t2 <- read.delim(paste0(data_dir, subj, "_DF_T2.txt"), header=T)
  
  # make combined data frame (for distance calc)
  df_comb <- as.data.frame(matrix(NA, nrow = (dim(df_t1)[1]+dim(df_t2)[1]), ncol = 4))
  colnames(df_comb) <- c("Test", "Stim", "Beh", "Type")
  df_comb$Test <- c(rep("1", dim(df_t1)[1]), rep("2", dim(df_t2)[1]))
  df_comb$Stim <- c(df_t1$Image, df_t2$Image)
  df_comb$Beh <- c(df_t1$Behavior, df_t2$Behavior)
  df_comb$Type <- c(df_t1$StimType, rep("Both", dim(df_t2)[1]))
  ind_t1 <- grep(1, df_comb$Test)
  ind_t2 <- grep(2, df_comb$Test)
  
  # make recency data frame
  df_rec <- as.data.frame(matrix(NA, nrow = dim(df_t2)[1], ncol = 4))
  colnames(df_rec) <- c("Stim", "T1type", "Dist", "Beh")
  
  # find t1 correspondence of t2 stim
  #   determine distance, type, and behavior
  for(i in 1:length(ind_t2)){
    h_ind <- ind_t2[i]
    t2_stim <- df_comb$Stim[h_ind]
    df_rec$Stim[i] <- t2_stim
    t1_loc <- which(df_comb$Stim[ind_t1] == t2_stim)
    df_rec$T1type[i] <- df_comb$Type[t1_loc]
    df_rec$Dist[i] <- as.numeric(h_ind) - as.numeric(t1_loc)
    df_rec$Beh[i] <- paste0(df_comb$Beh[t1_loc], "-", df_comb$Beh[h_ind])
  }
  
  # find position, number of relevant behaviors
  pos_CrFa <- which(df_rec$Beh == "CR-FA")
  num_CrFa <- as.numeric(length(pos_CrFa))
  
  pos_CrHit <- which(df_rec$Beh == "CR-Hit")
  num_CrHit <- as.numeric(length(pos_CrHit))

  # find average stimulus distance for e/beh
  avg_CrFa <- round(mean(df_rec$Dist[pos_CrFa]))
  avg_CrHit <- round(mean(df_rec$Dist[pos_CrHit]))
  
  # append out
  df_group[out_row,] <- c(subj, num_CrFa, avg_CrFa, num_CrHit, avg_CrHit)
  out_row <- out_row+1
}

# do proportion testing
Prop.Function <- function(a,b){
  func.df <- matrix(0,nrow=length(a),ncol=2)
  func.df[,1] <- a
  func.df[,2] <- b
  input <- colMeans(func.df)
  output <- prop.test(input[1],input[2],p=0.5,correct=TRUE)
  return(output)
}
prop_out <- Prop.Function(as.numeric(df_group$Num.CrFa), (as.numeric(df_group$Num.CrFa) + as.numeric(df_group$Num.CrHit)))

# write out
out_table <- paste0(out_dir, "Group_table_recency.txt")
write.table(df_group, out_table, quote = F, sep = "\t", row.names = F)

out_stats <- paste0(out_dir, "Stats_prop_recency.txt")
h_out <- capture.output(prop_out)
write.table(h_out, out_stats, row.names = F, quote = F)

