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

# start output df
df_out <- as.data.frame(matrix(NA, nrow = dim(ref_table)[1], ncol = 7))
colnames(df_out) <- c("Subj", "AvgFF", "AvgFH", "AvgHF", "AvgHH", "Sig.FFH", "Sig.HFH")

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
  
  # find position of relevant behaviors
  pos_FaFa <- which(df_rec$Beh == "FA-FA")
  pos_FaHit <- which(df_rec$Beh == "FA-Hit")
  pos_HitHit <- which(df_rec$Beh == "Hit-Hit")
  pos_HitFa <- which(df_rec$Beh == "Hit-FA")
  
  # find average distance for e/beh
  avg_FaFa <- round(mean(df_rec$Dist[pos_FaFa]))
  avg_FaHit <- round(mean(df_rec$Dist[pos_FaHit]))
  avg_HitFa <- round(mean(df_rec$Dist[pos_HitFa]))
  avg_HitHit <- round(mean(df_rec$Dist[pos_HitHit]))
  
  # test for difference, account for unequal num (test for equal variance)
  lure_vout <- var.test(df_rec$Dist[pos_FaFa], df_rec$Dist[pos_FaHit])
  if(lure_vout$p.value > 0.05){
    lure_out <- t.test(df_rec$Dist[pos_FaFa], df_rec$Dist[pos_FaHit], var.equal=T)
  }else{
    lure_out <- t.test(df_rec$Dist[pos_FaFa], df_rec$Dist[pos_FaHit], var.equal=F)
  }
  
  # symbolize significance
  if(lure_out$p.value < 0.05){
    lure_sig <- "*"
  }else{
    lure_sig <- "--"
  }
  
  targ_vout <- var.test(df_rec$Dist[pos_HitFa], df_rec$Dist[pos_HitHit])
  if(targ_vout$p.value > 0.05){
    targ_out <- t.test(df_rec$Dist[pos_HitFa], df_rec$Dist[pos_HitHit], var.equal=T)
  }else{
    targ_out <- t.test(df_rec$Dist[pos_HitFa], df_rec$Dist[pos_HitHit], var.equal=F)
  }
  if(targ_out$p.value < 0.05){
    targ_sig <- "*"
  }else{
    targ_sig <- "--"
  }
  
  # append out
  df_out[out_row,] <- c(subj, avg_FaFa, avg_FaHit, avg_HitFa, avg_HitHit, lure_sig, targ_sig)
  out_row <- out_row+1
}

# write out
out_file <- paste0(out_dir, "Group_table_recency.txt")
write.table(df_out, out_file, quote = F, sep = "\t", row.names = F)
