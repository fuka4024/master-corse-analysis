BB_blast_result <- read.table("BB_blast_result.txt",as.is=TRUE)
colnames(BB_blast_result) <- c("query","subject","identity","aligne-len","mismatch",
                                "gap_open","q_start","q_end","s_start","s_end","evalue","bit-score")

List_new <- unique(BB_blast_result$query)
new_old <- data.frame(matrix(NA,ncol=2,nrow=length(List_new)))
colnames(new_old) <- c("new","old")
for(VI in 1:length(List_new)){
  new_old[VI,] <- BB_blast_result[match(List_new[VI],BB_blast_result$query),1:2]
}
write.csv(x=new_old,file = "C:/Users/monom/Documents/new_old.csv")
