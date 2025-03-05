blast_result <- read.table("BmToll_Dm_blast_result.txt")
colnames(blast_result) <- c("query","subject","identity","aligne-len","mismatch",
                                  "gap_open","q_start","q_end","s_start","s_end","evalue","bit-score")

sort <- blast_result[sort(as.numeric(blast_result$evalue),decreasing = F,index= T)$ix,]
List <- unique(blast_result$query)
top <- sort[match(as.character(List),as.character(sort$query)),]
rownames(top) <- top$query
best_rows <- unlist(lapply(as.character(List),
                                        function(x) intersect(which(blast_result$query==x),
                                                              which(blast_result$evalue==top[x,"evalue"])
                                        )))
best <- blast_result[best_rows,]

library(data.table)
row_table <- fread("fbgn_NAseq_Uniprot_fb_2023_04.tsv",showProgress = FALSE,encoding = "UTF-8",data.table = FALSE)
row_table[,1:2] <- list(NULL)
row_table[,2:5] <- list(NULL)
is_blank <- function(x){is.na(x) | x == ""}
delete_row <- apply(as.matrix(row_table[,2]),1,
                    function(x){
                      all(is_blank(x))
                    })
Drosophila_idlist <- row_table[!delete_row,]
write.csv(x=Drosophila_idlist,file = "/Users/monom/Documents/Drosophila_idlist.csv")

fbgn_symbal_list <- readLines("fbgn_annotation_ID_fb_2024_05.tsv")
fbgn_symbal_list <- fbgn_symbal_list[6:17877]
fbgn_symbal_2 <- strsplit(fbgn_symbal_list,"\tDmel\t")
fbgn_symbol <- c()
for(RI in 1:length(fbgn_symbal_2)){
  fbgn_symbol <- c(fbgn_symbol,fbgn_symbal_2[[RI]][1])
}
fbgn_line <- c()
for(RI in 1:length(fbgn_symbal_2)){
  fbgn_line <- c(fbgn_line,fbgn_symbal_2[[RI]][2])
}
fbgn_line_2 <- strsplit(fbgn_line,"\t")
fbgn <- c()
for(RI in 1:length(fbgn_line_2)){
  fbgn <- c(fbgn,fbgn_line_2[[RI]][1])
}
fbgn_matrix <- as.data.frame(matrix(NA,ncol=3,nrow=length(fbgn)))
colnames(fbgn_matrix) <- c("Fb_Symbol","Fbgn","Td")
fbgn_matrix$Fb_Symbol <- fbgn_symbol
fbgn_matrix$Fbgn <- fbgn
Fbgn_ref <- read.csv("Drosophila_idlist.csv")
colnames(Fbgn_ref) <- c("Fb_Symbol","Fbgn","Td","Pd")
for(RI in 1:nrow(Fbgn_ref)){
  Hit <- which(Fbgn_ref$Fbgn[RI] ==fbgn_matrix$Fbgn)
  if(length(Hit) == 0){
    next
  }else{
    Fbgn_ref$Fb_Symbol[RI] <- fbgn_matrix$Fb_Symbol[Hit]
  }
}
write.csv(x=Fbgn_ref,file = "/Users/monom/Documents/Drosophila_symbol.csv")

best$Dm_Symbol <- rep(NA,nrow(best))

for(RI in 1:nrow(best)){
  Hit <- which(str_replace_all(best$subject[RI],paste0("\\.","[0-9]"),"") == Fbgn_ref$Pd)
  if(length(Hit) == 0){
    next
  }else{
    best$Dm_Symbol[RI] <- Fbgn_ref$Fb_Symbol[Hit]
  }
}
write.csv(x=best,file = "/Users/monom/Documents/BmToll_Dm_best_symbol.csv")

