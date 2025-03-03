Toll_list <- data.frame(symbol = c("BmToll-1","BmToll-2","BmToll-3","BmToll-4","BmToll-5",
                                  "BmToll-6","BmToll-7","BmToll-8","BmToll-9","BmToll-10","Bmtoll-11"),
                       ID = c("NM_001123349.1","XM_004921670.4","XM_038015525.1","XM_004927259.4","XM_021348482.2",
                              "XM_012697885.3","XM_004921675.3","XM_004921685.4","XM_012691451.3","XM_004921681.4","XM_004921682.4"))

FB_EP <-read.table("FB_EP.txt",head=TRUE)
rownames(FB_EP) <- FB_EP[,1]
FB_EP <- FB_EP[,-1]
colnames(FB_EP) <- c("EP_4_3_1","EP_4_3_2","EP_4_3_3",
                     "EP_5_0_1","EP_5_0_2","EP_5_0_3",
                     "EP_5_3_1","EP_5_3_2","EP_5_3_3",
                     "FB_4_3_1","FB_4_3_2","FB_4_3_3",
                     "FB_5_0_1","FB_5_0_2","FB_5_0_3",
                     "FB_5_3_1","FB_5_3_2","FB_5_3_3")

library(edgeR)

FB_EP_2 <- FB_EP[rownames(FB_EP) %in% c("NM_001123349.1","XM_004921670.4","XM_038015525.1",
                                      "XM_012697885.3","XM_004921675.3","XM_004921685.4","XM_012691451.3","XM_004921681.4"), ]

FB_EP_2$num <- rep(NA,nrow(FB_EP_2))
for(RI in 1:nrow(Toll_list)){
  Hit <- which(Toll_list$ID[RI]==rownames(FB_EP_2))
  FB_EP_2$num[Hit] <- RI
  
}

FB_EP_3 <- FB_EP_2[order(FB_EP_2$num),]

Toll_family <- c("BmToll-1","BmToll-2","BmToll-3",
                 "BmToll-6","BmToll-7","BmToll-8","BmToll-9","BmToll-10")

FB_EP_4 <- data.frame(matrix(rep(NA,48*3),nrow = 3,ncol=48))
colnames(FB_EP_4) <- c(unlist(lapply(Toll_family,function(x) paste0("EP_",x,"_",rep(c(1,2,3))))),
                       unlist(lapply(Toll_family,function(x) paste0("FB_",x,"_",rep(c(1,2,3))))))
  

rownames(FB_EP_4) <- c("4_3", "5_0", "5_3")

a <- 1
b <- 25
for(RI in 1:length(Toll_family)){
  Toll <- FB_EP_3[RI,1:18]
  FB_EP_4[1,as.numeric(a):as.numeric(a+2)] <- Toll[1:3]
  FB_EP_4[2,as.numeric(a):as.numeric(a+2)] <- Toll[4:6]
  FB_EP_4[3,as.numeric(a):as.numeric(a+2)] <- Toll[7:9]
  FB_EP_4[1,as.numeric(b):as.numeric(b+2)] <- Toll[10:12]
  FB_EP_4[2,as.numeric(b):as.numeric(b+2)] <- Toll[13:15]
  FB_EP_4[3,as.numeric(b):as.numeric(b+2)] <- Toll[16:18]
  a <- a + 3
  b <- b + 3
}


tissue <- factor(rep(c("EP", "FB"), each = 24))  
Toll_f <- factor(rep(rep(Toll_family,each = 3),2))
Rep <- factor(rep(c(1, 2, 3), times = 16))  

tissue <- relevel(tissue,ref="FB")


y <- DGEList(FB_EP_4, group=tissue)
y <-calcNormFactors(y, method="TMM")

design <- model.matrix(~tissue+Toll_f+Rep)
rownames(design) <- colnames(y)

y <- estimateDisp(y, design)
fit <- glmFit(y, design)

