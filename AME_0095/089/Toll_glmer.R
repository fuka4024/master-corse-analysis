##NO.0:install.packages
install.packages("lme4")
install.packages("ComplexHeatmap")

##NO.1:call.packages
library(lme4)
library(edgeR)


Toll_list <- data.frame(symbol = c("BmToll-1","BmToll-2","BmToll-3","BmToll-4","BmToll-5",
                                   "BmToll-6","BmToll-7","BmToll-8","BmToll-9","BmToll-10","Bmtoll-11"),
                        ID = c("NM_001123349.1","XM_004921670.4","XM_038015525.1","XM_004927259.4","XM_021348482.2",
                               "XM_012697885.3","XM_004921675.3","XM_004921685.4","XM_012691451.3","XM_004921681.4","XM_004921682.4"))
Toll_family <- c("BmToll-1","BmToll-2",
                 "BmToll-6","BmToll-7","BmToll-8","BmToll-9","BmToll-10","Bmtoll-11")
FB_EP <-read.table("FB_EP.txt",head=TRUE)
rownames(FB_EP) <- FB_EP[,1]
FB_EP <- FB_EP[,-1]
colnames(FB_EP) <- c("EP_4_3_1","EP_4_3_2","EP_4_3_3",
                     "EP_5_0_1","EP_5_0_2","EP_5_0_3",
                     "EP_5_3_1","EP_5_3_2","EP_5_3_3",
                     "FB_4_3_1","FB_4_3_2","FB_4_3_3",
                     "FB_5_0_1","FB_5_0_2","FB_5_0_3",
                     "FB_5_3_1","FB_5_3_2","FB_5_3_3")



FB_EP_2 <- FB_EP[rownames(FB_EP) %in% Toll_list$ID, ]

FB_EP_2$num <- rep(NA,nrow(FB_EP_2))
for(RI in 1:nrow(Toll_list)){
  Hit <- which(Toll_list$ID[RI]==rownames(FB_EP_2))
  FB_EP_2$num[Hit] <- RI
  
}

FB_EP_3 <- FB_EP_2[order(FB_EP_2$num),]

dat <- data.frame(matrix(rep(NA,198*5),nrow = 2*3*11*3,ncol=5))
  colnames(dat) <- c("tissue", "rep", "Toll", "stage","count")


star_4 <- c()
star_5_0 <- c()
star_5_3 <- c()
for(RI in 1:6){
  star_4 <- c(star_4,FB_EP_3[, grepl("4_3", colnames(FB_EP_3))][,RI])
  star_5_0 <- c(star_5_0,FB_EP_3[, grepl("5_0", colnames(FB_EP_3))][,RI])
  star_5_3 <- c(star_5_3,FB_EP_3[, grepl("5_3", colnames(FB_EP_3))][,RI])
}
  dat$rep   <- as.factor(rep(rep(c(1,2,3),each=11),6))
  dat$tissue <- relevel(as.factor(rep(rep(c("EP", "FB"), each = 33),3)),ref="EP")
  dat$Toll <- as.factor(rep(Toll_list$symbol,18))
  dat$stage <- relevel(as.factor(c(rep("s4d3",length(star_4)),rep("s5d0",length(star_5_0)),rep("s5d3",length(star_5_3)))),ref="s5d0")
  dat$count <- round(c(star_4,star_5_0,star_5_3))


  design_matrix <- model.matrix(~ stage + tissue + Toll + stage:tissue + stage:Toll + tissue:Toll +  stage:tissue:Toll, data=dat)
  m1 <- glmer.nb(count ~ -1+ design_matrix+ (1|rep), data = dat)
  AIC(m1)
  f <- summary(m1 <- glmer.nb(count ~ -1+ design_matrix+ (1|rep), data = dat))
write.csv(f$coefficients,"C:/Users/monom/Documents/Toll/stage_Toll_pzalue.csv")
