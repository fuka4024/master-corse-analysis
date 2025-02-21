library(edgeR)
library(stringr)

FB_EP <-read.table("FB_EP.txt",head=TRUE)
rownames(FB_EP) <- FB_EP[,1]
FB_EP <- FB_EP[,-1]
colnames(FB_EP) <- c("EP_4_3_1","EP_4_3_2","EP_4_3_3",
                     "EP_5_0_1","EP_5_0_2","EP_5_0_3",
                     "EP_5_3_1","EP_5_3_2","EP_5_3_3",
                     "FB_4_3_1","FB_4_3_2","FB_4_3_3",
                     "FB_5_0_1","FB_5_0_2","FB_5_0_3",
                     "FB_5_3_1","FB_5_3_2","FB_5_3_3")

FB_EP_2 <- FB_EP[-(which(rowSums(FB_EP) == 0)),]


treat <- factor(rep(c("EP","FB"),each = 9))
treat <- relevel(treat,ref="EP")

Time <- factor(c(rep(c("4_3","5_0","5_3"),each=3),rep(c("4_3","5_0","5_3"),each=3)))
Rep <- factor(rep(c(1,2,3),6))

y <- DGEList(FB_EP_2, group=treat)
y <-calcNormFactors(y, method="TMM")

design <- model.matrix(~Rep+treat+Time+treat:Time)
rownames(design) <- colnames(y)

y <- estimateDisp(y, design)
fit <- glmFit(y, design)

dat <- as.data.frame(cpm(FB_EP_2,normalized.lib.sizes = TRUE))

Tp <- c("4_3","5_0","5_3")

colnames(design) <- make.names(colnames(design))

contrast_4_3_vs_5_0 <- makeContrasts(
  contrast = Time5_0,
  levels = design
)

contrast_5_0_vs_5_3 <- makeContrasts(
  contrast = Time5_3 - Time5_0,
  levels = design
)

# 4_3 vs 5_0
lrt_4_3_vs_5_0 <- glmLRT(fit, contrast = contrast_4_3_vs_5_0)
tab <- as.data.frame(topTags(lrt_4_3_vs_5_0, n=nrow(lrt_4_3_vs_5_0$table)))
Y <- tab[which(tab$FDR < 0.05),]
dat$NAlist <- rep(NA,nrow(dat))
for(II in 1:nrow(Y)){
  Hit <- which(rownames(Y)[II] == rownames(dat))
  if(length(Hit)==0){
    next
  }else{
    dat$NAlist[Hit] <- "Hit"
  }}
dat_2 <- subset(dat,subset=!is.na(NAlist))
write.csv(x=dat_2,file = paste0("C:/Users/monom/Documents/DE_5/4_3_vs_5_0_count.csv"))
write.csv(x=tab,file = paste0("C:/Users/monom/Documents/DE_5/4_3_vs_5_0_tab.csv"))
write.csv(x=Y,file = paste0("C:/Users/monom/Documents/DE_5/4_3_vs_5_0_tab_2.csv"))

# 5_0 vs 5_3
lrt_5_0_vs_5_3 <- glmLRT(fit, contrast = contrast_5_0_vs_5_3)

tab <- as.data.frame(topTags(lrt.treat, n=nrow(lrt_5_0_vs_5_3$table)))
Y <- tab[which(tab$FDR < 0.05),]
dat$NAlist <- rep(NA,nrow(dat))
for(II in 1:nrow(Y)){
  Hit <- which(rownames(Y)[II] == rownames(dat))
  if(length(Hit)==0){
    next
  }else{
    dat$NAlist[Hit] <- "Hit"
  }}
dat_2 <- subset(dat,subset=!is.na(NAlist))
write.csv(x=dat_2,file = paste0("C:/Users/monom/Documents/DE_5/5_0_vs_5_3_count.csv"))
write.csv(x=tab,file = paste0("C:/Users/monom/Documents/DE_5/5_0_vs_5_3_tab.csv"))
write.csv(x=Y,file = paste0("C:/Users/monom/Documents/DE_5/5_0_vs_5_3_tab_2.csv"))

pvalue <- c()
  X <- read.csv(paste0("C:/Users/monom/Documents/DE_5/4_3_vs_5_0_tab.csv"),header = TRUE)
  y <- read.csv(paste0("C:/Users/monom/Documents/DE_5/5_0_vs_5_3_tab.csv"),header = TRUE)
  pvalue <- c(pvalue,X$PValue,y$PValue)
  pvalue_2 <- p.adjust(pvalue,method="holm")
  
  
  X <- read.csv(paste0("C:/Users/monom/Documents/DE_5/4_3_vs_5_0_tab.csv"),header = TRUE)
  X$PValue <- pvalue_2[1:26600]
  Y <- X[which(X$FDR < 0.05 & X$PValue < 0.05),]
  dat$NAlist <- rep(NA,nrow(dat))
  for(II in 1:nrow(Y)){
    Hit <- which(rownames(Y)[II] == rownames(dat))
    if(length(Hit)==0){
      next
    }else{
      dat$NAlist[Hit] <- "Hit"
    }}
  dat_2 <- subset(dat,subset=!is.na(NAlist))
  write.csv(x=dat_2,file = paste0("C:/Users/monom/Documents/DE_6/4_3_vs_5_0_count.csv"))
  write.csv(x=X,file = paste0("C:/Users/monom/Documents/DE_6/4_3_vs_5_0_tab.csv"))
  write.csv(x=Y,file = paste0("C:/Users/monom/Documents/DE_6/4_3_vs_5_0_tab_2.csv"))
  
  X <- read.csv(paste0("C:/Users/monom/Documents/DE_5/5_0_vs_5_3_tab.csv"),header = TRUE)
  X$PValue <- pvalue_2[26601:53200]
  Y <- X[which(X$FDR < 0.05 & X$PValue < 0.05),]
  dat$NAlist <- rep(NA,nrow(dat))
  for(II in 1:nrow(Y)){
    Hit <- which(rownames(Y)[II] == rownames(dat))
    if(length(Hit)==0){
      next
    }else{
      dat$NAlist[Hit] <- "Hit"
    }}
  dat_2 <- subset(dat,subset=!is.na(NAlist))
  write.csv(x=dat_2,file = paste0("C:/Users/monom/Documents/DE_6/5_0_vs_5_3_count.csv"))
  write.csv(x=X,file = paste0("C:/Users/monom/Documents/DE_6/5_0_vs_5_3_tab.csv"))
  write.csv(x=Y,file = paste0("C:/Users/monom/Documents/DE_6/5_0_vs_5_3_tab_2.csv"))
  
X <- read.csv(paste0("C:/Users/monom/Documents/DE_6/4_3_vs_5_0_tab_2.csv"),header = TRUE)
Y <- read.csv(paste0("C:/Users/monom/Documents/DE_6/5_0_vs_5_3_tab_2.csv"),header = TRUE)
TR <- c()
TR <- c(TR,X$X,Y$X)
TR <- unique(TR)

dat <- as.data.frame(cpm(FB_EP_2,normalized.lib.sizes = TRUE,log = TRUE))
dat$NAlist <- rep(NA,nrow(dat))
for(II in 1:length(TR)){
  Hit <- which(TR[II] == rownames(dat))
  if(length(Hit)==0){
    next
  }else{
    dat$NAlist[Hit] <- "Hit"
  }}
dat_2 <- subset(dat,subset=!is.na(NAlist))

# EPに該当する列を抽出
EP_data <- dat_2[, grepl("EP", colnames(FB_EP))]
# タイムポイントごとに抽出
EP_4_3 <- EP_data[, grepl("4_3", colnames(EP_data))]
EP_5_0 <- EP_data[, grepl("5_0", colnames(EP_data))]
EP_5_3 <- EP_data[, grepl("5_3", colnames(EP_data))]
# 各タイムポイントでの遺伝子ごとの平均値
EP_means <- data.frame(
  "4_3" = rowMeans(EP_4_3),
  "5_0" = rowMeans(EP_5_0),
  "5_3" = rowMeans(EP_5_3)
)

# タイムポイントを数値に変換
time_points <- c(4.3, 5.0, 5.3)

# 傾きを格納するベクトルを初期化
slopes <- numeric(nrow(EP_means))

# 各遺伝子ごとに線形回帰を適用
for (i in 1:nrow(EP_means)) {
  # 線形回帰モデルを適用
  fit <- lm(as.numeric(EP_means[i, ]) ~ time_points)
  
  # 傾きを保存
  slopes[i] <- coef(fit)[2]  # 傾きは回帰係数の2番目
}

# 結果をデータフレーム化
results <- data.frame(Gene = rownames(EP_means), Slope = slopes)
UD <- results[-which(results$Slope < 1 & results$Slope > -1),]

UD_Gene <- UD$Gene
dat_2$slope <- rep(NA,nrow(dat_2))

for(II in 1:length(UD_Gene)){
  Hit <- which(UD_Gene[II] == rownames(dat_2))
  if(length(Hit)==0){
    next
  }else{
    dat_2$slope[Hit] <- UD$Slope[II]
  }}
dat_3 <- subset(dat_2,subset=!is.na(slope))

tidy_data_2 <- data.frame(matrix(rep(NA,nrow(dat_3)*4),nrow = nrow(dat_3),ncol=4))
colnames(tidy_data_2) <- c("symbol","Mean","Sqrt","CV")
tidy_data_2$symbol <- rownames(dat_3)
for(RI in 1:nrow(dat_3)){
  tidy_data_2$Mean[RI] <- rowMeans(dat_3[RI, grepl("FB", colnames(FB_EP))])
  tidy_data_2$Sqrt[RI] <- sd(dat_3[RI, grepl("FB", colnames(FB_EP))])
  tidy_data_2$CV[RI] <- tidy_data_2$Sqrt[RI]/tidy_data_2$Mean[RI]
}

png(paste0("C:/Users/monom/Documents/DE_6/hist.png")
    ,width=3000, height=4000, res=200 )
hist(tidy_data_2$CV[tidy_data_2$CV >= -1 & tidy_data_2$CV <= 1], 
     breaks = seq(-1, 1, 0.1), 
     main = "CV Between -1 and 1", 
     xlab = "CV", 
     col = "green")
dev.off() 

CVG <- tidy_data_2[tidy_data_2$CV >= -1 & tidy_data_2$CV <= 1,]
dat_3$CV <- rep(NA,nrow(dat_3))
for(RI in 1:nrow(CVG)){
  Hit <- which(rownames(dat_3) == CVG$symbol[RI])
  if(length(Hit) == 0){
    next
  }else{
    dat_3$CV[Hit] <- CVG$CV[RI]
  }
}
dat_4 <- subset(dat_3,subset=!is.na(CV))

Ortholog_GO_D <- read.csv("C:/Users/monom/Documents/Ortholog_GO_Drosophila.csv",header = TRUE)
BmoriID_s <- read.csv("C:/Users/monom/Documents/Symbol_BmoriID.csv",header = TRUE)


  DEGs <- rownames(dat_4)
  Ortholog_GO_D$transcript <- rep(NA,nrow(Ortholog_GO_D))
  Ortholog_GO_D$Model <- rep(NA,nrow(Ortholog_GO_D))
  Ortholog_GO_D$FDR <- rep(NA,nrow(Ortholog_GO_D))
  for(SI in 1:length(DEGs)){
    Hit <- which(DEGs[SI] == BmoriID_s$Tb)
    if(length(Hit)==0){
      next
    }else{
      SecondHit <- which(BmoriID_s$Gene[Hit] == Ortholog_GO_D$Gb_p50)
      if(length(SecondHit)==0){
        next
      }else{
        Ortholog_GO_D$transcript[SecondHit] <- DEGs[SI]
        Ortholog_GO_D$Model[SecondHit] <- "Hit" 
      }
    }
  }
  X <- subset(Ortholog_GO_D,subset=!is.na(Model))
  write.csv(x=X,file = paste0("C:/Users/monom/Documents/DE_6/GO.csv"))


GOterm <- read.csv("C:/Users/monom/Documents/Tb_Goterm_list.csv",header = TRUE)

  XX <- X
  XX$Model <- rep(NA,nrow(XX))
  for(QI in 1:nrow(XX)){
    Hit <- which(XX$GO[QI] == GOterm$Goid)
    if(length(Hit)==0){
      next
    }else{
      XX$Model[QI] <- GOterm$Goterm[Hit]
    }
    
  }
  XX <- XX[order(XX$FDR),]
  write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_6/GOterm.csv"))


  XXX <- XX
  XXX$symbol <- rep(NA,nrow(XXX))
  for(QI in 1:nrow(XXX)){
    Hit <- which(XXX$Gb_p50[QI] == BmoriID_s$Gene)
    if(length(Hit)==0){
      next
    }else{
      XXX$symbol[QI] <- BmoriID_s$Fb_Symbol[Hit[1]]
    }
    
  }
  write.csv(x=XXX,file = paste0("C:/Users/monom/Documents/DE_6/symbol.csv"))

  TR_list <- data.frame(symbol = c("mamo","lwr"),
                         ID = c("XM_021346602.2","NM_001046942.1"))
  dat_4$NAlist <- rep(NA,nrow(dat_4))
  for(RI in 1:nrow(TR_list)){
    Hit <- which(TR_list$ID[RI]==rownames(dat_4))
    if(length(Hit)==0){
      next
    }else{
      dat_4$NAlist[Hit] <- TR_list$symbol[RI]
    }
  }
  dat_5 <- subset(dat_4,subset=!is.na(NAlist))
  
  EP_data <- dat_5[, grepl("EP", colnames(FB_EP))]
  # タイムポイントごとに抽出
  EP_4_3 <- EP_data[, grepl("4_3", colnames(EP_data))]
  EP_5_0 <- EP_data[, grepl("5_0", colnames(EP_data))]
  EP_5_3 <- EP_data[, grepl("5_3", colnames(EP_data))]
  
  FB_data <- dat_5[, grepl("FB", colnames(FB_EP))]
  # タイムポイントごとに抽出
  FB_4_3 <- FB_data[, grepl("4_3", colnames(FB_data))]
  FB_5_0 <- FB_data[, grepl("5_0", colnames(FB_data))]
  FB_5_3 <- FB_data[, grepl("5_3", colnames(FB_data))]
  
  tidy_data <- data.frame(matrix(rep(NA,36*5),nrow = 36,ncol=5))
  colnames(tidy_data) <- c("count_mean","time","symbol","tissue","Rep")
  tidy_data$count_mean[1:3] <- EP_4_3[1,]
  tidy_data$count_mean[4:6] <- EP_4_3[2,]
  tidy_data$count_mean[7:9] <- EP_5_0[1,]
  tidy_data$count_mean[10:12] <- EP_5_0[2,]
  tidy_data$count_mean[13:15] <- EP_5_3[1,]
  tidy_data$count_mean[16:18] <- EP_5_3[2,]
  tidy_data$count_mean[19:21] <- FB_4_3[1,]
  tidy_data$count_mean[22:24] <- FB_4_3[2,]
  tidy_data$count_mean[25:27] <- FB_5_0[1,]
  tidy_data$count_mean[28:30] <- FB_5_0[2,]
  tidy_data$count_mean[31:33] <- FB_5_3[1,]
  tidy_data$count_mean[34:36] <- FB_5_3[2,]
  tidy_data$Rep <- rep(c(1,2,3),12)
  tidy_data$time <- rep(rep(c("4_3","5_0","5_3"),each=3),4)
  tidy_data$symbol <- rep(c("lwr","mamo"),each=3,6)
  tidy_data$tissue <- rep(c("EP","FB"),each=18)
  
  
  png(paste0("C:/Users/monom/Documents/DE_6/EP_FB_Time.png")
      ,width=6000, height=2000, res=500 )
  par(mfrow=c(1,2))
  q <- ggplot(tidy_data,aes(x=time,y=count_mean,color=tissue,group = symbol))+
    geom_line(size=2) +
    labs(x = "time",y = "log(cpm)",colour="tissue") +
    scale_color_manual(values = c("tomato4","turquoise4"))+
    scale_y_continuous(breaks = seq(1,10,1),limits = c(1,10))+
    theme(text = element_text(size = 60))+
    geom_smooth()+
    theme_bw()+
    facet_wrap(~factor(symbol,levels = c("lwr","mamo")))
  print(q)
  dev.off()
  
  library(ggplot2)
  
  png(paste0("C:/Users/monom/Documents/DE_6/EP_FB_Time.png")
      ,width=6000, height=2000, res=500 )
  par(mfrow=c(1,2))
  # プロットの作成
  p <- ggplot(tidy_data, aes(x = time, y = count_mean, color = tissue, group = tissue)) +
    geom_line(size = 1.5) + # ラインプロット
    geom_point(size = 3) +  # 各データポイント
    facet_wrap(~symbol) +   # 遺伝子ごとに分割
    labs(title = "Count Mean Over Time by Tissue", 
         x = "Time", 
         y = "Count Mean", 
         color = "Tissue") +
    theme_minimal(base_size = 15) + # シンプルなテーマ
    scale_color_manual(values = c("tomato4", "turquoise4")) # 色の指定
  
  # プロットの表示
  print(p)
  dev.off()
  
  library(ggplot2)
  png(paste0("C:/Users/monom/Documents/DE_6/EP_FB_Time.png")
      ,width=6000, height=2000, res=500 )
  par(mfrow=c(1,2))
  # プロットの作成
  p <- ggplot(tidy_data, aes(x = time, y = count_mean, color = tissue)) +
    geom_point(aes(shape = as.factor(Rep)), size = 3, alpha = 0.7) + # リプリケートの点
    stat_summary(fun = mean, geom = "line", aes(group = tissue), size = 1.5) + # 平均値のライン
    facet_wrap(~symbol) + # 遺伝子ごとに分割
    labs(title = "Count Mean Over Time by Tissue and Symbol",
         x = "Time",
         y = "Count Mean",
         color = "Tissue",
         shape = "Replicate") +
    theme_minimal(base_size = 15) +
    scale_color_manual(values = c("tomato4", "turquoise4")) # 組織ごとの色分け
  
  # プロットの表示
  print(p)
  dev.off()
  
  
  
  
  
  
  Ortholog_GO_D <- read.csv("C:/Users/monom/Documents/Ortholog_GO_Drosophila.csv",header = TRUE)
  BmoriID_s <- read.csv("C:/Users/monom/Documents/Symbol_BmoriID.csv",header = TRUE)
  

    X <- read.csv(paste0("C:/Users/monom/Documents/DE_8/FB_5_0_vs_5_3_tab_2.csv"),header = TRUE)
    DEGs <- X$X
    Ortholog_GO_D$transcript <- rep(NA,nrow(Ortholog_GO_D))
    Ortholog_GO_D$Model <- rep(NA,nrow(Ortholog_GO_D))
    Ortholog_GO_D$FDR <- rep(NA,nrow(Ortholog_GO_D))
    for(SI in 1:length(DEGs)){
      Hit <- which(DEGs[SI] == BmoriID_s$Tb)
      if(length(Hit)==0){
        next
      }else{
        SecondHit <- which(BmoriID_s$Gene[Hit] == Ortholog_GO_D$Gb_p50)
        if(length(SecondHit)==0){
          next
        }else{
          Ortholog_GO_D$FDR[SecondHit] <- X$FDR[SI]
          Ortholog_GO_D$transcript[SecondHit] <- DEGs[SI]
          Ortholog_GO_D$Model[SecondHit] <- "Hit" 
        }
      }
    }
    X <- subset(Ortholog_GO_D,subset=!is.na(Model))
    write.csv(x=X,file = paste0("C:/Users/monom/Documents/DE_8/FB_5_0_vs_5_3_GO.csv"))
  
  
  GOterm <- read.csv("C:/Users/monom/Documents/Tb_Goterm_list.csv",header = TRUE)
    XX <- read.csv(paste0("C:/Users/monom/Documents/DE_8/FB_5_0_vs_5_3_GO.csv"),header = TRUE)
    XX$Model <- rep(NA,nrow(XX))
    for(QI in 1:nrow(XX)){
      Hit <- which(XX$GO[QI] == GOterm$Goid)
      if(length(Hit)==0){
        next
      }else{
        XX$Model[QI] <- GOterm$Goterm[Hit]
      }
      
    }
    XX <- XX[order(XX$FDR),]
    write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_8/FB_5_0_vs_5_3_GOterm.csv"))
  
  

    XX <- read.csv(paste0("C:/Users/monom/Documents/DE_8/FB_5_0_vs_5_3_GOterm.csv"),header = TRUE)
    XX$symbol <- rep(NA,nrow(XX))
    for(QI in 1:nrow(XX)){
      Hit <- which(XX$Gb_p50[QI] == BmoriID_s$Gene)
      if(length(Hit)==0){
        next
      }else{
        XX$symbol[QI] <- BmoriID_s$Fb_Symbol[Hit[1]]
      }
      
    }
    write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_8/FB_5_0_vs_5_3_symbol.csv"))
  
  