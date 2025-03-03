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

contrasts <- makeContrasts(                         
  EP_5_0_vs_5_3 = Time5_3 - Time5_0,                  # EPで5_0と5_3を比較 
  FB_5_0_vs_5_3 = treatFB.Time5_3 - treatFB.Time5_0,  # FBで5_0と5_3を比較
  levels = design
)


for (RI in colnames(contrasts)) {
  lrt <- glmLRT(fit, contrast = contrasts[, RI])
  tab <- as.data.frame(topTags(lrt, n=nrow(lrt$table)))
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
  write.csv(x=dat_2,file = paste0("C:/Users/monom/Documents/DE_7/",RI,"_count.csv"))
  write.csv(x=tab,file = paste0("C:/Users/monom/Documents/DE_7/",RI,"_tab.csv"))
  write.csv(x=Y,file = paste0("C:/Users/monom/Documents/DE_7/",RI,"_tab_2.csv"))

}

CoEf <- c(5,7)
Tp <- c("EP_4_3_vs_5_0","FB_4_3_vs_5_0")
for(RI in 1:length(CoEf)){
  lrt.treat <- glmLRT(fit, coef= CoEf[RI])
  tab <- as.data.frame(topTags(lrt.treat, n=nrow(lrt.treat$table)))
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
  write.csv(x=dat_2,file = paste0("C:/Users/monom/Documents/DE_7/",Tp[RI],"_count.csv"))
  write.csv(x=tab,file = paste0("C:/Users/monom/Documents/DE_7/",Tp[RI],"_tab.csv"))
  write.csv(x=Y,file = paste0("C:/Users/monom/Documents/DE_7/",Tp[RI],"_tab_2.csv"))
}


##補正後
pvalue <- c()
for(RI in colnames(contrasts)){
X <- read.csv(paste0("C:/Users/monom/Documents/DE_7/",RI,"_tab.csv"),header = TRUE)
pvalue <- c(pvalue,X$PValue)
}
pvalue_2 <- p.adjust(pvalue,method="holm")

a <- 1
for(RI in colnames(contrasts)){
X <- read.csv(paste0("C:/Users/monom/Documents/DE_7/",RI,"_tab.csv"),header = TRUE)
X$PValue <- pvalue_2[as.numeric(a):as.numeric(a+26599)]
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
write.csv(x=dat_2,file = paste0("C:/Users/monom/Documents/DE_8/",RI,"_count.csv"))
write.csv(x=X,file = paste0("C:/Users/monom/Documents/DE_8/",RI,"_tab.csv"))
write.csv(x=Y,file = paste0("C:/Users/monom/Documents/DE_8/",RI,"_tab_2.csv"))
a <- a+26600
}


X <- read.csv(paste0("C:/Users/monom/Documents/DE_8/",Tp[1],"_tab_2.csv"),header = TRUE)
Y <- read.csv(paste0("C:/Users/monom/Documents/DE_8/",colnames(contrasts)[1],"_tab_2.csv"),header = TRUE)
X$logFC_Y <- rep(NA,nrow(X))
X$UP_DOWN <- rep(NA,nrow(X))
for(RI in 1:nrow(X)){
  Hit <- which(X$X[RI] == Y$X)
  if(length(Hit)==0){
    next
  }else{
    X$logFC_Y[RI] <- Y$logFC[Hit]
    if(X$logFC[RI] > 0 & X$logFC_Y[RI] > 0){
    X$UP_DOWN[RI] <- "UP"
    }else{
    if(X$logFC[RI] < 0 & X$logFC_Y[RI] < 0){
      X$UP_DOWN[RI] <- "DOWN"
    }else{
      next
    }
    }}
}
XX <- subset(X,subset=!is.na(UP_DOWN))
write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_8/EP_UPDOWN_tab_3.csv"))

X <- read.csv(paste0("C:/Users/monom/Documents/DE_8/",Tp[2],"_tab_2.csv"),header = TRUE)
Y <- read.csv(paste0("C:/Users/monom/Documents/DE_8/",colnames(contrasts)[2],"_tab_2.csv"),header = TRUE)
X$logFC_Y <- rep(NA,nrow(X))
X$UP_DOWN <- rep(NA,nrow(X))
for(RI in 1:nrow(X)){
  Hit <- which(X$X[RI] == Y$X)
  if(length(Hit)==0){
    next
  }else{
    X$logFC_Y[RI] <- Y$logFC[Hit]
    if(X$logFC[RI] > 0 & X$logFC_Y[RI] > 0){
      X$UP_DOWN[RI] <- "UP"
    }else{
      if(X$logFC[RI] < 0 & X$logFC_Y[RI] < 0){
        X$UP_DOWN[RI] <- "DOWN"
      }else{
        next
      }
    }}
}
XX <- subset(X,subset=!is.na(UP_DOWN))
write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_8/FB_UPDOWN_tab_3.csv"))

X <- read.csv(paste0("C:/Users/monom/Documents/DE_8/",Tp[1],"_tab_2.csv"),header = TRUE)
Y <- read.csv(paste0("C:/Users/monom/Documents/DE_8/",colnames(contrasts)[1],"_tab_2.csv"),header = TRUE)
X$logFC_Y <- rep(NA,nrow(X))
X$UP_DOWN <- rep(NA,nrow(X))
for(RI in 1:nrow(X)){
  Hit <- which(X$X[RI] == Y$X)
  if(length(Hit)==0){
    next
  }else{
    X$logFC_Y[RI] <- Y$logFC[Hit]
    if(X$logFC[RI] > 0 & X$logFC_Y[RI] > 0){
      X$UP_DOWN[RI] <- "UP"
    }else{
      if(X$logFC[RI] < 0 & X$logFC_Y[RI] < 0){
        X$UP_DOWN[RI] <- "DOWN"
      }else{
        next
      }
    }}
}
XX <- subset(X,subset=!is.na(UP_DOWN))
write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_8/EP_UPDOWN_tab_3.csv"))

X <- read.csv(paste0("C:/Users/monom/Documents/DE_8/",Tp[2],"_tab_2.csv"),header = TRUE)
Y <- read.csv(paste0("C:/Users/monom/Documents/DE_8/",colnames(contrasts)[2],"_tab_2.csv"),header = TRUE)
XX$FBlogFC <- rep(NA,nrow(XX))
for(RI in 1:nrow(XX)){
  Hit <- which(XX$X[RI] == X$X)
  if(length(Hit)==0){
    next
  }else{
    XX$FBlogFC[RI] <- X$logFC[Hit]
    }
}

XX$FBlogFC_2 <- rep(NA,nrow(XX))
for(RI in 1:nrow(XX)){
  Hit <- which(XX$X[RI] == Y$X)
  if(length(Hit)==0){
    next
  }else{
    XX$FBlogFC_2[RI] <- Y$logFC[Hit]
  }
}
XX$symbol <- c("serpin-23","nrf-6","FPI-F","obst-E","Lp-c12","LRR protein")
write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_8/EP_UPDOWN_FBNN_tab_3.csv"))


dat <- as.data.frame(cpm(FB_EP_2,normalized.lib.sizes = TRUE,log = TRUE))
dat$NAlist <- rep(NA,nrow(dat))
for(II in 1:nrow(XX)){
  Hit <- which(XX$X[II] == rownames(dat))
  if(length(Hit)==0){
    next
  }else{
    dat$NAlist[Hit] <- XX$symbol[II]
  }}
dat_2 <- subset(dat,subset=!is.na(NAlist))


EP_data <- dat_2[, grepl("EP", colnames(FB_EP))]
# タイムポイントごとに抽出
EP_4_3 <- EP_data[, grepl("4_3", colnames(EP_data))]
EP_5_0 <- EP_data[, grepl("5_0", colnames(EP_data))]
EP_5_3 <- EP_data[, grepl("5_3", colnames(EP_data))]

FB_data <- dat_2[, grepl("FB", colnames(FB_EP))]
# タイムポイントごとに抽出
FB_4_3 <- FB_data[, grepl("4_3", colnames(FB_data))]
FB_5_0 <- FB_data[, grepl("5_0", colnames(FB_data))]
FB_5_3 <- FB_data[, grepl("5_3", colnames(FB_data))]

tidy_data <- data.frame(matrix(rep(NA,18*6*5),nrow = 18*6,ncol=5))
colnames(tidy_data) <- c("count","time","symbol","tissue","Rep")

a <- 1
for(RI in 1:6){
  tidy_data$count[as.numeric(a):as.numeric(a+2)] <- EP_4_3[RI,]
  a <- a + 3
}

for(RI in 1:6){
  tidy_data$count[as.numeric(a):as.numeric(a+2)] <- EP_5_0[RI,]
  a <- a + 3
}

for(RI in 1:6){
  tidy_data$count[as.numeric(a):as.numeric(a+2)] <- EP_5_3[RI,]
  a <- a + 3
}

for(RI in 1:6){
  tidy_data$count[as.numeric(a):as.numeric(a+2)] <- FB_4_3[RI,]
  a <- a + 3
}

for(RI in 1:6){
  tidy_data$count[as.numeric(a):as.numeric(a+2)] <- FB_5_0[RI,]
  a <- a + 3
}

for(RI in 1:6){
  tidy_data$count[as.numeric(a):as.numeric(a+2)] <- FB_5_3[RI,]
  a <- a + 3
}

tidy_data$count <- unlist(tidy_data$count)
tidy_data$Rep <- rep(c(1,2,3),6*6)
tidy_data$time <- rep(rep(c("4_3","5_0","5_3"),each=18),2)
tidy_data$symbol <- rep(rep(dat_2$NAlist,each=3),6)
tidy_data$tissue <- rep(c("EP","FB"),each=3*3*6)

summary_data <- tidy_data %>%
  group_by(symbol, tissue, time) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    ymin = min(count, na.rm = TRUE),
    ymax = max(count, na.rm = TRUE)
  ) %>%
  ungroup()


library(ggplot2)
png(paste0("C:/Users/monom/Documents/DE_8/EP_FB_Time.png")
    ,width=4000, height=4000, res=700 )
par(mfrow=c(2,3))
q <- ggplot()+
  geom_ribbon(data = summary_data, aes(x = time, ymin = ymin, ymax = ymax, fill = tissue, group = interaction(symbol, tissue)), alpha = 0.2) +
  geom_line(data = summary_data,aes(x=time,y=mean_count,color=tissue,group = interaction(symbol, tissue)),size=2) +
  labs(x = "time",y = "log(cpm)",colour="tissue",fill = "tissue") +
  scale_color_manual(values = c("tomato4","turquoise4"))+
  scale_fill_manual(values = c("tomato4", "turquoise4")) +
  scale_y_continuous(breaks = seq(-5,16,1),limits = c(-5,16))+
  theme(text = element_text(size = 60))+
  theme(
    # 軸ラベルのサイズ
    axis.title = element_text(size = 60),
    # 軸目盛りのサイズ
    axis.text = element_text(size = 60),
    # ファセットタイトルのサイズ
    strip.text = element_text(size = 60, face = "bold"),
    # 凡例のサイズ
    legend.title = element_text(size = 60),
    legend.text = element_text(size = 60)
  )+
  geom_smooth()+
  theme_bw()+
  facet_wrap(~factor(symbol,levels = XX$symbol))
print(q)
dev.off()


EV <- c("EP_4_3_vs_5_0","EP_5_0_vs_5_3","FB_4_3_vs_5_0","FB_5_0_vs_5_3")
plots <- list()
for(RI in 1:length(EV)){
  X <- read.csv(paste0("C:/Users/monom/Documents/DE_8/",EV[RI],"_tab.csv"),header = TRUE)
  X <- mutate(X, thcolor = ifelse(FDR >= 0.05,0,
                                  ifelse(logFC >= 1,1,
                                         ifelse(logFC <= -1,2,
                                                0))))
  for(II in 1:nrow(dat_2)){
    Hit <- which(rownames(dat_2)[II] == X$X)
    X$thcolor[Hit] <- dat_2$NAlist[II]
  }
  
  write.csv(x=X,paste0("C:/Users/monom/Documents/DE_8/",EV[RI],"_FC_signal.csv"))

  # Yのデータフレームに、条件に基づいたフィルターを適用
  Y <- X %>%
    filter(!thcolor %in% c("0", "1", "2") & FDR < 0.05  & PValue < 0.05) %>%
    mutate(FDR_2 = FDR)  # -log10(FDR)を新しい列として追加
  
  
  p <- ggplot()+
    geom_point(data=X,mapping =aes(x = logFC, y = -log10(FDR)),color = "gray88",alpha=0.5 )+
    geom_point(data = X %>% filter(thcolor == "2"),mapping = aes(x = logFC, y = -log10(FDR)),color = "lightsteelblue",alpha=0.3,size = 2.0)+
    geom_point(data = X %>% filter(thcolor == "1"),mapping = aes(x = logFC, y = -log10(FDR)),color = "lightsalmon",alpha=0.3,size = 2.0)+
    geom_point(data = Y,mapping = aes(x = logFC, y = -log10(FDR)),color = "red",size = 1.0)+
    geom_text_repel(data = Y,aes(x = logFC, y = -log10(FDR),label = thcolor),size = 4,max.overlaps = 20, box.padding = 0.5,vjust = -1, size = 3, color = "black")+
    geom_hline(yintercept = -log10(0.05),size = 0.2,color = "dark green") +
    geom_vline(xintercept = -1 ,size = 0.2,color = "dark green") +
    geom_vline(xintercept = 1 ,size = 0.2,color = "dark green") +
    theme(legend.position = "none", panel.grid = element_blank(),
          # 軸ラベルのサイズ
          axis.title = element_text(size = 60),
          # 軸目盛りのサイズ
          axis.text = element_text(size = 60),
          # ファセットタイトルのサイズ
          strip.text = element_text(size = 60, face = "bold"),
          # 凡例のサイズ
          legend.title = element_text(size = 60),
          legend.text = element_text(size = 60))+
    labs(x = "log2 (Fold Change)",y = "-log10 (FDR)",title = EV[RI]) +
    theme_bw()
  
  plots[[RI]] <- p
  
}

png("C:/Users/monom/Documents/DE_8/EP_FB_DE_Time.png", width = 4000, height = 4000, res = 700)
plot_grid(plotlist = plots, ncol = 2)
dev.off()

