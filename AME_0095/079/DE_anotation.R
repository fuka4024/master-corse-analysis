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
Tp <- c("4_3","5_0","5_3")
A <- 1
for(RI in 1:length(Tp)){
  W <- cbind(FB_EP_2[,as.numeric(A):as.numeric(A+2)],FB_EP_2[,as.numeric(A+9):as.numeric(A+11)])
  treat <- factor(rep(c("EP","FB"),each = 3))
  treat <- relevel(treat,ref="EP")
  Rep <- factor(rep(c(1,2,3),2))
  y <- DGEList(W, group=treat)
  y <-calcNormFactors(y, method="TMM")

design <- model.matrix(~Rep+treat)
rownames(design) <- colnames(y)

y <- estimateDisp(y, design)
fit <- glmFit(y, design)

dat <- as.data.frame(cpm(W,normalized.lib.sizes = TRUE))

  lrt.treat <- glmLRT(fit, coef= 4)
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
  write.csv(x=dat_2,file = paste0("C:/Users/monom/Documents/DE_2/",Tp[RI],"_count.csv"))
  write.csv(x=tab,file = paste0("C:/Users/monom/Documents/DE_2/",Tp[RI],"_tab.csv"))
  write.csv(x=Y,file = paste0("C:/Users/monom/Documents/DE_2/",Tp[RI],"_tab_2.csv"))
  A <- A+3
}

Ortholog_GO_D <- read.csv("C:/Users/monom/Documents/Ortholog_GO_Drosophila.csv",header = TRUE)
BmoriID_s <- read.csv("C:/Users/monom/Documents/Symbol_BmoriID.csv",header = TRUE)

for(RI in 1:length(Tp)){
  X <- read.csv(paste0("C:/Users/monom/Documents/DE_2/",Tp[RI],"_tab_2.csv"),header = TRUE)
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
  write.csv(x=X,file = paste0("C:/Users/monom/Documents/DE_2/",Tp[RI],"_GO.csv"))
}

GOterm <- read.csv("C:/Users/monom/Documents/Tb_Goterm_list.csv",header = TRUE)
for(RI in 1:length(Tp)){
  XX <- read.csv(paste0("C:/Users/monom/Documents/DE_2/",Tp[RI],"_GO.csv"),header = TRUE)
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
  write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_2/",Tp[RI],"_GOterm.csv"))
}

for(RI in 1:length(Tp)){
  XX <- read.csv(paste0("C:/Users/monom/Documents/DE_2/",Tp[RI],"_GOterm.csv"),header = TRUE)
  XX$symbol <- rep(NA,nrow(XX))
  for(QI in 1:nrow(XX)){
    Hit <- which(XX$Gb_p50[QI] == BmoriID_s$Gene)
    if(length(Hit)==0){
      next
    }else{
      XX$symbol[QI] <- BmoriID_s$Fb_Symbol[Hit[1]]
    }
    
  }
  write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_2/",Tp[RI],"_symbol.csv"))
}

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



for(RI in 1:8){
  lrt.treat <- glmLRT(fit, coef= as.numeric(RI))
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
  write.csv(x=dat_2,file = paste0("C:/Users/monom/Documents/DE_3/",RI,"_count.csv"))
  write.csv(x=tab,file = paste0("C:/Users/monom/Documents/DE_3/",RI,"_tab.csv"))
  write.csv(x=Y,file = paste0("C:/Users/monom/Documents/DE_3/",RI,"_tab_2.csv"))
}

Ortholog_GO_D <- read.csv("C:/Users/monom/Documents/Ortholog_GO_Drosophila.csv",header = TRUE)
BmoriID_s <- read.csv("C:/Users/monom/Documents/Symbol_BmoriID.csv",header = TRUE)

for(RI in 1:8){
  X <- read.csv(paste0("C:/Users/monom/Documents/DE_3/",RI,"_tab_2.csv"),header = TRUE)
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
  write.csv(x=X,file = paste0("C:/Users/monom/Documents/DE_3/",RI,"_GO.csv"))
}

GOterm <- read.csv("C:/Users/monom/Documents/Tb_Goterm_list.csv",header = TRUE)
for(RI in 1:8){
  XX <- read.csv(paste0("C:/Users/monom/Documents/DE_3/",RI,"_GO.csv"),header = TRUE)
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
  write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_3/",RI,"_GOterm.csv"))
}

for(RI in 1:8){
  XX <- read.csv(paste0("C:/Users/monom/Documents/DE_3/",RI,"_GOterm.csv"),header = TRUE)
  XX$symbol <- rep(NA,nrow(XX))
  for(QI in 1:nrow(XX)){
    Hit <- which(XX$Gb_p50[QI] == BmoriID_s$Gene)
    if(length(Hit)==0){
      next
    }else{
      XX$symbol[QI] <- BmoriID_s$Fb_Symbol[Hit[1]]
    }
    
  }
  write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_3/",RI,"_symbol.csv"))
}


plots <- list()
for(RI in 1:8){
  X <- read.csv(paste0("C:/Users/monom/Documents/DE_3/",RI,"_tab.csv"),header = TRUE)
  X <- mutate(X, thcolor = ifelse(FDR >= 0.05,0,
                                  ifelse(logFC >= 1,1,
                                         ifelse(logFC <= -1,2,
                                                0))))
  for(II in 1:nrow(mel_list)){
    Hit <- which(mel_list$ID[II] == X$X)
    X$thcolor[Hit] <- mel_list$symbol[II]
  }
  
  write.csv(x=X,paste0("C:/Users/monom/Documents/DE_3/",RI,"_FC_signal.csv"))
  
  # Yのデータフレームに、条件に基づいたフィルターを適用
  Y <- X %>%
    filter(!thcolor %in% c("0", "1", "2") & FDR < 0.05) %>%
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
    theme(legend.position = "none", panel.grid = element_blank())+
    labs(x = "log2 (Fold Change)",y = "-log10 (FDR)") +
    theme_bw()
  
  plots[[RI]] <- p
  
}

png("C:/Users/monom/Documents/DE_3/combined_volcano_plot.png", width = 10000, height = 10000, res = 300)
plot_grid(plotlist = plots, ncol = 3)
dev.off()
