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

CoEf <- c(4,7,8)
Tp <- c("4_3","5_0","5_3")

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



pvalue <- c()
for(RI in 1:8){
  X <- read.csv(paste0("C:/Users/monom/Documents/DE_3/",RI,"_tab.csv"),header = TRUE)
  pvalue <- c(pvalue,X$PValue)
}
##Bonferroni型の補正
pvalue_2 <- p.adjust(pvalue,method="holm")

a <- 1
for(RI in 1:8){
  X <- read.csv(paste0("C:/Users/monom/Documents/DE_3/",RI,"_tab.csv"),header = TRUE)
  X$PValue <- pvalue_2[as.numeric(a):as.numeric(a+26599)]
  Y <- X[which(X$FDR < 0.05 & X$PValue < 0.05),]
  dat$NAlist <- rep(NA,nrow(dat))
  for(II in 1:nrow(Y)){
    Hit <- which(Y$X[II] == rownames(dat))
    if(length(Hit)==0){
      next
    }else{
      dat$NAlist[Hit] <- "Hit"
    }}
  dat_2 <- subset(dat,subset=!is.na(NAlist))
  write.csv(x=dat_2,file = paste0("C:/Users/monom/Documents/DE_4/",RI,"_count.csv"))
  write.csv(x=X,file = paste0("C:/Users/monom/Documents/DE_4/",RI,"_tab.csv"))
  write.csv(x=Y,file = paste0("C:/Users/monom/Documents/DE_4/",RI,"_tab_2.csv"))
  a <- a+26600
}

Ortholog_GO_D <- read.csv("C:/Users/monom/Documents/Ortholog_GO_Drosophila.csv",header = TRUE)
BmoriID_s <- read.csv("C:/Users/monom/Documents/Symbol_BmoriID.csv",header = TRUE)

for(RI in 1:8){
  X <- read.csv(paste0("C:/Users/monom/Documents/DE_4/",RI,"_tab_2.csv"),header = TRUE)
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
  write.csv(x=X,file = paste0("C:/Users/monom/Documents/DE_4/",RI,"_GO.csv"))
}

GOterm <- read.csv("C:/Users/monom/Documents/Tb_Goterm_list.csv",header = TRUE)
for(RI in 1:8){
  XX <- read.csv(paste0("C:/Users/monom/Documents/DE_4/",RI,"_GO.csv"),header = TRUE)
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
  write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_4/",RI,"_GOterm.csv"))
}

for(RI in 1:8){
  XX <- read.csv(paste0("C:/Users/monom/Documents/DE_4/",RI,"_GOterm.csv"),header = TRUE)
  XX$symbol <- rep(NA,nrow(XX))
  for(QI in 1:nrow(XX)){
    Hit <- which(XX$Gb_p50[QI] == BmoriID_s$Gene)
    if(length(Hit)==0){
      next
    }else{
      XX$symbol[QI] <- BmoriID_s$Fb_Symbol[Hit[1]]
    }
    
  }
  write.csv(x=XX,file = paste0("C:/Users/monom/Documents/DE_4/",RI,"_symbol.csv"))
}

library("ggrepel")
library("cowplot")
mel_list <- data.frame(symbol = c("apontic-like","TH","Ddc","lac-2",
                                  "yellow-b", "yellow-d2", "yellow-f2", "yellow-h", "yellow-e","yellow-c",
                                  "cactus"),
                       ID = c("XM_021347184.2","NM_001145322.1","XM_004930959.4","XM_038013104.1",
                              "XM_038012239.1","NM_001043957.1","NM_001043959.1","XM_012690108.3","XM_038016621.1","NM_001043961.1",
                              "NM_001172720.1"))

  EV <- c(4,5)
  plots <- list()
  for(RI in 1:length(EV)){
    X <- read.csv(paste0("C:/Users/monom/Documents/DE_4/",EV[RI],"_tab.csv"),header = TRUE)
    X <- mutate(X, thcolor = ifelse(FDR >= 0.05,0,
                                    ifelse(logFC >= 1,1,
                                           ifelse(logFC <= -1,2,
                                                  0))))
    for(II in 1:nrow(mel_list)){
      Hit <- which(mel_list$ID[II] == X$X)
      X$thcolor[Hit] <- mel_list$symbol[II]
    }
    
    write.csv(x=X,paste0("C:/Users/monom/Documents/DE_4/",RI,"_FC_signal.csv"))
    
    # Yのデータフレームに、条件に基づいたフィルターを適用
    Y <- X %>%
      filter(!thcolor %in% c("0", "1", "2") & FDR < 0.05  & PValue < 0.05) %>%
      mutate(FDR_2 = FDR)  # -log10(FDR)を新しい列として追加
    
    png(paste0("C:/Users/monom/Documents/DE_4/EP_FB_",EV[RI],".png"), width = 2000, height = 1500, res = 300)
    p <- ggplot()+
      geom_point(data=X,mapping =aes(x = logFC, y = -log10(FDR)),color = "gray88",alpha=0.5 )+
      geom_point(data = X %>% filter(thcolor == "2"),mapping = aes(x = logFC, y = -log10(FDR)),color = "#B2ABD3",alpha=0.3,size = 2.0)+
      geom_point(data = X %>% filter(thcolor == "1"),mapping = aes(x = logFC, y = -log10(FDR)),color = "#FEB863",alpha=0.3,size = 2.0)+
      geom_point(data = Y,mapping = aes(x = logFC, y = -log10(FDR)),color = "red",size = 1.0)+
      geom_text_repel(data = Y,aes(x = logFC, y = -log10(FDR),label = thcolor),size = 5,max.overlaps = 20, box.padding = 0.5,vjust = -1, color = "black")+
      geom_hline(yintercept = -log10(0.05),size = 0.2,color = "dark green") +
      geom_vline(xintercept = -1 ,size = 0.2,color = "dark green") +
      geom_vline(xintercept = 1 ,size = 0.2,color = "dark green") +
      theme(legend.position = "none", panel.grid = element_blank(),text = element_text(size = 25))+
      labs(x = "log2 (Fold Change)",y = "-log10 (FDR)") +
      theme_bw()
   plots[[RI]] <- p
    
  
    
  }

png("C:/Users/monom/Documents/DE_4/EP_FB_Time.png", width = 2000, height = 1500, res = 300)
plot_grid(plotlist = plots, ncol = 2)
dev.off()
