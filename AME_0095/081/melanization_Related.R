##NO.0 install.packages##
install.packages("pheatmap")
install.packages("dplyr")
install.packages("ggrepel")
install.packages("cowplot")

##call packages##
library(edgeR)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)

##set directory
setwd("C:/Users/monom/Documents/")

##count matrix##
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
dat <- as.data.frame(cpm(FB_EP_2,normalized.lib.sizes = TRUE,log = TRUE))
mel_list <- data.frame(symbol = c("apontic-like","TH","Ddc","lac-2","Tak1",
                                  "yellow-b", "yellow-d2", "yellow-f2", "yellow-h", "yellow-e","yellow-c",
                                  "IKKB","cactus"),
                       ID = c("XM_021347184.2","NM_001145322.1","XM_004930959.4","XM_038013104.1","XM_021351619.2",
                              "XM_038012239.1","NM_001043957.1","NM_001043959.1","XM_012690108.3","XM_038016621.1","NM_001043961.1",
                              "XM_004928112.2","NM_001172720.1"))

dat$NAlist <- rep(NA,nrow(dat))
for(RI in 1:nrow(mel_list)){
  Hit <- which(mel_list$ID[RI]==rownames(dat))
  if(length(Hit)==0){
    next
  }else{
    dat$NAlist[Hit] <- mel_list$symbol[RI] 
  }
}

dat_2 <- subset(dat,subset=!is.na(NAlist))
dat_3 <- dat_2
rownames(dat_3) <- dat_2$NAlist
dat_4 <- dat_3[,-19]

s4d3 <- dat_4[, grepl("4_3", colnames(dat_4))]
s5d0 <- dat_4[, grepl("5_0", colnames(dat_4))]
s5d3 <- dat_4[, grepl("5_3", colnames(dat_4))]

dat_5 <- cbind(s4d3,s5d0,s5d3)
colnames(dat_5) <- c(rep("EP_4_3", 3), rep("FB_4_3", 3), rep("EP_5_0", 3),
                       rep("FB_5_0", 3), rep("EP_5_3", 3), rep("FB_5_3", 3))

breaks <- seq(-4,10, length.out = 100) # 色の区切りを全体に基づいて作成
colors <- colorRampPalette(c("lightyellow","gold", "darkorange", "darkred"))(length(breaks) - 1)
png(paste0("C:/Users/monom/Documents/DE/melanization_related.png")
    ,width=5000, height=5000, res=700 )
par(mar = c(3,3,3,3))
par(oma = c(10,10,10,10))
pheatmap(as.matrix(dat_5),
         color = colors,
         breaks = breaks,
         cellwidth = 15, cellheight = 12,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         gaps_col = c(6,12),
         display_numbers = FALSE)
dev.off() 

pheatmap(as.matrix(FB_EP_7),
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "PuOr")))(100),
         cellwidth = 15, cellheight = 12,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = 9,
         scale = "none",
         display_numbers = FALSE)


library(ggplot2)
library(dplyr)
RI <- 1
W <- read.csv(paste0("C:/Users/monom/Documents/DE/",instar_day[RI],"_tab.csv"),header = TRUE)
W$Group <- rep("4_3",nrow(W))
RI <- 2
Y <- read.csv(paste0("C:/Users/monom/Documents/DE/",instar_day[RI],"_tab.csv"),header = TRUE)
Y$Group <- rep("5_0",nrow(W))
RI <- 3
Z <- read.csv(paste0("C:/Users/monom/Documents/DE/",instar_day[RI],"_tab.csv"),header = TRUE)
Z$Group <- rep("5_3",nrow(W))
X <- rbind(W,Y,Z)

plots <- list()
for(RI in 1:length(instar_day)){
  X <- read.csv(paste0("C:/Users/monom/Documents/DE/",instar_day[RI],"_tab.csv"),header = TRUE)


  X <- mutate(X, thcolor = ifelse(FDR >= 0.05,0,
                                  ifelse(logFC >= 1,1,
                                         ifelse(logFC <= -1,2,
                                                0))))
  for(II in 1:nrow(mel_list)){
    Hit <- which(mel_list$ID[II] == X$X)
    X$thcolor[Hit] <- mel_list$symbol[II]
  }

  write.csv(x=X,paste0("C:/Users/monom/Documents/DE/FC_signal.csv"))
  
  Y <- X %>%
    filter(!thcolor %in% c("0", "1", "2") & FDR < 0.05) %>%
    mutate(FDR_2 = FDR)  

 p <- ggplot(X,aes(x = logFC, y = -log10(FDR)))+
    geom_point(data=X,mapping =aes(x = logFC, y = -log10(FDR)),color = "gray88",alpha=0.5 )+
    geom_point(data = X %>% filter(thcolor == "2"),mapping = aes(x = logFC, y = -log10(FDR)),color = "lightsteelblue",alpha=0.3,size = 2.0)+
    geom_point(data = X %>% filter(thcolor == "1"),mapping = aes(x = logFC, y = -log10(FDR)),color = "lightsalmon",alpha=0.3,size = 2.0)+
    geom_point(data = Y,mapping = aes(x = logFC, y = -log10(FDR)),color = "red",size = 1.0)+
    geom_text_repel(data = Y,aes(x = logFC, y = -log10(FDR),label = thcolor),size = 4,max.overlaps = 30, box.padding = 1.5,vjust = -1, size = 3, color = "black")+
    geom_hline(yintercept = -log10(0.05),size = 0.2,color = "dark green") +
    geom_vline(xintercept = -1 ,size = 0.2,color = "dark green") +
    geom_vline(xintercept = 1 ,size = 0.2,color = "dark green") +
    theme(legend.position = "none", panel.grid = element_blank())+
    labs(x = "log2 (Fold Change)",y = "-log10 (FDR)") +
    theme_bw()
 plots[[RI]] <- p
}

  


png("C:/Users/monom/Documents/DE/combined_volcano_plot.png", width = 3000, height = 1500, res = 300)
plot_grid(plotlist = plots, ncol = 3)
dev.off()



