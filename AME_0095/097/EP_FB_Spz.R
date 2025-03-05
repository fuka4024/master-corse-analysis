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

Hit <- grep("XM_012695233.2",rownames(dat))
spz3 <- dat[Hit,]
Hit_2 <- grep("XM_038020186.1",rownames(dat))
ea <- dat[Hit_2,]
colnames(spz3) <- rep(c(rep(c("s4d3", "s5d0", "s5d3"),each=3)),2)
colnames(ea) <- rep(c(rep(c("s4d3", "s5d0", "s5d3"),each=3)),2)
spz_toll_ea <- rbind(Toll8,spz3[,1:9],spz3[,10:18],ea[,1:9],ea[,10:18])


breaks <- seq(-4,9, length.out = 100) # 色の区切りを全体に基づいて作成
colors <- colorRampPalette(c("lightyellow","gold", "darkorange", "darkred"))(length(breaks) - 1)
library(pheatmap)
png(paste0("C:/Users/monom/Documents/spz3.png")
    ,width=1000, height=500, res=150 )
par(mar = c(3,3,3,3))
par(oma = c(10,10,10,10))
pheatmap(as.matrix(spz_toll_ea),
         color = colors,
         breaks = breaks,
         cellwidth = 15, cellheight = 12,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_row = c(2,4),
         scale = "none",
         display_numbers = FALSE)
dev.off()