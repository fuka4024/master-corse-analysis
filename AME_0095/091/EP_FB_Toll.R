FB_EP <-read.table("FB_EP.txt",head=TRUE)
rownames(FB_EP) <- FB_EP[,1]
FB_EP <- FB_EP[,-1]
colnames(FB_EP) <- c("EP_4_3_1","EP_4_3_2","EP_4_3_3",
                     "EP_5_0_1","EP_5_0_2","EP_5_0_3",
                     "EP_5_3_1","EP_5_3_2","EP_5_3_3",
                     "FB_4_3_1","FB_4_3_2","FB_4_3_3",
                     "FB_5_0_1","FB_5_0_2","FB_5_0_3",
                     "FB_5_3_1","FB_5_3_2","FB_5_3_3")

Toll_list <- data.frame(symbol = c("BmToll-1","BmToll-2","BmToll-3","BmToll-4","BmToll-5",
                                   "BmToll-6","BmToll-7","BmToll-8","BmToll-9","BmToll-10","Bmtoll-11"),
                        ID = c("NM_001123349.1","XM_004921670.4","XM_038015525.1","XM_004927259.4","XM_021348482.2",
                               "XM_012697885.3","XM_004921675.3","XM_004921685.4","XM_012691451.3","XM_004921681.4","XM_004921682.4"))
Toll_family <- c("BmToll-1","BmToll-2",
                 "BmToll-6","BmToll-7","BmToll-8","BmToll-9","BmToll-10","Bmtoll-11")

library(edgeR)
FB_EP <- as.data.frame(cpm(FB_EP,normalized.lib.sizes = TRUE,log=TRUE))
FB_EP_2 <- FB_EP[rownames(FB_EP) %in% c("NM_001123349.1","XM_004921670.4",
                                        "XM_012697885.3","XM_004921675.3","XM_004921685.4","XM_012691451.3",
                                        "XM_004921681.4","XM_004921682.4"), ]

FB_EP_2$num <- rep(NA,nrow(FB_EP_2))
for(RI in 1:nrow(Toll_list)){
  Hit <- which(Toll_list$ID[RI]==rownames(FB_EP_2))
  FB_EP_2$num[Hit] <- RI
  
}

FB_EP_3 <- FB_EP_2[order(FB_EP_2$num),]


FB_EP_4 <- FB_EP_3[,10:18]

dat <- data.frame(matrix(rep(NA,180*5),nrow = 180,ncol=5))
colnames(dat) <- c("tissue", "rep", "Toll", "stage","count")


star_4 <- c()
star_5_0 <- c()
star_5_3 <- c()
for(RI in 1:6){
  star_4 <- c(star_4,FB_EP_3[, grepl("4_3", colnames(FB_EP_3))][,RI])
  star_5_0 <- c(star_5_0,FB_EP_3[, grepl("5_0", colnames(FB_EP_3))][,RI])
  star_5_3 <- c(star_5_3,FB_EP_3[, grepl("5_3", colnames(FB_EP_3))][,RI])
}
dat$rep   <- as.factor(rep(rep(c(1,2,3),each=10),6))
dat$tissue <- relevel(as.factor(rep(rep(c("EP", "FB"), each = 30),3)),ref="EP")
dat$Toll <- as.factor(rep(Toll_family,18))
dat$stage <- relevel(as.factor(c(rep("s4d3",length(star_4)),rep("s5d0",length(star_5_0)),rep("s5d3",length(star_5_3)))),ref="s5d0")
dat$count <- round(c(star_4,star_5_0,star_5_3))

# 必要なパッケージを読み込み
library(ggplot2)
library(dplyr)

# データの準備
# レプリケイトごとの平均値を計算
data_avg <- dat %>%
  group_by(tissue, stage, Toll) %>%
  summarise(mean_count = mean(count, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(stage = factor(stage, levels = c("s4d3", "s5d0", "s5d3")))


# EPのプロット
png("C:/Users/monom/Documents/Toll/EP_Toll.png", width = 5000, height = 5000, res = 700)

plot_EP <- ggplot(data_avg %>% filter(tissue == "EP") , aes(x = stage, y = mean_count, group = Toll, color = Toll)) +
  geom_line(linewidth = 5) + # 線の太さを調整
  geom_point(size = 2) + # 点を追加
  scale_color_brewer(palette = "Set3") + # 色のセットを変更可能
  labs(title = "EP Tissue", x = "Stage", y = "Average CPM", color = "Toll") +
  theme(legend.position = "none", panel.grid = element_blank(),text = element_text(size = 100))+
  theme_bw()

print(plot_EP)

dev.off()

# FBのプロット
png("C:/Users/monom/Documents/Toll/FB_Toll.png", width = 5000, height = 5000, res = 700)
plot_FB <- ggplot(data_avg %>% filter(tissue == "FB"), aes(x = stage, y = mean_count, group = Toll, color = Toll)) +
  geom_line(linewidth = 5) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set3") +
  labs(title = "FB Tissue", x = "Stage", y = "Average CPM", color = "Toll") +
  theme(legend.position = "none", panel.grid = element_blank(),text = element_text(size = 100))+
  theme_bw()

print(plot_FB)

dev.off()
# プロットを表示
print(plot_EP)
print(plot_FB)





# データの準備
# heatmap_data は元のデータフレーム
rownames(FB_EP_4) <- Toll_family  # 遺伝子名を行名にする
colnames(FB_EP_4) <- c(rep("FB_4_3", 3), rep("FB_5_0", 3), rep("FB_5_3", 3))

library(pheatmap)
library(RColorBrewer)
# Heatmap オブジェクトを作成
png(paste0("C:/Users/monom/Documents/Toll/heatmap.png")
    ,width=5000, height=5000, res=1000 )
par(mar = c(3,3,3,3))
par(oma = c(10,10,10,10))
pheatmap(as.matrix(FB_EP_4),
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "PuOr")))(100),
         cellwidth = 15, cellheight = 12,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "none",
         display_numbers = FALSE)
dev.off()

FB_EP_5 <- FB_EP[rownames(FB_EP) %in% c("NM_001123349.1","XM_004921670.4","XM_038015525.1","XM_004927259.4","XM_021348482.2",
                                        "XM_012697885.3","XM_004921675.3","XM_004921685.4","XM_012691451.3","XM_004921681.4","XM_004921682.4"), ]

FB_EP_5$num <- rep(NA,nrow(FB_EP_5))
for(RI in 1:nrow(Toll_list)){
  Hit <- which(Toll_list$ID[RI]==rownames(FB_EP_5))
  FB_EP_5$num[Hit] <- RI
  
}

FB_EP_5 <- FB_EP_5[order(FB_EP_5$num),1:18]

rownames(FB_EP_5) <- Toll_list$symbol  # 遺伝子名を行名にする
colnames(FB_EP_5) <- c(rep("EP_4_3", 3), rep("EP_5_0", 3), rep("EP_5_3", 3),
                       rep("FB_4_3", 3), rep("FB_5_0", 3), rep("FB_5_3", 3))
png(paste0("C:/Users/monom/Documents/Toll/heatmapEP_FB.png")
    ,width=10000, height=5000, res=1000 )
par(mar = c(3,3,3,3))
par(oma = c(10,10,10,10))
pheatmap(as.matrix(FB_EP_5[,1:18]),
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "PuOr")))(100),
         cellwidth = 15, cellheight = 12,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "none",
         display_numbers = FALSE)
dev.off()


FB_EP_6 <- FB_EP_3[,1:9]
rownames(FB_EP_6) <- Toll_family  # 遺伝子名を行名にする
colnames(FB_EP_6) <- c(rep("EP_4_3", 3), rep("EP_5_0", 3), rep("EP_5_3", 3))

png(paste0("C:/Users/monom/Documents/Toll/heatmap_EP.png")
    ,width=5000, height=5000, res=1000 )
par(mar = c(3,3,3,3))
par(oma = c(10,10,10,10))
pheatmap(as.matrix(FB_EP_6),
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "PuOr")))(100),
         cellwidth = 15, cellheight = 12,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "none",
         display_numbers = FALSE)
dev.off()

FB_EP_7 <- FB_EP_3[,1:18]
rownames(FB_EP_7) <- Toll_family  # 遺伝子名を行名にする
colnames(FB_EP_7) <- c(rep("EP_4_3", 3), rep("EP_5_0", 3), rep("EP_5_3", 3),
                       rep("FB_4_3", 3), rep("FB_5_0", 3), rep("FB_5_3", 3))

png(paste0("C:/Users/monom/Documents/Toll/heatmap_EP_FB_gap.png")
    ,width=5000, height=5000, res=1000 )
par(mar = c(3,3,3,3))
par(oma = c(10,10,10,10))
pheatmap(as.matrix(FB_EP_7),
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "PuOr")))(100),
         cellwidth = 15, cellheight = 12,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = 9,
         scale = "none",
         display_numbers = FALSE)
dev.off()



coef_data <- f$coefficients

# 信頼区間の計算
conf_int <- confint(m1)

# ggplotを使用したグラフ作成
ggplot(coef_data, aes(x = stage, y = Estimate, ymin = conf_int[,1], ymax = conf_int[,2])) +
  geom_col(fill = "skyblue") +
  geom_errorbar(width = 0.2, color = "blue") +
  labs(title = "Stage Effect on Response",
       x = "Stage",
       y = "Coefficient Estimate") +
  theme_minimal()

