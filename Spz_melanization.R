d1 <- matrix(c(0, 7, 9, 2), nrow=2, byrow=T)
fisher.test(d1)

# データフレームの準備
data <- matrix(c(78,0,22,100), nrow = 2, byrow = TRUE,
               dimnames = list(c("(+)","(-)"), c("pPIG-A3G-Spätzle3","pPIG-A3GR")))

# 棒グラフの作成
png("C:/Users/monom/Documents/melanization.png", width = 5000, height = 5000, res = 500)
par(mar = c(10, 10, 8, 6))
barplot_data <- barplot(
  data,
  beside = FALSE,                # 群ごとに並べて表示
  col = c("#FEB863", "#B2ABD3"), # カラー
  ylim = c(0, max(data) + 3),   # y軸の上限を調整
  legend.text = rownames(data), # 凡例を追加
  bty = "n",
  args.legend = list(title = "melanization", x = "topright",bty = "n",cex=2),
  main = "",
  xlab = "plasmid",
  ylab = "個体数(%)",
  cex.lab = 2,
  cex.axis = 2,
  width = 1,
  space = c(0.1, 0.5),
  xlim = c(0,5)
)

dev.off()
