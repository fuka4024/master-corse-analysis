library(ggplot2)

# 条件リストを作成
conditions <- list(Area = Area, Area_2 = Area_2, Area_3 = Area_3, Area_4 = Area_4)

# 結果を保存するリスト
results <- list()

# 1000回のリサンプリングとウィルコクソン検定
n_resamples <- 1000
for (condition_name in names(conditions)) {
  data <- conditions[[condition_name]]
  n_sample <- floor(length(data) * 0.6) # サンプル数の6割
  
  res <- replicate(n_resamples, {
    sample_data <- sample(data, size = n_sample, replace = TRUE)
    test <- wilcox.test(sample_data, data)
    c(w = test$statistic, p = test$p.value)
  })
  
  results[[condition_name]] <- t(res) # 転置して行:サンプル、列:wとp値
}

# ヒストグラムを描画
plot_histograms <- function(data, title) {
  data <- as.data.frame(data)
  colnames(data) <- c("w", "p_value")
  
  # wのヒストグラム
  p1 <- ggplot(data, aes(x = w)) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
    labs(title = paste(title, "- W (Test Statistic)"), x = "W", y = "Frequency") +
    theme_minimal()
  
  # p値のヒストグラム
  p2 <- ggplot(data, aes(x = p_value)) +
    geom_histogram(bins = 30, fill = "red", alpha = 0.7) +
    labs(title = paste(title, "- P-values"), x = "P-value", y = "Frequency") +
    theme_minimal()
  
  list(w_plot = p1, p_value_plot = p2)
}

# ヒストグラムの作成と表示
histograms <- lapply(names(results), function(name) {
  plot_histograms(results[[name]], name)
})

# 結果をプロット
for (RI in 1:4) {
  png(paste0("C:/Users/monom/Documents/wilcox/w_",RI,".png")
      ,width=2000, height=2000, res=300 )
  print(histograms[[RI]]$w_plot)
  dev.off()
  
  png(paste0("C:/Users/monom/Documents/wilcox/p_",RI,".png")
      ,width=2000, height=2000, res=300 )
  print(histograms[[RI]]$p_value_plot)
  dev.off()
}



# 実際のWの値（ウィルコクソン検定の結果から）
actual_w <- list(
  Area_vs_Area_2 = 641659,
  Area_3_vs_Area_2 = 110379,
  Area_4_vs_Area_2 = 147053
)

# 各リサンプリング結果のW分布
null_distributions <- list(
  Area = results$Area[, "w.W"],
  Area_2 = results$Area_2[, "w.W"],
  Area_3 = results$Area_3[, "w.W"],
  Area_4 = results$Area_4[, "w.W"]
)

# 実験条件ごとのリサンプリング比較マップ
comparison_map <- list(
  Area_vs_Area_2 = c("Area", "Area_2"),
  Area_3_vs_Area_2 = c("Area_3", "Area_2"),
  Area_4_vs_Area_2 = c("Area_4", "Area_2")
)

# プロットを生成する関数
plot_actual_vs_null <- function(actual_w, null_w, title, condition_name) {
  ggplot() +
    geom_histogram(aes(x = null_w), bins = 30, fill = "blue", alpha = 0.7) +
    geom_vline(aes(xintercept = actual_w), color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste0(title, "compared to", condition_name,"'"),
      x = "W Value",
      y = "Frequency"
    ) +
    theme_minimal()
}

# プロットの作成
QI <- 1
for (RI in 1:length(comparison_map)) {
  actual_value <- actual_w[[RI]]  # 実際のWの値
  conditions <- comparison_map[[RI]]  # リサンプリング対象の条件リスト
  for (II in 1:length(conditions)) {
    Hit <- which(names(null_distributions) == conditions[II])
    null_w <- null_distributions[[Hit]]
    # 対応する条件のリサンプリング結果
    png(
      filename = paste0("C:/Users/monom/Documents/wilcox/",QI, "_2.png"),
      width = 2000, height = 2000, res = 300
    )
    a <- plot_actual_vs_null(
      actual_value, null_w, title = conditions[II], condition_name = conditions[II]
    )
    print(a)
    dev.off()
    QI <- QI + 1
  }
}


for (RI in 1:length(plots)) {
  png(
    filename = paste0("C:/Users/monom/Documents/wilcox/",RI, "_2.png"),
    width = 2000, height = 2000, res = 300
  )
  print(plots[[RI]])  # 対応するプロットを描画
  dev.off()
}
