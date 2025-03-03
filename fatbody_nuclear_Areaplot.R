# 必要なライブラリ
library(ggplot2)
library(viridis)
library(png)
library(grid)
library(RColorBrewer)
library(dplyr)

# 画像の幅と高さを指定
image_width <- 1024  # 画像の幅
image_height <- 1024  # 画像の高さ

# 二値化画像を読み込む
binary_image <- readPNG("fatbody_DAPI_6/merge_3.png")  # 二値化画像ファイル名を指定
data <- read.csv("fatbody_DAPI_6/Results_3_p.csv")
data$Y <- image_height - data$Y


# 背景画像をプロット
p <- ggplot(data, aes(x = X, y = Y)) +
  # 背景として二値化画像を追加
  annotation_custom(
    rasterGrob(binary_image, 
               width = unit(1, "npc"), 
               height = unit(1, "npc")),
    xmin = 0, xmax = image_width,
    ymin = 0, ymax = image_height
  ) +
  # カラープロットを重ねる
  geom_point(aes(color = Area), size = 3) +
  scale_color_viridis_c(option = "turbo")  +  # カラースケール
  labs(title = "",
       x = "",
       y = "",
       color = "Area") +
  coord_fixed(ratio = 1, xlim = c(0, image_width), ylim = c(0, image_height)) +  # 座標範囲を設定
  theme_minimal() +  # シンプルなテーマ
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # タイトルのスタイル
    axis.title = element_text(size = 14),  # 軸ラベルのスタイル
    axis.text = element_text(size = 12)    # 軸目盛りのスタイル
  )



png(
  filename = paste0("C:/Users/monom/Documents/fatbody_DAPI_6/overlay_plot_3.png"),
  width = 2000, height = 2000, res = 300
)
print(p)
dev.off()



