tab <- read.csv("C:/Users/monom/Documents/Toll/stage_Toll_pzalue.csv")

star_4 <- tab[which("design_matrixstages4d3:tissueFB" == tab$X),]
star_5_0 <- tab[which("design_matrixtissueFB" == tab$X),]
star_5_3 <- tab[which("design_matrixstages5d3:tissueFB" == tab$X),]
stage <- rbind(star_4,star_5_0,star_5_3)
stage_2 <- stage
stage_2$Estimate[1] <- tab$Estimate[1] + tab$Estimate[2] + tab$Estimate[4] + stage$Estimate[1]
stage_2$Estimate[2] <- tab$Estimate[1] +  tab$Estimate[4]
stage_2$Estimate[3] <- tab$Estimate[1] + tab$Estimate[3] + tab$Estimate[4] + stage$Estimate[3]

star_4_all <- tab[grepl("design_matrixstages4d3:tissueFB" , tab$X),]
star_4_all$Toll <- c("Toll-family", "Toll-10", "Toll-11","Toll-2",
                     "Toll-3","Toll-4","Toll-5","Toll-6","Toll-7","Toll-8","Toll-9")
star_5_0_all <- tab[grepl("design_matrixtissueFB" , tab$X),]
star_5_0_all$Toll <- c("Toll-family", "Toll-10", "Toll-11","Toll-2",
                       "Toll-3","Toll-4","Toll-5","Toll-6","Toll-7","Toll-8","Toll-9")
star_5_3_all <- tab[grepl("design_matrixstages5d3:tissueFB" , tab$X),]
star_5_3_all$Toll <- c("Toll-family", "Toll-10", "Toll-11","Toll-2",
                       "Toll-3","Toll-4","Toll-5","Toll-6","Toll-7","Toll-8","Toll-9")
stage_all <- rbind(star_4_all,star_5_0_all,star_5_3_all)

p <- c(stage[,5],stage_all[,5])
p_a <- p.adjust(p,method = "holm")

stage[,5] <- p_a[1:3]
stage_all[,5] <- p_a[4:length(p_a)]

write.csv(stage,"C:/Users/monom/Documents/Toll/stage_Toll_3.csv")
write.csv(stage_all,"C:/Users/monom/Documents/Toll/stage_Toll_all.csv")

library(ggplot2)


stage_2 <- stage_2 %>%
  mutate(
    X = factor(X, levels = c(
      "design_matrixstages4d3:tissueFB",
      "design_matrixtissueFB",
      "design_matrixstages5d3:tissueFB"
    )),
    X_label = factor(
      X,
      labels = c("day 3 of \n 4th instar larvae", "day 0 of \n 5th instar larvae", "day 3 of \n 5th instar larvae")
    ),
    Absolute_Estimate = abs(Estimate), # 絶対値
    Sign = ifelse(Estimate >= 0, "> 0", "< 0") # 正負を示す列
  )


png(paste0("C:/Users/monom/Documents/Toll/stage_col_2.png")
    ,width=5000, height=5000, res=1000 )
ggplot(stage_2, aes(x = X_label, y = Absolute_Estimate)) +
  geom_col(width = 0.6, color = "black",fill= "#F8F2E9") + # 棒グラフ
  geom_errorbar(
    aes(ymin = Absolute_Estimate - Std..Error, ymax = Absolute_Estimate + Std..Error),
    width = 0.2
  ) + # エラーバー # 色分け
  labs(
    title = "",
    y = "Estimate",
    x="",
    fill="Estimate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 15))# X軸のラベルを傾けて調整
dev.off()




# 必要なパッケージの読み込み
library(ggplot2)
library(dplyr)

# サンプルデータフレームの作成（実際のデータに置き換えてください）
# 例:
# stage_all <- data.frame(
#   X = c("design_matrixstages4d3:tissueFB", "design_matrixtissueFB", "design_matrixstages5d3:tissueFB"),
#   Toll = c("Fatbody", "Toll-2", "Toll-3"),
#   Estimate = c(-1.4704626, 0.4992212, -1.9370492),
#   Std..Error = c(0.5347263, 0.3814260, 0.5517856)
# )

# ラベルリスト
SS <- c("design_matrixstages4d3:tissueFB",
        "design_matrixtissueFB",
        "design_matrixstages5d3:tissueFB")
st <- c("stage4", "stage5_0", "stage5_3")

# 各ステージに「Stage」列を追加
star_4_all <- star_4_all %>%
  mutate(Stage = "stage4")

star_5_0_all <- star_5_0_all %>%
  mutate(Stage = "stage5_0")

star_5_3_all <- star_5_3_all %>%
  mutate(Stage = "stage5_3")

# 全データを結合
stage_all <- bind_rows(star_4_all, star_5_0_all, star_5_3_all)

# データ整形：絶対値と正負で色分け用の列を追加
stage_part <- stage_all %>%
  mutate(
    # Tollの順序を指定: "Fatbody"を先頭に、その後"Toll-2"から"Toll-11"まで順に
    Toll = factor(
      Toll,
      levels = c("Toll-family", "Toll-2", 
                 "Toll-6", "Toll-7", "Toll-8", "Toll-9", "Toll-10", "Toll-11")
    ),
    Absolute_Estimate = abs(Estimate), # 絶対値
    Sign = ifelse(Estimate >= 0, "> 0", "< 0"), # 正負の区別
    Stage = factor(Stage, levels = c("stage4", "stage5_0", "stage5_3")) # ステージ順序
  )

# ステージごとにプロットを生成
for (RI in 1:length(SS)) {
  # PNGの保存設定
  png(
    filename = paste0("C:/Users/monom/Documents/Toll/stage_", st[RI], ".png"),
    width = 5000,
    height = 5000,
    res = 700 # 解像度を調整（1000は非常に高すぎる場合があります）
  )
  
  # フィルタリングされたデータを使用してプロット作成
  plot <- ggplot(
    data = stage_part %>% filter(grepl(SS[RI], X)), 
    aes(x = Toll, y = Absolute_Estimate, fill = Sign)
  ) +
    geom_col(width = 0.6, color = "black") + # 棒グラフ
    geom_errorbar(
      aes(ymin = Absolute_Estimate - Std..Error, ymax = Absolute_Estimate + Std..Error),
      width = 0.2
    ) + # エラーバー
    scale_fill_manual(values = c("> 0" = "#FEB863", "< 0" ="#B2ABD3" )) + # 色分け
    labs(
      title =st[RI],
      x = "Toll",
      y = "Absolute Estimate",
      fill = "estimate"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), # X軸ラベルを傾ける
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 25)# タイトルを中央揃え
    )
  
  # プロットの保存
  print(plot)
  dev.off()
}
