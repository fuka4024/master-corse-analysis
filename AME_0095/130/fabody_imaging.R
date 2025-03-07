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
binary_image <- readPNG("merge.png")  # 二値化画像ファイル名を指定
data <- read.csv("Results_test.csv")
data$Y <- image_height - data$Y

breaks <- seq(-4,9, length.out = 100)

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

# プロットを保存
ggsave("C:/Users/monom/Documents/overlay_plot.png", width = 12, height = 8, dpi = 300)  # 保存時にサイズと解像度を指定

png(
  filename = paste0("C:/Users/monom/Documents/overlay_plot.png"),
  width = 2000, height = 2000, res = 300
)
print(p)
dev.off()

##
data <- read.csv("fatbody_DAPI_6/Results_1_p.csv")
data_2 <- read.csv("fatbody_DAPI_6/Results_1_m.csv")
data$Group <- rep("Non",nrow(data))
Spz_num <- c(32,33,41,57,61,55,65,70,74,76,77,94,101,120,128,129,122,117,105,79,80,83,87,91,104,121,119,141,144,
             162,179,186,176,189,219,224,241,238,100,125,132,138,131,152,147,148,157,177,178,170,154,170,194,199)

data$Group[Spz_num] <- rep("Spz",length(Spz_num))
data$Group_2 <- unlist(lapply( 1:nrow(data),function(x){paste0(data$Group[x],x)}))
data$Xm <- data_2$X
data$Ym <- data_2$Y
data$Area <- data_2$Area


group_Spz <- data %>% filter(Group == "Spz")
group_Non <- data %>% filter(Group == "Non")  


dat <- data.frame(matrix(rep(NA,nrow(group_Spz)*nrow(group_Non)*12),nrow = nrow(group_Spz)*nrow(group_Non),ncol=12))
colnames(dat) <- c("Spz","Non","Spz_X","Spz_Y","Non_X","Non_Y","Spz_Xm","Spz_Ym","Non_Xm","Non_Ym","distance","Non_Area")
dat$Spz <- rep(group_Spz$Group_2,nrow(group_Non))
dat$Non <- rep(group_Non$Group_2,each = nrow(group_Spz))
dat$Non_Area <- rep(group_Non$Area,each = nrow(group_Spz))
dat$Spz_X <- rep(group_Spz$X,nrow(group_Non))
dat$Spz_Y <- rep(group_Spz$Y,nrow(group_Non))
dat$Non_X <- rep(group_Non$X,each = nrow(group_Spz))
dat$Non_Y <- rep(group_Non$Y,each = nrow(group_Spz))
dat$Spz_Xm <- rep(group_Spz$Xm,nrow(group_Non))
dat$Spz_Ym <- rep(group_Spz$Ym,nrow(group_Non))
dat$Non_Xm <- rep(group_Non$Xm,each = nrow(group_Spz))
dat$Non_Ym <- rep(group_Non$Ym,each = nrow(group_Spz))


for (II in 1:nrow(dat)) {
    dat$distance[II] <- sqrt((dat$Spz_Xm[II] - dat$Non_Xm[II])^2 + (dat$Spz_Ym[II] - dat$Non_Ym[II])^2)
}

dat_sort <- dat[sort(as.numeric(dat$distance),decreasing = F,index= T)$ix,]
dat_top <- dat_sort[match(as.character(group_Non$Group_2),as.character(dat_sort$Non)),]
rownames(dat_top) <- dat_top$Non
dat_best_rows <- unlist(lapply(as.character(group_Non$Group_2),
                                        function(x) intersect(which(dat$Non==x),
                                                              which(dat$distance==dat_top[x,"distance"])
                                        )))
dat_best <- dat[dat_best_rows,]

dat_2 <- data.frame(matrix(rep(NA,nrow(group_Spz)*2),nrow = nrow(group_Spz),ncol=2))
colnames(dat_2) <- c("distance","Non_Area")
dat_2$distance <- rep(0,nrow(dat_2))
dat_2$Non_Area <- group_Spz$Area
dat_3 <- rbind(dat_best[,11:12],dat_2)
A <- sort(as.numeric(dat_3$Non_Area),decreasing = T,index= F)
N <- mean(A[239:319])
AN <- dat_3$Non_Area/N
test <- dat_3
dat_3$Non_Area <- AN

ggplot(dat_3, aes(x = distance, y = Non_Area)) +
  geom_point( size = 3) +
  labs(title = "",
       x = "",
       y = "") +
  theme_minimal() +  # シンプルなテーマ
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # タイトルのスタイル
    axis.title = element_text(size = 14),  # 軸ラベルのスタイル
    axis.text = element_text(size = 12)    # 軸目盛りのスタイル
  )






##
data2 <- read.csv("fatbody_DAPI_6/Results_2_p.csv")
data2_2 <- read.csv("fatbody_DAPI_6/Results_2_m.csv")
data2$Group <- rep("Non",nrow(data2))
Spz_num <- c(47,54,53,66,56,78,70,84,81,91,103,118,124,144,165,173,186,195,218,227,245,269,279,302,326,339,361,372,
             392,408,426,447,63,62,65,76,102,126,175,202,224,248,256,281,317,337,350,367,374,398,418,441,450,
             94,93,96,82,98,99,110,116,123,135,157,166,179,149,150,169,142,208,194,187,198,188,184,216,237,258,289,285,307,
             323,356,363,387,388,385,362,378,399,404,413,428,440,449,435,419,420,405,397)

data2$Group[Spz_num] <- rep("Spz",length(Spz_num))
data2$Group_2 <- unlist(lapply( 1:nrow(data2),function(x){paste0(data2$Group[x],x)}))
data2$Xm <- data2_2$X
data2$Ym <- data2_2$Y
data2$Area <- data2_2$Area


group_Spz_2 <- data2 %>% filter(Group == "Spz")
group_Non_2 <- data2 %>% filter(Group == "Non")  


dat2 <- data.frame(matrix(rep(NA,nrow(group_Spz_2)*nrow(group_Non_2)*12),nrow = nrow(group_Spz_2)*nrow(group_Non_2),ncol=12))
colnames(dat2) <- c("Spz","Non","Spz_X","Spz_Y","Non_X","Non_Y","Spz_Xm","Spz_Ym","Non_Xm","Non_Ym","distance","Non_Area")
dat2$Spz <- rep(group_Spz_2$Group_2,nrow(group_Non_2))
dat2$Non <- rep(group_Non_2$Group_2,each = nrow(group_Spz_2))
dat2$Non_Area <- rep(group_Non_2$Area,each = nrow(group_Spz_2))
dat2$Spz_X <- rep(group_Spz_2$X,nrow(group_Non_2))
dat2$Spz_Y <- rep(group_Spz_2$Y,nrow(group_Non_2))
dat2$Non_X <- rep(group_Non_2$X,each = nrow(group_Spz_2))
dat2$Non_Y <- rep(group_Non_2$Y,each = nrow(group_Spz_2))
dat2$Spz_Xm <- rep(group_Spz_2$Xm,nrow(group_Non_2))
dat2$Spz_Ym <- rep(group_Spz_2$Ym,nrow(group_Non_2))
dat2$Non_Xm <- rep(group_Non_2$Xm,each = nrow(group_Spz_2))
dat2$Non_Ym <- rep(group_Non_2$Ym,each = nrow(group_Spz_2))


for (II in 1:nrow(dat2)) {
  dat2$distance[II] <- sqrt((dat2$Spz_Xm[II] - dat2$Non_Xm[II])^2 + (dat2$Spz_Ym[II] - dat2$Non_Ym[II])^2)
}

dat2_sort <- dat2[sort(as.numeric(dat2$distance),decreasing = F,index= T)$ix,]
dat2_top <- dat2_sort[match(as.character(group_Non_2$Group_2),as.character(dat2_sort$Non)),]
rownames(dat2_top) <- dat2_top$Non
dat2_best_rows <- unlist(lapply(as.character(group_Non_2$Group_2),
                               function(x) intersect(which(dat2$Non==x),
                                                     which(dat2$distance==dat2_top[x,"distance"])
                               )))
dat2_best <- dat2[dat2_best_rows,]

dat2_2 <- data.frame(matrix(rep(NA,nrow(group_Spz_2)*2),nrow = nrow(group_Spz_2),ncol=2))
colnames(dat2_2) <- c("distance","Non_Area")
dat2_2$distance <- rep(0,nrow(dat2_2))
dat2_2$Non_Area <- group_Spz_2$Area
dat2_3 <- rbind(dat2_best[,11:12],dat2_2)
A <- sort(as.numeric(dat2_3$Non_Area),decreasing = T,index= F)
N <- mean(A[345:length(A)])
AN <- dat2_3$Non_Area/N
dat2_3$Non_Area <- AN

ggplot(dat2_3, aes(x = distance, y = Non_Area)) +
  geom_point( size = 3) +
  labs(title = "",
       x = "",
       y = "") +
  theme_minimal() +  # シンプルなテーマ
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # タイトルのスタイル
    axis.title = element_text(size = 14),  # 軸ラベルのスタイル
    axis.text = element_text(size = 12)    # 軸目盛りのスタイル
  )





##
data3 <- read.csv("fatbody_DAPI_6/Results_3_p.csv")
data3_2 <- read.csv("fatbody_DAPI_6/Results_3_m.csv")
data3$Group <- rep("Non",nrow(data3))
Spz_num <- c(27,35,45,47,67,76,91,116,112,139,125,134,130,124,115,104,92,97,90,74,80,68,61,59,87,105,107,108)

data3$Group[Spz_num] <- rep("Spz",length(Spz_num))
data3$Group_2 <- unlist(lapply( 1:nrow(data3),function(x){paste0(data3$Group[x],x)}))
data3$Xm <- data3_2$X
data3$Ym <- data3_2$Y
data3$Area <- data3_2$Area


group_Spz_3 <- data3 %>% filter(Group == "Spz")
group_Non_3 <- data3 %>% filter(Group == "Non")  



dat3 <- data.frame(matrix(rep(NA,nrow(group_Spz_3)*nrow(group_Non_3)*12),nrow = nrow(group_Spz_3)*nrow(group_Non_3),ncol=12))
colnames(dat3) <- c("Spz","Non","Spz_X","Spz_Y","Non_X","Non_Y","Spz_Xm","Spz_Ym","Non_Xm","Non_Ym","distance","Non_Area")
dat3$Spz <- rep(group_Spz_3$Group_2,nrow(group_Non_3))
dat3$Non <- rep(group_Non_3$Group_2,each = nrow(group_Spz_3))
dat3$Non_Area <- rep(group_Non_3$Area,each = nrow(group_Spz_3))
dat3$Spz_X <- rep(group_Spz_3$X,nrow(group_Non_3))
dat3$Spz_Y <- rep(group_Spz_3$Y,nrow(group_Non_3))
dat3$Non_X <- rep(group_Non_3$X,each = nrow(group_Spz_3))
dat3$Non_Y <- rep(group_Non_3$Y,each = nrow(group_Spz_3))
dat3$Spz_Xm <- rep(group_Spz_3$Xm,nrow(group_Non_3))
dat3$Spz_Ym <- rep(group_Spz_3$Ym,nrow(group_Non_3))
dat3$Non_Xm <- rep(group_Non_3$Xm,each = nrow(group_Spz_3))
dat3$Non_Ym <- rep(group_Non_3$Ym,each = nrow(group_Spz_3))


for (II in 1:nrow(dat3)) {
  dat3$distance[II] <- sqrt((dat3$Spz_Xm[II] - dat3$Non_Xm[II])^2 + (dat3$Spz_Ym[II] - dat3$Non_Ym[II])^2)
}

dat3_sort <- dat3[sort(as.numeric(dat3$distance),decreasing = F,index= T)$ix,]
dat3_top <- dat3_sort[match(as.character(group_Non_3$Group_2),as.character(dat3_sort$Non)),]
rownames(dat3_top) <- dat3_top$Non
dat3_best_rows <- unlist(lapply(as.character(group_Non_3$Group_2),
                                function(x) intersect(which(dat3$Non==x),
                                                      which(dat3$distance==dat3_top[x,"distance"])
                                )))
dat3_best <- dat3[dat3_best_rows,]

dat3_2 <- data.frame(matrix(rep(NA,nrow(group_Spz_3)*2),nrow = nrow(group_Spz_3),ncol=2))
colnames(dat3_2) <- c("distance","Non_Area")
dat3_2$distance <- rep(0,nrow(dat3_2))
dat3_2$Non_Area <- group_Spz_3$Area
dat3_3 <- rbind(dat3_best[,11:12],dat3_2)
A <- sort(as.numeric(dat3_3$Non_Area),decreasing = T,index= F)
N <- mean(A[230:length(A)])
AN <- dat3_3$Non_Area/N
dat3_3$Non_Area <- AN

ggplot(dat3_3, aes(x = distance, y = Non_Area)) +
  geom_point( size = 3) +
  labs(title = "",
       x = "",
       y = "") +
  theme_minimal() +  # シンプルなテーマ
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # タイトルのスタイル
    axis.title = element_text(size = 14),  # 軸ラベルのスタイル
    axis.text = element_text(size = 12)    # 軸目盛りのスタイル
  )


dat_4 <-rbind(dat_3,dat2_3,dat3_3)
dat_7 <- data.frame(matrix(rep(NA,nrow(dat_4)*3),nrow = nrow(dat_4),ncol=3))
colnames(dat_7) <- c("distance","Non_Area","individual")
dat_7[,1:2] <- dat_4
dat_7$individual <- c(rep(1,nrow(dat_3)),rep(2,nrow(dat2_3)),rep(3,nrow(dat3_3)))


# 各種回帰モデルを適用
linear_model <- lm(Non_Area ~ distance, data=dat_4)
quadratic_model <- lm(Non_Area ~ poly(distance, 2), data=dat_4)
exponential_model <- lm(log(Non_Area+1) ~ distance, data=dat_4)  # Non_Areaを対数変換して線形回帰
logarithmic_model <- lm(Non_Area ~ log(distance+1), data=dat_4)

models <- list(linear_model, quadratic_model, exponential_model, logarithmic_model,glm_poisson,glm_gamma)
model_names <- c("Linear", "Quadratic", "Exponential", "Logarithmic","p","g")
aic_values <- sapply(models, AIC)

aic_df <- data.frame(Model = model_names, AIC = aic_values)
print(aic_df)


library(lmerTest)
sm <- summary(m1 <- lmer(Non_Area ~  distance + (1|individual) , data = dat_7))
predicted_values <- predict(m1, re.form = NA)
dat_7$predicted <- predicted_values
dat_7$individual <- as.factor(dat_7$individual)

scale_factor <- 316 / 1024
png(
  filename = paste0("C:/Users/monom/Documents/distance_Area_merge.png"),
  width = 4000, height = 2000, res = 300
)

 p1 <-   # ≈ 0.308 μm / pixel
 ggplot(dat_image_3, aes(x = X * scale_factor, y = Y * scale_factor)) +
   # 背景として二値化画像を追加
   annotation_custom(
     rasterGrob(binary_image,
                width = unit(1, "npc"),
                height = unit(1, "npc")),
     xmin = 0, xmax = image_width * scale_factor,
     ymin = 0, ymax = image_height * scale_factor
   ) +
   # カラープロットを重ねる
   geom_point(aes(color = distance), size = 3) +
   scale_color_viridis_c(option = "turbo")  +  # カラースケール
   labs(title = "",
        x = "μm",
        y = "μm",
        color = "distance") +
   coord_fixed(ratio = 1,
               xlim = c(0, image_width * scale_factor),
               ylim = c(0, image_height * scale_factor)) +  # 座標範囲を設定
   theme_minimal() +  # シンプルなテーマ
   theme(
     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # タイトルのスタイル
     axis.ticks = element_blank(),  # tickの線を消す
     axis.text = element_blank(),   # tickの数字を消す
     axis.title = element_blank(),  # 軸のラベルを消す
     axis.line = element_blank()    # 軸目盛りのスタイル
   )
 
 
p2 <- ggplot(dat_7, aes(distance, Non_Area,colour = distance)) +
  geom_point(size = 3.8,alpha=0.23) +
  geom_line(aes(y = predicted), color = "#934EC2", size = 2) + 
  scale_color_viridis_c(option = "turbo")  +
  labs(title = "",
       x = "distance (μm)",
       y = "Area") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12),    
    legend.position = "right"  
  ) 



# 3. 散布図とカラースケールを結合
plot_grid(p1 + theme(legend.position = "none"),  # 凡例なしの散布図
          p2,  # カラースケールのみ
          nrow = 1, rel_heights = c(4, 4))  # 散布図と凡例の高さ比

dev.off()






png(
  filename = paste0("C:/Users/monom/Documents/distance_Area.png"),
  width = 2000, height = 2000, res = 300
)

ggplot(dat_4, aes(x = distance, y = Non_Area)) +
  geom_point(color = "#377EB8", size = 3) +
  labs(title = "",
       x = "distance(μm)",
       y = "Area") +
  theme_minimal() +  # シンプルなテーマ
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # タイトルのスタイル
    axis.title = element_text(size = 14),  # 軸ラベルのスタイル
    axis.text = element_text(size = 12)    # 軸目盛りのスタイル
  )

dev.off()

dat_5 <- rbind(dat_3,dat2_3,dat3_3)
A <- sort(as.numeric(dat_5$Non_Area),decreasing = T,index= F)
N <- mean(A[814:length(A)])
AN <- dat_5$Non_Area/N
dat_5$Non_Area <- AN

png(
  filename = paste0("C:/Users/monom/Documents/distance_Area_2.png"),
  width = 2000, height = 2000, res = 300
)

  ggplot(dat_5, aes(x = distance, y = Non_Area)) +
  geom_point(color = "#377EB8", size = 3) +
  labs(title = "",
       x = "distance(μm)",
       y = "Area") +
  theme_minimal() +  # シンプルなテーマ
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # タイトルのスタイル
    axis.title = element_text(size = 14),  # 軸ラベルのスタイル
    axis.text = element_text(size = 12)    # 軸目盛りのスタイル
  )



dev.off()

image_width <- 1024  # 画像の幅
image_height <- 1024  # 画像の高さ

# 二値化画像を読み込む
binary_image <- readPNG("fatbody_DAPI_6/merge_1.png")  # 二値化画像ファイル名を指定


dat_image <-  as.data.frame(cbind(dat_best$Non_X,dat_best$Non_Y,dat_best$distance))
colnames(dat_image) <- c("X","Y","distance")
dat_image_2 <- as.data.frame(cbind(group_Spz$X,group_Spz$Y))
dat_image_2$distance <- rep(0,nrow(dat_image_2))
colnames(dat_image_2) <- c("X","Y","distance")
dat_image_3 <- rbind(dat_image,dat_image_2)
dat_image_3$Y <- image_height - dat_image_3$Y


# 背景画像をプロット
png(
  filename = paste0("C:/Users/monom/Documents/distance_scale.png"),
  width = 2000, height = 2000, res = 300
)
ggplot(dat_image_3, aes(x = X, y = Y)) +
  # 背景として二値化画像を追加
  annotation_custom(
    rasterGrob(binary_image, 
               width = unit(1, "npc"), 
               height = unit(1, "npc")),
    xmin = 0, xmax = image_width,
    ymin = 0, ymax = image_height
  ) +
  # カラープロットを重ねる
  geom_point(aes(color = distance), size = 3) +
  scale_color_viridis_c(option = "turbo")  +  # カラースケール
  labs(title = "",
       x = "",
       y = "",
       color = "distance") +
  coord_fixed(ratio = 1, xlim = c(0, image_width), ylim = c(0, image_height)) +  # 座標範囲を設定
  theme_minimal() +  # シンプルなテーマ
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # タイトルのスタイル
    axis.title = element_text(size = 14),  # 軸ラベルのスタイル
    axis.text = element_text(size = 12)    # 軸目盛りのスタイル
  )
dev.off()

##layer0
layer0 <- c()
layer0 <- c(layer0,dat_3$Non_Area[which(dat_3$distance == 0)],dat2_3$Non_Area[which(dat2_3$distance == 0)],
            dat3_3$Non_Area[which(dat3_3$distance == 0)])

##layer1
layer1 <- c()
layer1 <- c(layer1,dat_3$Non_Area[which(dat_3$distance > 0 & dat_3$distance <= 75)],
            dat2_3$Non_Area[which(dat2_3$distance > 0 & dat2_3$distance <= 75)],
            dat3_3$Non_Area[which(dat3_3$distance > 0 & dat3_3$distance <= 75)])

##layer2
layer2 <- c()
layer2 <- c(layer2,dat_3$Non_Area[which(dat_3$distance > 75 & dat_3$distance <= 150)],
            dat2_3$Non_Area[which(dat2_3$distance > 75 & dat2_3$distance <= 150)],
            dat3_3$Non_Area[which(dat3_3$distance > 75 & dat3_3$distance <= 150)])

##layer3
layer3 <- c()
layer3 <- c(layer3,dat_3$Non_Area[which(dat_3$distance > 150)],dat2_3$Non_Area[which(dat2_3$distance > 150)],
            dat3_3$Non_Area[which(dat3_3$distance > 150)])

layer <- c(rep("layer0",length(layer0)),rep("layer1",length(layer1)),rep("layer2",length(layer2)),rep("layer3",length(layer3)))
datdat <- data.frame(matrix(rep(NA,length(layer)*2),nrow = length(layer),ncol=2))
colnames(datdat) <- c("Area","layer")
datdat$layer <- layer
datdat$Area <- c(layer0,layer1,layer2,layer3)

png("C:/Users/monom/Documents/fatbody_DAPI_6/n_area.png", width = 5000, height = 3000, res = 500)

ggplot(datdat, aes(x = layer, y = Area, colour = layer)) +
  geom_boxplot(fill = "white", color = "black", alpha = 0.1, outlier.shape = NA, width = 0.5) + # 箱ひげ図
  geom_jitter(width = 0.2, size = 2, alpha = 0.3) + # ジッター
  scale_color_manual(values = c("#FF71FF", "#3A78D4", "#87E7C1", "#EE9E21")) + # 色設定
  theme_minimal() +
  labs(
    x = "",
    y = "面積 log10 (Area) μm²"
  ) +
  theme_bw() +
  theme(text = element_text(size = 25))

dev.off()