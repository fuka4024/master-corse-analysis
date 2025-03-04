##install.packages
install.packages("exactRankTests")
install.packages("ggplot2")
install.packages("ggsignif")

##load.packages
library(exactRankTests)
library(ggplot2)
library(ggsignif)

##set directory
setwd("C:/Users/monom/Documents/Spatzle_area_5/")

##pPIG-A3G-Spätzle3画像解析結果csvデータの格納(tidy data)
plotsize <- c()
Anum <- c()
med <- c()
v <- c()
q <- c()
rownum <- 0
for(RI in 1:10){
  area <- read.csv(paste0("C:/Users/monom/Documents/Spatzle_area_5/result/Results_",RI+1,".csv"))
  med <- c(med,median(log10(area$Area)))
  q <- c(q,quantile(log10(area$Area),0.75)/nrow(area))
  v <- c(v,median(log10(area$Area))/nrow(area))
  plotsize <- c(plotsize,area$Area)
  Anum <- c(Anum,nrow(area))
  rownum <- rownum+nrow(area)
}
tidy_frame <- as.data.frame(matrix(rep(NA,rownum*3),nrow=rownum,ncol=3))
colnames(tidy_frame) <- c("Area","Number","Spätzle3")
tidy_frame$Area <- plotsize

##individual Number
Kotai <- c()
for(RI in 1:10){
  Kotai <- c(Kotai,rep(RI,Anum[RI]))
}
tidy_frame$Number <- Kotai
tidy_frame$Spatzle <- rep("pPIG-A3G-Spätzle3",nrow(tidy_frame))

##pPIG-A3GR画像解析結果csvデータの格納(tidy data)
plotsize_2 <- c()
Anum_2 <- c()
v_2 <- c()
med_2 <- c()
q_2 <- c()
rownum_2 <- 0
for(RI in 11:17){
  area <- read.csv(paste0("C:/Users/monom/Documents/Spatzle_area_5/result/Results_",RI+1,".csv"))
  plotsize_2 <- c(plotsize_2,area$Area)
  q_2 <- c(q_2,quantile(log10(area$Area),0.75)/nrow(area))
  v_2 <- c(v_2,median(log10(area$Area))/nrow(area))
  med_2 <- c(med_2,median(log10(area$Area)))
  Anum_2 <- c(Anum_2,nrow(area))
  rownum_2 <- rownum_2+nrow(area)
}
tidy_frame_2 <- as.data.frame(matrix(rep(NA,rownum_2*3),nrow=rownum_2,ncol=3))
colnames(tidy_frame_2) <- c("Area","Number","Spätzle3")
tidy_frame_2$Area <- plotsize_2
Kotai_2 <- c()
for(RI in 1:7){
  Kotai_2 <- c(Kotai_2,rep(RI+10,Anum_2[RI]))
}
tidy_frame_2$Number <- Kotai_2
tidy_frame_2$Spatzle <- rep("pPIG-A3GR",nrow(tidy_frame_2))

tidy_frame_3 <- rbind(tidy_frame,tidy_frame_2)
tidy_frame_3$Area <- log10(tidy_frame_3$Area)
tidy_frame_4 <- tidy_frame_3[-which(tidy_frame_3$Area < 2.4),]


##random sampling & shapiro.test
normal_distribution_test_Spatzle <- list()
for(RI in 1:17){
  d1 <- tidy_frame_3$Area[which(tidy_frame_3$Number == RI)]
  random_sample <- sample(d1, size = 10, replace = FALSE) 
  normal_distribution_test_Spatzle[[RI]] <- shapiro.test(d1)
  print(RI)
}
random_sample <- sample(d1, size = 5, replace = FALSE) 

##box plot & jitter plot
png("C:/Users/monom/Documents/Spatzle_area_5/Spatzle_area.png",width=3000, height=3000, res=500)
  ggplot(tidy_frame_3,aes(x=Spatzle,y=Area))+
    geom_boxplot(fill = "white",color = "black", alpha = 0.5, outlier.shape = NA,width=0.5) + # 箱ひげ図（アウトライヤー非表示）
    geom_jitter(width = 0.2, size = 2,color = "gray40") + # 散布図を重ねる
    theme_minimal() +
    labs( 
         x = "Spatzle", 
         y = "面積 (Area)Area") +
    labs(x=element_blank(),y="log10(cell size) μm²")+
    theme_bw()+
    theme(text = element_text(size = 25))
  dev.off()
  
 
##individual plot
  tidy_frame_3$Number <- as.factor(tidy_frame_3$Number)
  png("C:/Users/monom/Documents/Spatzle_area_5/spz_1.png",width=3000, height=3000, res=500) 
  ggplot(tidy_frame_3, aes(x = Number, y = Area, color = Number)) +
    geom_violin(fill = "white",color = "black", alpha = 0.5, outlier.shape = NA,width=0.5) + # 箱ひげ図（アウトライヤー非表示）
    geom_jitter(width = 0.2, size = 2) +
    scale_color_manual(values = scales::hue_pal()(17)) + # 10色
    facet_wrap(~ Spatzle) + # Spatzleグループごとに分ける
    labs(
      title = "Spatzle グループ別の個体別 Area",
      x = "個体番号 (Number)",
      y = "面積 (Area)",
      color = "個体番号"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5), # タイトルを中央揃え
      strip.text = element_text(size = 14)  # ファセットラベルのサイズ調整
    )+
  theme_bw()
  dev.off()
  
  
  ##quantile histogram
  png("C:/Users/monom/Documents/Spatzle_area_5/q.png",width=3000, height=3000, res=500)
  barplot(c(q,q_2))
  dev.off()