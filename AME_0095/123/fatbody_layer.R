install.packages("nortest")  
install.packages("MASS")
library(nortest)
library(MASS)
library(lmerTest)
library(exactRankTests)

##all spz size
Rown <- 0
Rowv <- c()
Area <- c()
for(RI in 1:3){
  tab <- read.csv(paste0("C:/Users/monom/Documents/fatbody_DAPI_5/Spz/Results_",RI,".csv"))
  Rown <- Rown + nrow(tab)
  Rowv <- c(Rowv,nrow(tab))
  Area <- c(Area,tab$Area)
}
Area <- log10(Area)


png(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/Spz/Area_spz_DAPI_hist.png")
    ,width=2000, height=2000, res=300 )
hist(Area,breaks = 150)
dev.off()

###statistical test###
#指数分布#
ks.test(Area, "pexp", rate = 1/mean(Area))

#ポアソン分布#
area <- trunc(Area)
observed <- table(area)
lambda <- mean(area)
expected <- dpois(as.numeric(names(observed)), lambda) * sum(observed)


chisq.test(observed, p = expected/sum(expected))

##一様分布##
ks.test(Area, "punif", min = min(Area), max = max(Area))

##ガンマ分布##

Area <- Area[Area > 0]  
Area <- na.omit(Area)   

mean_area <- mean(Area)
var_area <- var(Area)
fit <- fitdistr(Area, "gamma", start = list(shape = mean_area^2 / var_area, rate = mean_area / var_area))
ks.test(Area, "pgamma", shape = fit$estimate[1], rate = fit$estimate[2])



##spz expression size
Rown_3 <- 0
Rowv_3 <- c()
Area_3 <- c()
for(RI in 4:6){
  tab <- read.csv(paste0("C:/Users/monom/Documents/fatbody_DAPI_5/Spz/Results_",RI,".csv"))
  Rown_3 <- Rown_3 + nrow(tab)
  Rowv_3 <- c(Rowv_3,nrow(tab))
  Area_3 <- c(Area_3,tab$Area)
}
Area_3 <- log10(Area_3)

result <- c()
for(RI in 1:1000){
  test <- sample(Area_3,20)
  p <- ks.test(test, "pnorm", mean = mean(test), sd = sd(test))
  result <- c(result,p$p.value)
}

png(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/Spz/Area_p_hist_2.png")
    ,width=2000, height=2000, res=300 )
hist(result,breaks = 150)
dev.off()

png(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/Spz/Area_spz_DAPI_hist_2.png")
    ,width=2000, height=2000, res=300 )
hist(Area_3,breaks = 150)
dev.off()

##layer1 spz expression size
Rown_4 <- 0
Rowv_4 <- c()
Area_4 <- c()
for(RI in 7:9){
  tab <- read.csv(paste0("C:/Users/monom/Documents/fatbody_DAPI_5/Spz/Results_",RI,".csv"))
  Rown_4 <- Rown_4 + nrow(tab)
  Rowv_4 <- c(Rowv_4,nrow(tab))
  Area_4 <- c(Area_4,tab$Area)
}
Area_4 <- log10(Area_4)


png(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/Spz/Area_spz_DAPI_hist_3.png")
    ,width=2000, height=2000, res=300 )
hist(Area_4,breaks = 150)
dev.off()

##layer2 spz expression size
Rown_5 <- 0
Rowv_5 <- c()
Area_5 <- c()
for(RI in 10:12){
  tab <- read.csv(paste0("C:/Users/monom/Documents/fatbody_DAPI_5/Spz/Results_",RI,".csv"))
  Rown_5 <- Rown_5 + nrow(tab)
  Rowv_5 <- c(Rowv_5,nrow(tab))
  Area_5 <- c(Area_5,tab$Area)
}
Area_5 <- log10(Area_5)

##layer 3 spz 
Rown_6 <- 0
Rowv_6 <- c()
Area_6 <- c()
for(RI in 13:15){
  tab <- read.csv(paste0("C:/Users/monom/Documents/fatbody_DAPI_5/Spz/Results_",RI,".csv"))
  Rown_6 <- Rown_6 + nrow(tab)
  Rowv_6 <- c(Rowv_6,nrow(tab))
  Area_6 <- c(Area_6,tab$Area)
}
Area_6 <- log10(Area_6)

##layer 4 spz 
Rown_7 <- 0
Rowv_7 <- c()
Area_7 <- c()
for(RI in 16:17){
  tab <- read.csv(paste0("C:/Users/monom/Documents/fatbody_DAPI_5/Spz/Results_",RI,".csv"))
  Rown_7 <- Rown_7 + nrow(tab)
  Rowv_7 <- c(Rowv_7,nrow(tab))
  Area_7 <- c(Area_7,tab$Area)
}
Area_7 <- log10(Area_7)

##control DAPI size
Rown_2 <- 0
Rowv_2 <- c()
Area_2 <- c()
for(RI in 1:3){
  tab <- read.csv(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/DAPI/Results_",RI,".csv"))
  Rown_2 <- Rown_2 + nrow(tab)
  Rowv_2 <- c(Rowv_2,nrow(tab))
  Area_2 <- c(Area_2,tab$Area)
}
Area_2 <- log10(Area_2)


png(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/DAPI/Area_DAPI_hist.png")
    ,width=2000, height=2000, res=300 )
hist(Area_2,breaks = 150)
dev.off()


##spz all vs control
Rown_X <- Rown + Rown_2
dat <- data.frame(matrix(rep(NA,Rown_X*3),nrow = Rown_X,ncol=3))
colnames(dat) <- c("Area","Spz","individual")
dat$Area <- c(Area,Area_2)


dat$Spz <- relevel(as.factor(c(rep("Spz+",Rown),rep("Spz-",Rown_2))),ref = "Spz-")

individual <- c()
for(RI in 1:length(Rowv)){
  individual <- c(individual,rep(RI,Rowv[RI]))
}

individual_2 <- c()
for(RI in 1:length(Rowv_2)){
  individual_2 <- c(individual_2,rep(RI+10,Rowv_2[RI]))
}

dat$individual <- as.factor(c(individual,individual_2))

png("C:/Users/monom/Documents/fatbody_DAPI_2/Spz/n_area.png",width=3000, height=3000, res=500)
ggplot(dat,aes(x=Spz,y=Area))+
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

spz  <- summary(m1 <- lmer(Area ~  Spz + (1|individual), data = dat))

png(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/Spz/resid_hist.png")
    ,width=2000, height=2000, res=300 )
hist(resid(m1),breaks = 150)
dev.off()

result <- c()
for(RI in 1:1000){
  test <- sample(resid(m1),50)
  p <- ks.test(test, "pnorm", mean = mean(test), sd = sd(test))
  result <- c(result,p$p.value)
}

png(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/Spz/resid_p_hist_1.png")
    ,width=2000, height=2000, res=300 )
hist(result,breaks = 150)
dev.off()
##ウィルコクソン##
wilcox.test(Area_4,Area_2)


##spz expression vs control
Rown_X <- Rown_3 + Rown_2
dat <- data.frame(matrix(rep(NA,Rown_X*3),nrow = Rown_X,ncol=3))
colnames(dat) <- c("Area","Spz","individual")
dat$Area <- c(Area_3,Area_2)


dat$Spz <- relevel(as.factor(c(rep("Spz+",Rown_3),rep("Spz-",Rown_2))),ref = "Spz-")

individual <- c()
for(RI in 1:length(Rowv_3)){
  individual <- c(individual,rep(RI,Rowv_3[RI]))
}

individual_2 <- c()
for(RI in 1:length(Rowv_2)){
  individual_2 <- c(individual_2,rep(RI+10,Rowv_2[RI]))
}

dat$individual <- as.factor(c(individual,individual_2))


spz  <- summary(m1 <- lmer(Area ~  Spz + (1|individual), data = dat))
png(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/Spz/resid_hist_2.png")
    ,width=2000, height=2000, res=300 )
hist(resid(m1),breaks = 150)
dev.off()



##around spz expression vs control
Rown_X <- Rown_4 + Rown_2
dat <- data.frame(matrix(rep(NA,Rown_X*3),nrow = Rown_X,ncol=3))
colnames(dat) <- c("Area","Spz","individual")
dat$Area <- c(Area_4,Area_2)


dat$Spz <- relevel(as.factor(c(rep("Spz+",Rown_4),rep("Spz-",Rown_2))),ref = "Spz-")

individual <- c()
for(RI in 1:length(Rowv_4)){
  individual <- c(individual,rep(RI,Rowv_4[RI]))
}

individual_2 <- c()
for(RI in 1:length(Rowv_2)){
  individual_2 <- c(individual_2,rep(RI+10,Rowv_2[RI]))
}

dat$individual <- as.factor(c(individual,individual_2))


spz  <- summary(m1 <- lmer(Area ~  Spz + (1|individual), data = dat))
png(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/Spz/resid_hist_3.png")
    ,width=2000, height=2000, res=300 )
hist(resid(m1),breaks = 150)
dev.off()














result <- c()
for(RI in 1:1000){
  test <- sample(Area,20)
  p <- ks.test(test, "pnorm", mean = mean(test), sd = sd(test))
  result <- c(result,p$p.value)
}

png(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/Spz/Area_p_hist.png")
    ,width=2000, height=2000, res=300 )
hist(result,breaks = 150)
dev.off()


write.csv(spz$coefficients,"C:/Users/monom/Documents/Spatzle_area_5/spz.csv")

var_test <- var.test(Area,Area_2)
var_test_2 <- var.test(Area_3,Area_2)

t_test <- t.test(Area, Area_2, var.equal = TRUE)


Rown_X <- Rown + Rown_2
dat <- data.frame(matrix(rep(NA,Rown_X*3),nrow = Rown_X,ncol=3))
colnames(dat) <- c("Area","Spz","individual")
dat$Area <- c(Area,Area_2)


dat$Spz <- as.factor(c(rep("Spz+",Rown),rep("Spz-",Rown_2)))

individual <- c()
for(RI in 1:length(Rowv)){
  individual <- c(individual,rep(RI,Rowv[RI]))
}

individual_2 <- c()
for(RI in 1:length(Rowv_2)){
  individual_2 <- c(individual_2,rep(RI+3,Rowv_2[RI]))
}


dat$individual <- as.factor(c(individual,individual_2))
dat_2 <- dat[-which(dat$Area<0),]

png("C:/Users/monom/Documents/fatbody_DAPI_5/Spz/area_individual.png", 
    width = 3000, height = 3000, res = 500)

ggplot(dat_2, aes(x = individual, y = Area, colour = Spz)) +
  geom_boxplot(fill = "white", color = "black",alpha = 0.5, outlier.shape = NA, width = 0.4) +  # ボックス幅を少し小さく調整
  geom_jitter(width = 0.2, size = 2, alpha = 0.2) +
  scale_color_manual(values = c("#B2ABD3","#FEB863")) + 
  theme_minimal() +
  labs(x = "", 
       y = "面積 log10 (Area) μm²") +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +  # x軸の余白を広げる
  theme_bw() +
  theme(text = element_text(size = 25))

dev.off()


Rown_X <- Rown + Rown_2 + Rown_3 + Rown_4 + Rown_5 + Rown_6 + Rown_7
dat <- data.frame(matrix(rep(NA,Rown_X*3),nrow = Rown_X,ncol=3))
colnames(dat) <- c("Area","Spz")
dat$Area <- c(Area,Area_2,Area_3,Area_4,Area_5,Area_6,Area_7)


dat$Spz <- as.factor(c(rep("Spz+",Rown),rep("Spz-",Rown_2),rep("layer0",Rown_3),
                       rep("layer1",Rown_4),rep("layer2",Rown_5),rep("layer3",Rown_6),rep("layer4",Rown_7)))
dat$Spz <- factor(dat$Spz, levels = c("layer0", "layer1", "layer2", "layer3", "layer4", "Spz+", "Spz-"))

# 最小値以下のデータを除去
dat_2 <- dat[-which(dat$Area < min(Area_2)),]

# プロット作成
png("C:/Users/monom/Documents/fatbody_DAPI_5/Spz/n_area.png", width = 5000, height = 3000, res = 500)

ggplot(dat_2, aes(x = Spz, y = Area, colour = Spz)) +
  geom_boxplot(fill = "white", color = "black", alpha = 0.1, outlier.shape = NA, width = 0.5) + # 箱ひげ図
  geom_jitter(width = 0.2, size = 2, alpha = 0.3) + # ジッター
  scale_color_manual(values = c("#69FF69", "#72FFFF", "#6C6CFF", "#FFFF6F", "#FF71FF", "#FEB863", "#B2ABD3")) + # 色設定
  theme_minimal() +
  labs(
    x = "",
    y = "面積 log10 (Area) μm²"
  ) +
  theme_bw() +
  theme(text = element_text(size = 25))

dev.off()


Rown_X <- Rown_2 + Rown_3 + Rown_4 + Rown_5 + Rown_6 + Rown_7
dat <- data.frame(matrix(rep(NA,Rown_X*3),nrow = Rown_X,ncol=3))
colnames(dat) <- c("Area","Spz")
dat$Area <- c(Area_2,Area_3,Area_4,Area_5,Area_6,Area_7)


dat$Spz <- as.factor(c(rep("コントロール",Rown_2),rep("layer0",Rown_3),
                       rep("layer1",Rown_4),rep("layer2",Rown_5),rep("layer3",Rown_6),rep("layer4",Rown_7)))
dat$Spz <- factor(dat$Spz, levels = c("layer0", "layer1", "layer2", "layer3", "layer4",  "コントロール"))

# 最小値以下のデータを除去
dat_2 <- dat[-which(dat$Area < min(Area_2)),]

# プロット作成
png("C:/Users/monom/Documents/fatbody_DAPI_5/Spz/n_area_2.png", width = 5000, height = 3000, res = 500)

ggplot(dat_2, aes(x = Spz, y = Area, colour = Spz)) +
  geom_boxplot(fill = "white", color = "black", alpha = 0.1, outlier.shape = NA, width = 0.5) + # 箱ひげ図
  geom_jitter(width = 0.2, size = 2, alpha = 0.3) + # ジッター
  scale_color_manual(values = c("#69FF69", "#72FFFF", "#6C6CFF", "#FFFF6F", "#FF71FF", "#B2ABD3")) + # 色設定
  theme_minimal() +
  labs(
    x = "",
    y = "面積 log10 (Area) μm²"
  ) +
  theme_bw() +
  theme(text = element_text(size = 25))

dev.off()