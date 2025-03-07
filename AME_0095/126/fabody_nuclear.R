Rown <- 0
Rowv <- c()
Area <- c()
for(RI in 1:3){
  tab <- read.csv(paste0("C:/Users/monom/Documents/fatbody_DAPI_5/Spz/Results_",RI,".csv"))
  Rown <- Rown + nrow(tab)
  Rowv <- c(Rowv,nrow(tab))
  Area <- c(Area,tab$Area)
}


A <- sort(as.numeric(Area[1:Rowv[1]]),decreasing = T,index= F)
N <- mean(A[346:length(A)])
AN <- Area[1:Rowv[1]]/N
Area[1:Rowv[1]] <- AN

A <- sort(as.numeric(Area[as.numeric(Rowv[1]+1):as.numeric(Rowv[1]+Rowv[2])]),decreasing = T,index= F)
N <- mean(A[230:length(A)])
AN <- Area[as.numeric(Rowv[1]+1):as.numeric(Rowv[1]+Rowv[2])]/N
Area[as.numeric(Rowv[1]+1):as.numeric(Rowv[1]+Rowv[2])] <- AN

A <- sort(as.numeric(Area[as.numeric(Rowv[1]+Rowv[2]+1):as.numeric(Rowv[1]+Rowv[2]+Rowv[3])]),decreasing = T,index= F)
N <- mean(A[42:length(A)])
AN <- Area[as.numeric(Rowv[1]+Rowv[2]+1):as.numeric(Rowv[1]+Rowv[2]+Rowv[3])]/N
Area[as.numeric(Rowv[1]+Rowv[2]+1):as.numeric(Rowv[1]+Rowv[2]+Rowv[3])] <- AN


Rown_2 <- 0
Rowv_2 <- c()
Area_2 <- c()
for(RI in 1:3){
  tab <- read.csv(paste0("C:/Users/monom/Documents/fatbody_DAPI_2/DAPI/Results_",RI,".csv"))
  Rown_2 <- Rown_2 + nrow(tab)
  Rowv_2 <- c(Rowv_2,nrow(tab))
  Area_2 <- c(Area_2,tab$Area)
}

A <- sort(as.numeric(Area_2[1:Rowv_2[1]]),decreasing = T,index= F)
N <- mean(A[221:length(A)])
AN <- Area_2[1:Rowv_2[1]]/N
Area_2[1:Rowv_2[1]]<- AN

A <- sort(as.numeric(Area_2[as.numeric(Rowv_2[1]+1):as.numeric(Rowv_2[1]+Rowv_2[2])]),
          decreasing = T,index= F)
N <- mean(A[217:length(A)])
AN <- Area_2[as.numeric(Rowv_2[1]+1):as.numeric(Rowv_2[1]+Rowv_2[2])]/N
Area_2[as.numeric(Rowv_2[1]+1):as.numeric(Rowv_2[1]+Rowv_2[2])] <- AN

A <- sort(as.numeric(Area_2[as.numeric(Rowv_2[1]+Rowv_2[2]+1):
                              as.numeric(Rowv_2[1]+Rowv_2[2]+Rowv_2[3])]),
          decreasing = T,index= F)
N <- mean(A[180:length(A)])
AN <- Area_2[as.numeric(Rowv_2[1]+Rowv_2[2]+1):
               as.numeric(Rowv_2[1]+Rowv_2[2]+Rowv_2[3])]/N
Area_2[as.numeric(Rowv_2[1]+Rowv_2[2]+1):
         as.numeric(Rowv_2[1]+Rowv_2[2]+Rowv_2[3])] <- AN

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

png("C:/Users/monom/Documents/fatbody_DAPI_5/Spz/area_individual_n.png", 
    width = 3000, height = 3000, res = 500)

ggplot(dat, aes(x = individual, y = Area, colour = Spz)) +
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