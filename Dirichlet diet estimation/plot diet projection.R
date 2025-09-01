data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="") # data working directory

# plot predicted cormorant diet - month log-scale food
#####
cormFood <- read.table(paste(data_wd,"Cormorantfood_sim_quarter.csv",sep=""),header=TRUE,sep=';')
prey <- c("cod","flatfish")
cormFood <- cormFood %>% filter(species!="herring")
cormFood$species[!(cormFood$species %in% prey)] <- "flatfish"
cormFood <- aggregate(B_index~sampling+species,data=cormFood,FUN = sum)


cf <- aggregate(B_index~sampling+species,data=cormFood,FUN = sum)

meanC.B <- aggregate(B_index~species,data=cf,FUN=mean)$B_index
C.idx <- match(cf$species,prey)
cf$Bnorm <- (cf$B_index)/meanC.B[C.idx]

k <- exp(c(1.144524,1.540565,1.313913))
p <- matrix(0,ncol=3,nrow=33*4)
for(i in 1:(33*4)){
  food.i <- log(c(c(cf$Bnorm[i],cf$Bnorm[i+33*4])*100,100))
  
  p[i,] <- (food.i^k)/sum(food.i^k)
  
}
rowSums(p)

df1 <- as.data.frame(p)
names(df1) <- c("cod","flatfish","other")
df1$year <- rep(1991:2023,each=4)+(c(1/4,2/4,3/4,4/4)-1/8)

df_corm <- df1 %>%
  mutate(
    cod_cum = cod,
    flatfish_cum = cod + flatfish,
    other_cum = cod + flatfish +other
    
  )
library(ggplot2)

ggplot(df_corm, aes(x = year)) +
  geom_ribbon(aes(ymin = 0, ymax = cod_cum), fill = "indianred2", alpha = 0.5) +
  geom_ribbon(aes(ymin = cod_cum, ymax = flatfish_cum), fill = "olivedrab", alpha = 0.5) +
  geom_ribbon(aes(ymin = flatfish_cum, ymax = other_cum), fill = "lightskyblue", alpha = 0.5) +
  ylab("Diet proportion") +ggtitle("Cormorant food")+
  xlab("Year") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text=element_text(color="black", size=14,family="serif"),
        legend.position = "",
        legend.title = element_blank())


#####
# Plot cormorant food index
par(mfrow=c(1,1))
plot(rep(1991:2023,each=4)+(c(1/4,2/4,3/4,4/4)-1/8),cf$Bnorm[cf$species=="cod"],type='l',lwd=2,col="indianred2",
     ylab="Normalized biomass index",xlab="Year",ylim=c(0,4),main="Cormorant food")
lines(rep(1991:2023,each=4)+(c(1/4,2/4,3/4,4/4)-1/8),cf$Bnorm[cf$species=="flatfish"],type='l',lwd=2,col="olivedrab",ylab="Normalized biomass index",xlab="Year")
lines(c(1980,2030),c(1,1),col="lightskyblue",lwd=2,lty="dashed")
legend("top",               # legend position (or use x/y coordinates)
       legend = c("cod", "flatfish","other"),  # labels
       col = c("indianred2", "olivedrab","lightskyblue"),         # matching colors
       lty = c(1,1,2),                        # line type
       lwd = 2) 


# plot predicted grey seal diet - month log-scale food
#####
sealFood <- read.table(paste(data_wd,"Sealfood_sim_quarter.csv",sep=""),header=TRUE,sep=';')
prey <- c("cod","flatfish")
sealFood <- sealFood %>% filter(species!="herring")
sealFood$species[!(sealFood$species %in% prey)] <- "flatfish"
sealFood <- aggregate(B_index~sampling+species,data=sealFood,FUN = sum)

sf <- aggregate(B_index~sampling+species,data=sealFood,FUN = sum)

meanC.B <- aggregate(B_index~species,data=sf,FUN=mean)$B_index
C.idx <- match(sf$species,prey)
sf$Bnorm <- (sf$B_index)/meanC.B[C.idx]

k <- exp(c(0.3775862,0.2061222,-21.8424387))
p <- matrix(0,ncol=3,nrow=33*4)
for(i in 1:(33*4)){
  food.i <- log(c(c(sf$Bnorm[i],sf$Bnorm[i+33*4])*100,100))
  
  p[i,] <- (food.i^k)/sum(food.i^k)
  
}
rowSums(p)

df1 <- as.data.frame(p)
names(df1) <- c("cod","flatfish","other")
df1$year <- rep(1991:2023,each=4)+(c(1/4,2/4,3/4,4/4)-1/8)

df_seal <- df1 %>%
  mutate(
    cod_cum = cod,
    flatfish_cum = cod + flatfish,
    other_cum = cod + flatfish +other
    
  )
library(ggplot2)

ggplot(df_seal, aes(x = year)) +
  geom_ribbon(aes(ymin = 0, ymax = cod_cum), fill = "indianred2", alpha = 0.5) +
  geom_ribbon(aes(ymin = cod_cum, ymax = flatfish_cum), fill = "olivedrab", alpha = 0.5) +
  geom_ribbon(aes(ymin = flatfish_cum, ymax = other_cum), fill = "lightskyblue", alpha = 0.5) +
  ylab("Diet proportion") +ggtitle("Grey Seal food")+
  xlab("Year") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text=element_text(color="black", size=14,family="serif"),
        legend.position = "",
        legend.title = element_blank())


#####

# Plot greay seal food index
plot(rep(1991:2023,each=4)+(c(1/4,2/4,3/4,4/4)-1/8),sf$Bnorm[sf$species=="cod"],type='l',lwd=2,col="indianred2",
     ylab="Normalized biomass index",xlab="Year",ylim=c(0,4),main="Grey seal food")
lines(rep(1991:2023,each=4)+(c(1/4,2/4,3/4,4/4)-1/8),sf$Bnorm[sf$species=="flatfish"],type='l',lwd=2,col="olivedrab",ylab="Normalized biomass index",xlab="Year")
lines(c(1980,2030),c(1,1),col="lightskyblue",lwd=2,lty="dashed")
legend("top",               # legend position (or use x/y coordinates)
       legend = c("cod", "flatfish","other"),  # labels
       col = c("indianred2", "olivedrab","lightskyblue"),         # matching colors
       lty = c(1,1,2),                        # line type
       lwd = 2) 


n_year <- length(df_corm$year)
df_c <- data.frame(year=rep(floor(df_corm$year),3),
                   quarter = c("Q1","Q2","Q3","Q4"),
                   prey =rep(c("cod","flatfish","other"),each=n_year),
                   diet=c(df_corm$cod,df_corm$flatfish,df_corm$other),
                   predator=rep("cormorant",3*n_year))

df_s <- data.frame(year=rep(floor(df_seal$year),3),
                   quarter = c("Q1","Q2","Q3","Q4"),
                   prey =rep(c("cod","flatfish","other"),each=n_year),
                   diet=c(df_seal$cod,df_seal$flatfish,df_seal$other),
                   predator=rep("grey seal",3*n_year))


pred_diet <- rbind(df_c,df_s)
pred_diet <- aggregate(diet~prey+quarter+year+predator,data=pred_diet,FUN=mean)
#write.table(pred_diet,paste(data_wd,"pred_diet.csv",sep=""),row.names = FALSE,sep=';')


