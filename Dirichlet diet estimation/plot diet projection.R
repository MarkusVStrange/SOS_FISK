
# plot predicted cormorant diet - quarter normal scale food
#####
cormFood <- read.table("data/food_sim.csv",header=TRUE,sep=';')
prey <- c("cod","flatfish")
cormFood <- cormFood %>% filter(species!="herring")
cormFood$species[!(cormFood$species %in% prey)] <- "flatfish"
cormFood <- aggregate(biomass~sampling+species,data=cormFood,FUN = sum)


cf <- aggregate(biomass~sampling+species,data=cormFood,FUN = sum)

meanC.B <- aggregate(biomass~species,data=cf,FUN=mean)$biomass
C.idx <- match(cf$species,prey)
cf$Bnorm <- (cf$biomass)/meanC.B[C.idx]

k <- c(0.2290891,0.1069550)
p <- matrix(0,ncol=3,nrow=33)
for(i in 1:33){
  food.i <- c(cf$Bnorm[i],cf$Bnorm[i+33])
  other.prop <- f.other(sum(food.i*weights))
  p[i,] <- c(((food.i+1)^k)/(sum((food.i+1)^k)/(1-other.prop)),other.prop)
  
}

df1 <- as.data.frame(p)
names(df1) <- c("cod","flatfish","other")
df1$year <- 1991:2023

df <- df1 %>%
  mutate(
    cod_cum = cod,
    flatfish_cum = cod + flatfish,
    other_cum = cod + flatfish +other
    
  )
library(ggplot2)

ggplot(df, aes(x = year)) +
  geom_ribbon(aes(ymin = 0, ymax = cod_cum), fill = "blue", alpha = 0.5) +
  geom_ribbon(aes(ymin = cod_cum, ymax = flatfish_cum), fill = "green", alpha = 0.5) +
  geom_ribbon(aes(ymin = flatfish_cum, ymax = other_cum), fill = "red", alpha = 0.5) +
  ylab("Proportion") +
  xlab("Year") +
  theme_minimal()


#####

# plot predicted seal diet - quarter normal scale food
#####
cormFood <- read.table("data/food_sim.csv",header=TRUE,sep=';')
prey <- c("cod","flatfish")
cormFood <- cormFood %>% filter(species!="herring")
cormFood$species[!(cormFood$species %in% prey)] <- "flatfish"
cormFood <- aggregate(biomass~sampling+species,data=cormFood,FUN = sum)


cf <- aggregate(biomass~sampling+species,data=cormFood,FUN = sum)

meanC.B <- aggregate(biomass~species,data=cf,FUN=mean)$biomass
C.idx <- match(cf$species,prey)
cf$Bnorm <- (cf$biomass)/meanC.B[C.idx]
weights <- c(0.4756446, 0.5243554) # for 1. cod, 2. flatfish
cf$weights <- weights[C.idx]
cf$Bweight <- cf$Bnorm*cf$weights

fother <- aggregate(Bweight~sampling,data=cf,FUN=sum)
mean.other.diet <-  0.1143869
f.other <- function(x, m = mean.other.diet, R = 3) {
  1 / (1 + ((1/m - 1) * (x )^R))
}

k <- c(-3.612684,-3.646271)
p <- matrix(0,ncol=3,nrow=33)
for(i in 1:33){
  food.i <- c(cf$Bnorm[i],cf$Bnorm[i+33])
  other.prop <- f.other(sum(food.i*weights))
  p[i,] <- c(((food.i+1)^k)/(sum((food.i+1)^k)/(1-other.prop)),other.prop)
  
}

df1 <- as.data.frame(p)
names(df1) <- c("cod","flatfish","other")
df1$year <- 1991:2023

df <- df1 %>%
  mutate(
    cod_cum = cod,
    flatfish_cum = cod + flatfish,
    other_cum = cod + flatfish +other
    
  )
library(ggplot2)

ggplot(df, aes(x = year)) +
  geom_ribbon(aes(ymin = 0, ymax = cod_cum), fill = "blue", alpha = 0.5) +
  geom_ribbon(aes(ymin = cod_cum, ymax = flatfish_cum), fill = "green", alpha = 0.5) +
  geom_ribbon(aes(ymin = flatfish_cum, ymax = other_cum), fill = "red", alpha = 0.5) +
  ylab("Proportion") +
  xlab("X") +
  theme_minimal()


#####

# Plot food index
plot(1991:2023,log(cf$Bnorm[1:33]),type='l',lwd=2,col="indianred2",ylab="Normalized biomass index",xlab="Year")
lines(1991:2023,log(cf$Bnorm[34:66]),type='l',lwd=2,col="olivedrab",ylab="Normalized biomass index",xlab="Year")
legend("topright",               # legend position (or use x/y coordinates)
       legend = c("cod", "flatfish"),  # labels
       col = c("indianred2", "olivedrab"),         # matching colors
       lty = 1,                        # line type
       lwd = 2) 



# plot predicted cormorant diet - month normal scale food
#####
cormFood <- read.table("data/food_sim.csv",header=TRUE,sep=';')
prey <- c("cod","flatfish")
cormFood <- cormFood %>% filter(species!="herring")
cormFood$species[!(cormFood$species %in% prey)] <- "flatfish"
cormFood <- aggregate(biomass~sampling+species,data=cormFood,FUN = sum)


cf <- aggregate(biomass~sampling+species,data=cormFood,FUN = sum)

meanC.B <- aggregate(biomass~species,data=cf,FUN=mean)$biomass
C.idx <- match(cf$species,prey)
cf$Bnorm <- (cf$biomass)/meanC.B[C.idx]

k <- exp(c(-0.4770850,-0.1231021,-0.2052709))
p <- matrix(0,ncol=3,nrow=33)
for(i in 1:33){
  food.i <- c(c(cf$Bnorm[i],cf$Bnorm[i+33])*100,100)

  p[i,] <- (food.i^k)/sum(food.i^k)
  
}
rowSums(p)

df1 <- as.data.frame(p)
names(df1) <- c("cod","flatfish","other")
df1$year <- 1991:2023

df <- df1 %>%
  mutate(
    cod_cum = cod,
    flatfish_cum = cod + flatfish,
    other_cum = cod + flatfish +other
    
  )
library(ggplot2)

ggplot(df, aes(x = year)) +
  geom_ribbon(aes(ymin = 0, ymax = cod_cum), fill = "blue", alpha = 0.5) +
  geom_ribbon(aes(ymin = cod_cum, ymax = flatfish_cum), fill = "green", alpha = 0.5) +
  geom_ribbon(aes(ymin = flatfish_cum, ymax = other_cum), fill = "red", alpha = 0.5) +
  ylab("Proportion") +
  xlab("Year") +
  theme_minimal()


#####



# plot predicted cormorant diet - month log-scale food
#####
cormFood <- read.table("data/Cormfood_sim.csv",header=TRUE,sep=';')
prey <- c("cod","flatfish")
cormFood <- cormFood %>% filter(species!="herring")
cormFood$species[!(cormFood$species %in% prey)] <- "flatfish"
cormFood <- aggregate(biomass~sampling+species,data=cormFood,FUN = sum)


cf <- aggregate(biomass~sampling+species,data=cormFood,FUN = sum)

meanC.B <- aggregate(biomass~species,data=cf,FUN=mean)$biomass
C.idx <- match(cf$species,prey)
cf$Bnorm <- (cf$biomass)/meanC.B[C.idx]

k <- exp(c(1.144524,1.540565,1.313913))
p <- matrix(0,ncol=3,nrow=33)
for(i in 1:33){
  food.i <- log(c(c(cf$Bnorm[i],cf$Bnorm[i+33])*100,100))
  
  p[i,] <- (food.i^k)/sum(food.i^k)
  
}
rowSums(p)

df1 <- as.data.frame(p)
names(df1) <- c("cod","flatfish","other")
df1$year <- 1991:2023

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
  ylab("Diet proportion") +
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
plot(1991:2023,(cf$Bnorm[1:33]),type='l',lwd=2,col="indianred2",
     ylab="Normalized biomass index",xlab="Year",ylim=c(0,4))
lines(1991:2023,(cf$Bnorm[34:66]),type='l',lwd=2,col="olivedrab",ylab="Normalized biomass index",xlab="Year")
lines(c(1980,2030),c(1,1),col="lightskyblue",lwd=2,lty="dashed")
legend("top",               # legend position (or use x/y coordinates)
       legend = c("cod", "flatfish","other"),  # labels
       col = c("indianred2", "olivedrab","lightskyblue"),         # matching colors
       lty = c(1,1,2),                        # line type
       lwd = 2) 


# plot predicted grey seal diet - month log-scale food
#####
sealFood <- read.table("data/Sealfood_sim.csv",header=TRUE,sep=';')
prey <- c("cod","flatfish")
sealFood <- sealFood %>% filter(species!="herring")
sealFood$species[!(sealFood$species %in% prey)] <- "flatfish"
sealFood <- aggregate(biomass~sampling+species,data=sealFood,FUN = sum)


sf <- aggregate(biomass~sampling+species,data=sealFood,FUN = sum)

meanC.B <- aggregate(biomass~species,data=sf,FUN=mean)$biomass
C.idx <- match(sf$species,prey)
sf$Bnorm <- (sf$biomass)/meanC.B[C.idx]

k <- exp(c(0.3775862,0.2061222,-21.8424387))
p <- matrix(0,ncol=3,nrow=33)
for(i in 1:33){
  food.i <- log(c(c(sf$Bnorm[i],sf$Bnorm[i+33])*100,100))
  
  p[i,] <- (food.i^k)/sum(food.i^k)
  
}
rowSums(p)

df1 <- as.data.frame(p)
names(df1) <- c("cod","flatfish","other")
df1$year <- 1991:2023

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
  ylab("Diet proportion") +
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
plot(1991:2023,(sf$Bnorm[1:33]),type='l',lwd=2,col="indianred2",
     ylab="Normalized biomass index",xlab="Year",ylim=c(0,4))
lines(1991:2023,(sf$Bnorm[34:66]),type='l',lwd=2,col="olivedrab",ylab="Normalized biomass index",xlab="Year")
lines(c(1980,2030),c(1,1),col="lightskyblue",lwd=2,lty="dashed")
legend("top",               # legend position (or use x/y coordinates)
       legend = c("cod", "flatfish","other"),  # labels
       col = c("indianred2", "olivedrab","lightskyblue"),         # matching colors
       lty = c(1,1,2),                        # line type
       lwd = 2) 


n_year <- length(df_corm$year)
df_c <- data.frame(year=rep(rep(df_corm$year,each=4),3),
                   quarter = rep(c("Q1","Q2","Q3","Q4"),3*n_year),
                   prey =rep(c("cod","flatfish","other"),each=4*n_year),
                   diet=c(rep(df_corm$cod,each=4),rep(df_corm$flatfish,each=4),rep(df_corm$other,each=4)),
                   predator=rep("cormorant",3*4*n_year))

df_s <- data.frame(year=rep(rep(df_seal$year,each=4),3),
                   quarter = rep(c("Q1","Q2","Q3","Q4"),3*n_year),
                   prey =rep(c("cod","flatfish","other"),each=4*n_year),
                   diet=c(rep(df_seal$cod,each=4),rep(df_seal$flatfish,each=4),rep(df_seal$other,each=4)),
                   predator=rep("grey seal",3*4*n_year))


pred_diet <- rbind(df_c,df_s)
pred_diet <- aggregate(diet~prey+quarter+year+predator,data=pred_diet,FUN=mean)
write.table(pred_diet,"pred_diet.csv",row.names = FALSE,sep=';')


