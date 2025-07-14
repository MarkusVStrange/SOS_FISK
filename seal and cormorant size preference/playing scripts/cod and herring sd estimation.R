library(icesDatras)
library(icesVocab)
library(ggplot2)
library(dplyr)
library(minpack.lm)
library(lubridate)
library(tidyverse)

setwd("/Users/mavast/Desktop/Markus SOS/seal and cormorant size preference/")
# Read data from DATRAS
#####
#hLdata <- getDATRAS(record="HL",survey = "BITS", years = 1970:2024, quarters = 1:4)
#cadata <- getDATRAS(record="CA",survey = "BITS", years = 1970:2024, quarters = 1:4)
#hhdata <- getDATRAS(record="HH",survey = "BITS", years = 1970:2024, quarters = 1:4)
#write.csv(hLdata, "hLdata.csv",row.names = FALSE)
#write.csv(cadata, "cadata.csv",row.names = FALSE)
#write.csv(hhdata, "hhdata.csv",row.names = FALSE)


hldata <- read.table("data/hldata.csv",sep=',',header=TRUE)
cadata <- read.table("data/cadata.csv",sep=',',header=TRUE)
hhdata <- read.table("data/hhdata.csv",sep=',',header=TRUE)


#####
rm(list=setdiff(ls(),c('hldata','cadata','hhdata')))


#Prepare cod data
#####
#aphia.cod <- icesVocab::findAphia(c("cod"))
aphia.cod <- 126436

hh_unique <- hhdata %>%
  dplyr::select(Country,Quarter, Ship, StNo, HaulNo, Year,StatRec) 

cadata <- left_join(cadata,hh_unique)

WB_rectangles <- c(
  # Subdivision 22
  "37G1", "37G2", "37G3", "37G4", "38G1", "38G2", "38G3", "38G4", "39G1", "39G2", "39G3", "39G4",
  # Subdivision 23
  "40G1", "40G2", "40G3", "41G1", "41G2", "41G3",
  # Subdivision 24
  "38H1", "38H2", "38H3", "39H1", "39H2", "39H3", "40H1", "40H2", "40H3"
)
specs <- data.frame(Valid_Aphia = c(aphia.cod),
                    species = c('cod'))
cod.ca  <- cadata %>% filter(Valid_Aphia %in% specs$Valid_Aphia  & AreaCode %in% WB_rectangles)
cod.ca <- cod.ca %>% left_join(specs, by='Valid_Aphia')

hhWB <- hhdata %>% filter(StatRec %in% WB_rectangles)
hhWB$jday <- yday(paste(hhWB$Year,'-',hhWB$Month,'-',hhWB$Day,sep=""))
survey.day <- aggregate(jday~Year+Quarter,data=hhWB,FUN=mean)
survey.day$jday <- survey.day$jday/365
quarter.day <- aggregate(jday~Quarter,data=survey.day,FUN=mean)

#LngtCode: . = mm; 0 = 0.5 cm; 1 = cm
#Convert the length unit to cm
table(cod.ca$LngtCode)
cod.ca$LngtClass[cod.ca$LngtCode=="."] <- cod.ca$LngtClass[cod.ca$LngtCode=="."]/10

cod.ca <- cod.ca %>% left_join(quarter.day,by=c('Quarter'))
cod.ca <- cod.ca %>% 
  mutate(t_hatch = 61/365)
cod.ca$fish.age <- cod.ca$Age+cod.ca$jday-cod.ca$t_hatch
cod.ca <- cod.ca %>% filter(Age>1)

cod.ca <- cod.ca %>% filter(!(LngtClass<20 & Age>5)) # remove a 10 cm 8 year old cod. Must be a mistake
cod.age <- aggregate(LngtClass~Age+Year+Quarter+species+fish.age+jday+t_hatch,
                          data=cod.ca,FUN=mean)
#cod.sd <- aggregate(LngtClass~Age+Year+Quarter+species+fish.age+jday,
 #                    data=cod.ca,FUN=sd)
#####
rm(list=setdiff(ls(),c('hLdata','cadata','hhdata','cod.age','cod.sd','cod.ca')))

# Estimate cod variability
#####
# First find the sd pr. age for each sampling
cod.ca$qy <- paste(cod.ca$Quarter,cod.ca$Year)
qy.s <- unique(cod.ca$qy)
sd.s <- data.frame(fish.age = 0,
                   Year = 0,qy = 0, LngtClass = 0)
for (i in 1:length(qy.s)){
  dada <- cod.ca %>% filter(qy == qy.s[i])
  std <- aggregate(LngtClass~fish.age+Year+qy,data = dada,FUN=sd)
  sd.s <- rbind(sd.s,std)
}
sd.s <- sd.s %>% filter(!(is.na(LngtClass)) & LngtClass>0)
sd.mean <- aggregate(LngtClass~fish.age,data = sd.s,FUN=mean)
sd.mean$n <- as.numeric(table(sd.s$fish.age))

# Function to calculate sd of L as a function of t. Based on partial derivatives of L - 
calculate_sd <- function(t, k_sd) {
  L_hatch <- 0.4 # length at hatching, cm (Pepin et al, 1997)  t_hatch <- 32/365 # peak spawning mid march 
  fish.age <- t-t_hatch
  u <- 0.22 # Amplitude of seasonality
  w <- 0.90 # Timing of peak growth
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  k <- 0.11
  Linf <- 154.56
  
  # Partial derivatives
  dL_dLinf <- -exp(-k * (phi_t - phi_hatch +fish.age)) + 1
  dL_dk <- (L_hatch - Linf) * (phi_hatch - phi_t -fish.age ) * exp(-k * (phi_t +fish.age - phi_hatch))
  dL_dage <- -(L_hatch - Linf)*k*exp(-k*(phi_t + fish.age - phi_hatch))

  age_sd <- 0.1 # based on spawning duration
  Linf_sd <- 16.75
  # Variance
  var_L <- (dL_dLinf * Linf_sd)^2 + (dL_dk * k_sd)^2 + (dL_dage * age_sd)^2
  
  # Standard deviation
  sd_L <- sqrt(var_L)
  
  return(sd_L)
} #Original

# Fit parameter varibility to observed variability
y.c <- sd.mean$LngtClass
x.c <- sd.mean$fish.age
n.c <- sd.mean$n
cod.sd <- nlsLM(y.c~calculate_sd(x.c,k_sd),
                     start = list(k_sd=0.01),weights = n.c)
coef(cod.sd)
par(mfrow=c(1,1))
t.c <- ((0:2500)/100)[((0:2500)/100)>cod.age$t_hatch[1]]
plot(t.c, calculate_sd(t.c,k_sd=coef(cod.sd)),
     type = "l", col = "black", lwd = 2,
     ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,15),ylim=c(0,20))
points(sd.mean$fish.age,sd.mean$LngtClass,xlim=c(0,15),ylim=c(0,20))
#lines(t.c,vbgr.sd(t.c,1000),lwd=2,col='red')
#####

#Prepare herring data
#####
rm(list=setdiff(ls(),c('hLdata','cadata','hhdata')))
#aphia.cod <- icesVocab::findAphia(c("cod"))
aphia.herring <- 126417

hh_unique <- hhdata %>%
  dplyr::select(Country,Quarter, Ship, StNo, HaulNo, Year,StatRec) 

cadata <- left_join(cadata,hh_unique)

WB_rectangles <- c(
  # Subdivision 22
  "37G1", "37G2", "37G3", "37G4", "38G1", "38G2", "38G3", "38G4", "39G1", "39G2", "39G3", "39G4",
  # Subdivision 23
  "40G1", "40G2", "40G3", "41G1", "41G2", "41G3",
  # Subdivision 24
  "38H1", "38H2", "38H3", "39H1", "39H2", "39H3", "40H1", "40H2", "40H3"
)
specs <- data.frame(Valid_Aphia = c(aphia.herring),
                    species = c('herring'))
herring.ca  <- cadata %>% filter(Valid_Aphia %in% specs$Valid_Aphia  & AreaCode %in% WB_rectangles)
herring.ca <- herring.ca %>% left_join(specs, by='Valid_Aphia')

hhWB <- hhdata %>% filter(StatRec %in% WB_rectangles)
hhWB$jday <- yday(paste(hhWB$Year,'-',hhWB$Month,'-',hhWB$Day,sep=""))
survey.day <- aggregate(jday~Year+Quarter,data=hhWB,FUN=mean)
survey.day$jday <- survey.day$jday/365
quarter.day <- aggregate(jday~Quarter,data=survey.day,FUN=mean)

#LngtCode: . = mm; 0 = mm cm; 1 = cm
#Convert the length unit to cm
table(herring.ca$LngtCode)
herring.ca$LngtClass <- round(herring.ca$LngtClass/10)


plot(1:length(herring.ca$LngtClass),herring.ca$LngtClass)

herring.ca <- herring.ca %>% left_join(quarter.day,by=c('Quarter'))
herring.ca <- herring.ca %>% 
  mutate(t_hatch = 0.3078767)
herring.ca$fish.age <- herring.ca$Age+herring.ca$jday-herring.ca$t_hatch
plot(herring.ca$fish.age,herring.ca$LngtClass)
herring.ca <- herring.ca %>% filter(fish.age>1)

herring.age <- aggregate(LngtClass~Age+Year+Quarter+species+fish.age+jday+t_hatch,
                     data=herring.ca,FUN=mean)
#herring.sd <- aggregate(LngtClass~Age+Year+Quarter+species+fish.age+jday,
#                    data=herring.ca,FUN=sd)
#####
rm(list=setdiff(ls(),c('hldata','cadata','hhdata','herring.age','herring.sd','herring.ca')))


# Estimate herring variability
#####
# First find the sd pr. age for each sampling
herring.ca$qy <- paste(herring.ca$Quarter,herring.ca$Year)
qy.s <- unique(herring.ca$qy)
sd.s <- data.frame(fish.age = 0,
                   Year = 0,qy = 0, LngtClass = 0)
for (i in 1:length(qy.s)){
  dada <- herring.ca %>% filter(qy == qy.s[i])
  std <- aggregate(LngtClass~fish.age+Year+qy,data = dada,FUN=sd)
  sd.s <- rbind(sd.s,std)
}
sd.s <- sd.s %>% filter(!(is.na(LngtClass)) & LngtClass>0)
sd.mean <- aggregate(LngtClass~fish.age,data = sd.s,FUN=mean)
sd.mean$n <- as.numeric(table(sd.s$fish.age))

plot(sd.mean$fish.age,sd.mean$LngtClass)
# Function to calculate sd of L as a function of t. Based on partial derivatives of L - 
calculate_sd <- function(t,Linf_sd,k_sd) {
  #Parameters 2005-2010 for Western Baltic spring spawning herring (Gröhsler et al., 2013)
  
  L_hatch <- 0.645 # figure 2 in Bauer et al., 2014
  t_hatch <- mean(c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365) # peak spawning from figure 4 in Polte et al., 2021
  fish.age <- t-t_hatch
  u <- 0.22 # Amplitude of seasonality
  w <- 0.90 # Timing of peak growth
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  k <- 0.4705 #Parameters 2005-2010 for Western Baltic spring spawning herring (Gröhsler et al., 2013)
  Linf <- 30.57 #Parameters 2005-2010 for Western Baltic spring spawning herring (Gröhsler et al., 2013)
  
  # Partial derivatives
  dL_dLinf <- -exp(-k * (phi_t - phi_hatch +fish.age)) + 1
  dL_dk <- (L_hatch - Linf) * (phi_hatch - phi_t -fish.age ) * exp(-k * (phi_t +fish.age - phi_hatch))
  dL_dage <- -(L_hatch - Linf)*k*exp(-k*(phi_t + fish.age - phi_hatch))
  
  age_sd <- 0.05
 # Linf_sd <- 16.75
  # Variance
  var_L <- (dL_dLinf * Linf_sd)^2 + (dL_dk * k_sd)^2 + (dL_dage * age_sd)^2
  
  # Standard deviation
  sd_L <- sqrt(var_L)
  
  return(sd_L)
} #Original

# Fit parameter varibility to observed variability
y.h <- sd.mean$LngtClass
x.h <- sd.mean$fish.age
n.h <- sd.mean$n
herring.sd <- nlsLM(y.h~calculate_sd(x.h,Linf_sd,k_sd),
                start = list(Linf_sd=2,k_sd=0.09),weights = n.h)
coef(herring.sd)
par(mfrow=c(1,1))
t.h <- ((0:2500)/100)[((0:2500)/100)>herring.age$t_hatch[1]]
plot(t.h, calculate_sd(t.h,Linf_sd=coef(herring.sd)[1],k_sd=coef(herring.sd)[2]),
     type = "l", col = "black", lwd = 2,
     ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,15),ylim=c(0,4))
points(sd.mean$fish.age,sd.mean$LngtClass)
#lines(t.h,vbgr.sd(t.h,10000),lwd=2,col='red')
#####














# Simulation-based ground truthing
#####
t <- (0:250)/10
X_t <- t-floor(t) # time of year [0,1]
u <- 0.22 # Amplitude of seasonality
w <- 0.90 # Timing of peak growth
phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth


# Cod
t_hatch <- 32/365 # peak spawning mid march 
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
L_hatch <- 0.4
sd.0 <- matrix(rep(0,length(t)*100000),nrow = 100000)
for (i in 1:100000){
  k <- rnorm(1,mean=0.11,sd=coef(cod.sd))
  Linf <- rnorm(1,mean=154.56,sd=16.75)
  fish.age <- t-t_hatch+rnorm(1,0,0.1)
  L <- (L_hatch-Linf)*exp(-k*(phi_t+fish.age-phi_hatch))+Linf
  #lines(t,L,col='red')
  sd.0[i,] <- L
}
plot(t,apply(sd.0,2,sd),xlim=c(phi_hatch,15),ylim=c(0,20),lwd=2)
lines(c(phi_hatch,phi_hatch),c(-10,20))
lines(t,calculate_sd(t,k_sd=coef(cod.sd)),lwd=3,col='orange')

# herring
t_hatch <- 0.3078767 # peak spawning mid march 
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
L_hatch <- 0.645
sd.0 <- matrix(rep(0,length(t)*100000),nrow = 100000)
for (i in 1:100000){
  k <- rnorm(1,mean=0.453,sd=coef(herring.sd)[2])
  Linf <- rnorm(1,mean=30.57,sd=coef(herring.sd)[1])
  fish.age <- t-t_hatch+rnorm(1,0,0.05)
  L <- (L_hatch-Linf)*exp(-k*(phi_t+fish.age-phi_hatch))+Linf
  #lines(t,L,col='red')
  sd.0[i,] <- L
}
plot(t,apply(sd.0,2,sd),xlim=c(phi_hatch,15),ylim=c(0,5),lwd=2)
lines(c(phi_hatch,phi_hatch),c(-10,20))
lines(t,calculate_sd(t,Linf_sd=coef(herring.sd)[1],k_sd=coef(herring.sd)[2]),lwd=3,col='orange')
#####