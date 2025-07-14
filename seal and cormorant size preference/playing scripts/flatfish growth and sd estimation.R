library(icesDatras)
library(icesVocab)
library(ggplot2)
library(dplyr)
library(minpack.lm)
library(lubridate)
library(tidyverse)


setwd("/Users/mavast/Desktop/Markus SOS")
# Read data from DATRAS
#####
#hLdata <- getDATRAS(record="HL",survey = "BITS", years = 1970:2024, quarters = 1:4)
#cadata <- getDATRAS(record="CA",survey = "BITS", years = 1970:2024, quarters = 1:4)
#hhdata <- getDATRAS(record="HH",survey = "BITS", years = 1970:2024, quarters = 1:4)
#write.csv(hLdata, "hLdata.csv",row.names = FALSE)
#write.csv(cadata, "cadata.csv",row.names = FALSE)
#write.csv(hhdata, "hhdata.csv",row.names = FALSE)
#aphia.dab <- icesVocab::findAphia(c("dab"))
#aphia.plaice <- icesVocab::findAphia(c("plaice"))
#aphia.flounder <- icesVocab::findAphia(c("flounder"))

hldata <- read.table("hldata.csv",sep=',',header=TRUE)
cadata <- read.table("cadata.csv",sep=',',header=TRUE)
hhdata <- read.table("hhdata.csv",sep=',',header=TRUE)
aphia.dab <- 127139
aphia.plaice <- 127143
aphia.flounder <- 127141
#####
rm(list=setdiff(ls(),c('hldata','cadata','hhdata','aphia.dab','aphia.plaice','aphia.flounder')))

#Prepare data
#####
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
specs <- data.frame(Valid_Aphia = c(aphia.flounder,aphia.plaice,aphia.dab),
                    species = c('flounder','plaice','dab'))
flatfish.ca  <- cadata %>% filter(Valid_Aphia %in% specs$Valid_Aphia  & StatRec %in% WB_rectangles)
flatfish.ca <- flatfish.ca %>% left_join(specs, by='Valid_Aphia')

hhWB <- hhdata %>% filter(StatRec %in% WB_rectangles)
hhWB$jday <- yday(paste(hhWB$Year,'-',hhWB$Month,'-',hhWB$Day,sep=""))
survey.day <- aggregate(jday~Year+Quarter,data=hhWB,FUN=mean)
survey.day$jday <- survey.day$jday/365
quarter.day <- aggregate(jday~Quarter,data=survey.day,FUN=mean)
#names(quarter.day)[names(quarter.day)=='jday'] <- 'quarter.day'

# Fix that some lengths are in mm and some in cm
flatfish.ca$LngtClass[flatfish.ca$LngtCode==0] <- round(flatfish.ca$LngtClass[flatfish.ca$LngtCode==0]/10)
flatfish.ca$LngtClass[flatfish.ca$LngtCode=='.'] <- round(flatfish.ca$LngtClass[flatfish.ca$LngtCode=='.']/10)

flatfish.ca <- flatfish.ca %>% filter(LngtClass<90)
flatfish.ca <- flatfish.ca %>% filter(Age<50)

#flatfish.ca <- flatfish.ca %>% left_join(survey.day,by=c('Year','Quarter'))
flatfish.ca <- flatfish.ca %>% left_join(quarter.day,by=c('Quarter'))
flatfish.age <- aggregate(LngtClass~Age+Year+Quarter+species+jday,
                          data=flatfish.ca %>% filter(Age>1 & LngtClass<60),FUN=mean)
flatfish.sd <- aggregate(LngtClass~Age+Year+Quarter+species+jday,
                         data=flatfish.ca %>% filter(Age>1 & LngtClass<60),FUN=sd)
#####
rm(list=setdiff(ls(),c('hLdata','cadata','hhdata','flatfish.age','flatfish.sd','flatfish.ca')))

# Perform analysis
#####
#t_hatch <- 61/365 # time of hatching, Julian day / 365 - march 1st (early spring from Baden et al., 2022)
flatfish.age <- flatfish.age %>% 
  mutate(t_hatch = case_when(species=='plaice' ~ 15/365,
                             species=='flounder' ~ 75/365,
                             species=='dab' ~ 167/365))
#flatfish.age$t_hatch <- rep(61/365,length(flatfish.age$Age))
flatfish.age$fish.age <- flatfish.age$Age+flatfish.age$jday-flatfish.age$t_hatch

L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
u <- 0.22 # Amplitude of seasonality. From McQueen et al., 2019 (table 1) - Assuming similar seasonal growth between cod and flatfish
w <- 0.90 # Timing of peak growth. From McQueen et al., 2019 (table 1)
flatfish.age$X_t <- flatfish.age$jday # time of year [0,1]
flatfish.age$phi_t <- u*sin(2*pi*(flatfish.age$X_t-w))/(2*pi) # seasonal variability in growth
flatfish.age$phi_hatch <- u*sin(2*pi*(flatfish.age$t_hatch-w))/(2*pi) # seasonal variability in growth

#L_t <- (L_hatch-Linf)*exp(-k*(phi_t+t-phi_hatch-t_hatch))+Linf # vbgr estimated length
#flatfish.age <- flatfish.age %>% filter(!(species=='dab' & Age<3))
#fit <- data.frame(x=t,y=L_t)
df.fit <- aggregate(LngtClass~fish.age+X_t+phi_t+phi_hatch+species+t_hatch,data = flatfish.age,FUN=mean)

dab_age <- df.fit %>% filter(species=='dab')

plaice_age <- df.fit %>% filter(species=='plaice')
flounder_age <- df.fit %>% filter(species=='flounder')

fit.flounder <- nls(LngtClass~(L_hatch-Linf)*exp(-k*(phi_t+fish.age-phi_hatch))+Linf, # with seasonality
                  data=flounder_age,start = list(k=0.57,Linf=31))
fit.flounder2 <- nls(LngtClass~(L_hatch-Linf)*exp(-k*(fish.age))+Linf,
                    data=flounder_age,start = list(k=0.57,Linf=31))


fit.plaice <- nls(LngtClass~(L_hatch-Linf)*exp(-k*(phi_t+fish.age-phi_hatch))+Linf, # with seasonality
                  data=plaice_age,start = list(k=0.4,Linf=37))
fit.plaice2 <- nls(LngtClass~(L_hatch-Linf)*exp(-k*(fish.age))+Linf,
                   data=plaice_age,start = list(k=0.4,Linf=37))

fit.dab <- nls(LngtClass~(L_hatch-Linf)*exp(-k*(phi_t+fish.age-phi_hatch))+Linf, # with seasonality
               data=dab_age,start = list(k=0.51,Linf=27))
fit.dab2 <- nls(LngtClass~(L_hatch-Linf)*exp(-k*(fish.age))+Linf,
               data=dab_age,start = list(k=0.51,Linf=27))
fits <- list()
fits$flounder <- fit.flounder
fits$plaice <- fit.plaice
fits$dab <- fit.dab

AIC(fit.flounder)
AIC(fit.flounder2)
AIC(fit.plaice)
AIC(fit.plaice2)
AIC(fit.dab)
AIC(fit.dab2)

t <- (0:2500)/100
X_t <- t-floor(t) # time of year [0,1]
phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
phi_hatch <- c(dab_age$phi_hatch[1],plaice_age$phi_hatch[1],flounder_age$phi_hatch[1]) # Dab, Plaice, Flounder
t_hatch <- c(dab_age$t_hatch[1],plaice_age$t_hatch[1],flounder_age$t_hatch[1]) # Dab, Plaice, Flounder
#phi_hatch <- c(0.03481163,0.03481163,0.03481163) # Dab, Plaice, Flounder
#t_hatch <- c(61/365,61/365,61/365) # Dab, Plaice, Flounder
new_df <- data.frame(t = rep(t,3), X_t=rep(X_t,3),phi_t=rep(phi_t,3),
                     species=rep(c('dab','plaice','flounder'),each=length(t)),
                     fish.age = c(t-dab_age$t_hatch[1],t-plaice_age$t_hatch[1],t-flounder_age$t_hatch[1]),
                     phi_hatch = c(rep(dab_age$phi_hatch[1],length(t)),
                                    rep(plaice_age$phi_hatch[1],length(t)),
                                    rep(flounder_age$phi_hatch[1],length(t))))

dab_pred <- predict(fit.dab,newdata = new_df %>% filter(species=='dab'))
plaice_pred <- predict(fit.plaice,newdata = new_df %>% filter(species=='plaice')) 
flounder_pred <-  predict(fit.flounder,newdata = new_df %>% filter(species=='flounder'))


flatfish.pred <- data.frame(predictions = c(dab_pred,plaice_pred,flounder_pred),
                            age = rep(t,3),
                            species = rep(c('dab','plaice','flounder'),each=length(t)))

ggplot(data=flatfish.age,aes(x=fish.age,y=LngtClass,color=species))+geom_point()+
  facet_wrap(~species)+ylim(0,60)+xlim(0,15)+
  geom_line(data=flatfish.pred, aes(x = age, y = predictions),
            color='black',linewidth=1)+ theme_bw() +geom_abline(slope = 0,intercept = L_hatch)+
  labs(x='Age [years]',y='Fish length [cm]')+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')

#####
rm(list=setdiff(ls(),c('hLdata','cadata','hhdata','flatfish.age','flatfish.ca',
                       'flatfish.pred','fits')))

# Estimate variability
#####
# First find the sd pr. age for each sampling
flatfish.ca <- flatfish.ca %>% 
  mutate(t_hatch = case_when(species=='plaice' ~ 15/365,
                             species=='flounder' ~ 75/365,
                             species=='dab' ~ 167/365)) # Hatching times from FiskeAtlas

flatfish.ca$fish.age <- flatfish.ca$Age+flatfish.ca$jday-flatfish.ca$t_hatch
     
flatfish.ca <- flatfish.ca %>% filter(Age<30 & LngtClass<80 & Age>1)
flatfish.ca$qy <- paste(flatfish.ca$Quarter,flatfish.ca$Year)
qy.s <- unique(flatfish.ca$qy)
sd.s <- data.frame(fish.age = 0,species = '0',
                   Year = 0,qy = 0, LngtClass = 0)
for (i in 1:length(qy.s)){
  dada <- flatfish.ca %>% filter(qy == qy.s[i])
  std <- aggregate(LngtClass~fish.age+species+Year+qy,data = dada,FUN=sd)
  sd.s <- rbind(sd.s,std)
}
sd.s <- sd.s %>% filter(!(is.na(LngtClass)) & LngtClass>0)
sd.mean <- aggregate(LngtClass~fish.age+species,data = sd.s,FUN=mean)
sd.mean$n <- c(table(sd.s$fish.age[sd.s$species=='dab']),
               table(sd.s$fish.age[sd.s$species=='flounder']),
               table(sd.s$fish.age[sd.s$species=='plaice']))
flounder.mean <- sd.mean %>% filter(species=='flounder')
plaice.mean <- sd.mean %>% filter(species=='plaice')
dab.mean <- sd.mean %>% filter(species=='dab')

# Now find the variability in the parameters that fits the overall observed deviation
flatfish.age$fish.age <- flatfish.age$Age+flatfish.age$jday-flatfish.age$t_hatch
flatfish.age <- flatfish.age %>% filter(Age<30 & LngtClass<80)
L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
u <- 0.22 # Amplitude of seasonality. From McQueen et al., 2019 (table 1) - Assuming similar seasonal growth between cod and flatfish
w <- 0.90 # Timing of peak growth. From McQueen et al., 2019 (table 1)
flatfish.age$X_t <- flatfish.age$jday # time of year [0,1]
flatfish.age$phi_t <- u*sin(2*pi*(flatfish.age$X_t-w))/(2*pi) # seasonal variability in growth
flatfish.age$phi_hatch <- u*sin(2*pi*(flatfish.age$t_hatch-w))/(2*pi) # seasonal variability in growth

# Function to calculate sd of L as a function of t. Based on partial derivatives of L - 
calculate_sd <- function(t,Linf_sd, k_sd,species) {
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  t_hatch <- flatfish.age$t_hatch[flatfish.age$species==species][1]
  fish.age <- t-t_hatch
  u <- 0.15 # Amplitude of seasonality
  w <- 0.87 # Timing of peak growth
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  k <- coef(fits[[species]])[1]
  Linf <- coef(fits[[species]])[2]

  # Partial derivatives
  dL_dLinf <- -exp(-k * (phi_t - phi_hatch +fish.age)) + 1
  dL_dk <- (L_hatch - Linf) * (phi_hatch - phi_t -fish.age ) * exp(-k * (phi_t +fish.age - phi_hatch))
  dL_dage <- -(L_hatch - Linf)*k*exp(-k*(phi_t + fish.age - phi_hatch))
 
  # Partial derivatives at age 0
  dL_dLinf.0 <- -exp(-k * (phi_hatch - phi_hatch +0)) + 1
  dL_dk.0 <- (L_hatch - Linf) * ( -0 ) * exp(-k * (0))
  dL_dage.0 <- -(L_hatch - Linf)*k*exp(-k*(0))
  
  
  sd.age <- data.frame(species=c('flounder','plaice','dab'),
                       sd = c(0.05,0.08,0.08))
  age_sd <- sd.age$sd[sd.age$species==species]/dL_dage.0
  # Variance
  var_L <- (dL_dLinf * Linf_sd)^2 + (dL_dk * k_sd)^2 + (dL_dage * age_sd)^2
 
  # Standard deviation
  sd_L <- sqrt((dL_dLinf * Linf_sd)^2 +(dL_dk * k_sd)^2 +(dL_dage * age_sd)^2)
  
  return(sd_L)
} #Original
#par(mfrow=c(2,2))
#plot(fish.age[fish.age>0],L[fish.age>0],ylab = 'L')
#plot(fish.age[fish.age>0],dL_dLinf[fish.age>0],ylab = 'dL / dL_inf')
#plot(fish.age[fish.age>0],dL_dk[fish.age>0],ylab = 'dL / dk')
#plot(fish.age[fish.age>0],dL_dt_hatch[fish.age>0],ylab = 'dL / dage')

# Fit parameter varibility to observed variability
y.f <- sd.mean$LngtClass[sd.mean$species=='flounder']
x.f <- sd.mean$fish.age[sd.mean$species=='flounder']
n.f <- sd.mean$n[sd.mean$species=='flounder']
flounder.sd <- nlsLM(y.f~calculate_sd(x.f,Linf_sd,k_sd,species='flounder'),
                start = list(Linf_sd=2.5,k_sd=0.05))
par(mfrow=c(1,3))
t.f <- ((0:2500)/100)[((0:2500)/100)>flatfish.age$t_hatch[flatfish.age$species=='flounder'][1]]
plot(t.f, calculate_sd(t.f,Linf_sd=coef(flounder.sd)[1],k_sd=coef(flounder.sd)[2],
                       species='flounder'), type = "l", col = "black", lwd = 2,
     ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,15),ylim=c(0,4))
points(flounder.mean$fish.age,flounder.mean$LngtClass)

y.p <- sd.mean$LngtClass[sd.mean$species=='plaice']
x.p <- sd.mean$fish.age[sd.mean$species=='plaice']
n.p <- sd.mean$n[sd.mean$species=='plaice']
plaice.sd <- nlsLM(y.p~calculate_sd(x.p,Linf_sd,k_sd,species='plaice'),
                     start = list(Linf_sd=2.5,k_sd=0.025))
t.p <- ((0:2500)/100)[((0:2500)/100)>flatfish.age$t_hatch[flatfish.age$species=='plaice'][1]]
plot(t.p, calculate_sd(t.p,Linf_sd=coef(plaice.sd)[1],k_sd=coef(plaice.sd)[2],
                       species='plaice'), type = "l", col = "black", lwd = 2,
     ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,15),ylim=c(0,4))
points(plaice.mean$fish.age,plaice.mean$LngtClass)

y.d <- sd.mean$LngtClass[sd.mean$species=='dab']
x.d <- sd.mean$fish.age[sd.mean$species=='dab']
n.d <- sd.mean$n[sd.mean$species=='dab']
dab.sd <- nlsLM(y.d~calculate_sd(x.d,Linf_sd,k_sd,species='dab'),
                     start = list(Linf_sd=2.5,k_sd=0.05))
t.d <- ((0:2500)/100)[((0:2500)/100)>flatfish.age$t_hatch[flatfish.age$species=='dab'][1]]
plot(t.d, calculate_sd(t.d,Linf_sd=coef(dab.sd)[1],k_sd=coef(dab.sd)[2],
                       species='dab'), type = "l", col = "black", lwd = 2,
     ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,15),ylim=c(0,4))
points(dab.mean$fish.age,dab.mean$LngtClass)


par(mfrow=c(1,1))
t.f <- ((0:2500)/100)[((0:2500)/100)>flatfish.age$t_hatch[flatfish.age$species=='flounder'][1]]
t.p <- ((0:2500)/100)[((0:2500)/100)>flatfish.age$t_hatch[flatfish.age$species=='plaice'][1]]
t.d <- ((0:2500)/100)[((0:2500)/100)>flatfish.age$t_hatch[flatfish.age$species=='dab'][1]]
plot(t.f, calculate_sd(t.f,Linf_sd=coef(flounder.sd)[1],k_sd=coef(flounder.sd)[2],
                      species='flounder'), type = "l", col = "black", lwd = 2,
     ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,15),ylim=c(0,4))
lines(t.p, calculate_sd(t.p,Linf_sd=coef(plaice.sd)[1],k_sd=coef(plaice.sd)[2],
                     species='plaice'), type = "l", col = "blue", lwd = 2)
lines(t.d, calculate_sd(t.d,Linf_sd=coef(dab.sd)[1],k_sd=coef(dab.sd)[2],
                      species='dab'), type = "l", col = "orange", lwd = 2)
legend(10,1.5,legend = c('Flounder','Plaice','Dab'),col=c('black','blue','orange'),pch=19)
#####




