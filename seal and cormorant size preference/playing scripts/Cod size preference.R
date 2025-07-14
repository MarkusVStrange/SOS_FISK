setwd("/Users/mavast/Desktop/Markus SOS")
library(sicegar)
library(ggplot2)
library(lubridate)
library(minpack.lm)
set.seed(1)
##############################
# Cod size preference based on field observations
##############################
# Growth model, cohort spread by length, and individuals at by age and time, functions: vbgr.fixed(), vbgr.sd(), N_age
#####
vbgr <- function(age){ # von Bertalanffy growth rate. Parameters from McQueen et al., 2019 with fixed hatching time
  # define growth parameters
  Linf <- 97.9 # cm. McQueen et al., 2019
  L_hatch <- 0.4 # length at hatching, cm (Pepin et al, 1997)
  k <- 0.22 # vbgr growth rate, from McQueen
  u <- 0.23 # Amplitude of seasonality
  w <- 0.91 # Timing of peak growth
  t_hatch <- ((46+135)/2)/365 # time of hatching, Julian day / 365 - spawning from Feb. to May (Fiskeatlas), average set to April 1st
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf # vbgr estimated length
  return(L_t)
}


vbgr.sd <- function(age,n){ # von Bertalanffy growth rate. Parameters from McQueen et al., 2019 and estimated
  # t in fraction year, Julian day/365
  Linf <- rnorm(n,mean=97.9,sd=19.92531) # Assymptotic length - with stochasticity
  #Linf <- 154.56 # cm
  L_hatch <- 0.4 # length at hatching, cm (Pepin et al, 1997)
  t_hatch <- ((46+135)/2)/365
  k <- rnorm(n,mean=0.22,sd=4.689802e-07) # vbgr growth rate
  u <- 0.23 # Amplitude of seasonality
  w <- 0.91 # Timing of peak growth
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  sd <- rep(0,length(t))
  for (i in 1:length(sd)){
    fish.age <- age[i]+rnorm(n,mean=0,sd=4.862846e-02) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
    L_t <- (L_hatch-Linf)*exp(-k*(phi_t[i]+fish.age-phi_hatch))+Linf # vbgr estimated length
    sd[i] <- sd(L_t)
  }
  return(sd)
}

N_age <- function(dat,age,yrs){ # Number of individuals by age
  M_age$Z <- M_age$M_var+M_age$Fest
  Ns <- matrix(rep(0,(length(age)+1)*length(yrs)*365),ncol = length(age)+1)
  Ns[,length(age)+1] <- rep(yrs,each=365)
  x <- seq(1/365,1,length=365)
  for(i in 1:length(yrs)){
    for(j in 1:length(age)){
      di <- M_age %>% filter(ages==age[j] & years==yrs[i])
      dip1 <- M_age %>% filter(ages==(age[j]+1) & years==(yrs[i]+1))
      
      Ns[(365*(i-1)+1):(365*i),j] <- di$N*exp(-x*di$Z)
    }
  }
  return(Ns)
}
#####

#read and prepare cormorant and cod data
#####
setwd("/Users/mavast/Desktop/Cormorants")
d<- read.table("90s_cormorant_pellets.csv",header=TRUE,
               sep=";",as.is=TRUE)
info <- read.table("90s_sampling_information.csv",header=TRUE,
                   sep=";",as.is=TRUE)

month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
           "Sep","Oct","Nov","Dec")
# Fix species names
d$ART[d$ART=="SAN" | d$ART=="SOR"] <- rep('goby',length(d$ART[d$ART=="SAN" | d$ART=="SOR"]))
d$ART[d$ART=="TOR"] <- rep('cod',length(d$ART[d$ART=="TOR"]))
d$ART[d$ART=="ISI" | d$ART=="SKB" | d$ART=="R\xd8S"] <- rep('flatfish',length(d$ART[d$ART=="ISI" | d$ART=="SKB" | d$ART=="R\xd8S"]))
d$ART[d$ART=="\xc5LA"] <- rep('eelpout',length(d$ART[d$ART=="\xc5LA"]))
d$ART[d$ART=="ULK"] <- rep('sculpin',length(d$ART[d$ART=="ULK"]))
d$ART[d$ART=="TBI"] <- rep('sandeel',length(d$ART[d$ART=="TBI"]))
d$ART[d$ART=="SIL"] <- rep('herring',length(d$ART[d$ART=="SIL"]))
d$ART[-which(d$ART=="goby" | d$ART=="cod" | d$ART=="flatfish" | d$ART=="eelpout" |
               d$ART=="sculpin" | d$ART=="sandeel" | d$ART=="herring")] <- 
  rep('other',length(d$ART[-which(d$ART=="goby" | d$ART=="cod" | d$ART=="flatfish" | d$ART=="eelpout" |
                                    d$ART=="sculpin" | d$ART=="sandeel" | d$ART=="herring")]))


d$day <- rep(0,length(d$GYLP))
d$month <- rep(0,length(d$GYLP))
d$year <- rep(0,length(d$GYLP))
d$dmyl <- rep(0,length(d$GYLP))
d$location <- rep(0,length(d$GYLP))
#d$rel_abu <- rep(0,length(d$GYLP))
d$rel_bio <- rep(0,length(d$GYLP))
d <- d[-which(d$FISKLGD==0 & d$FISKVGT==0),]
d <- d[-which(is.na(d$FISKLGD) & is.na(d$FISKVGT)),]
d$size <- d$FISKLGD
d$size <- floor(d$size/10)*10

for (i in 1:length(d$GYLP)){
  gylp <- d$GYLP[i]
  inf <- subset(info,info$GYLP==gylp)
  
  d$day[i] <- inf$DAG
  d$month[i] <- month[inf$MON]
  d$year[i] <- inf$YEAR+1900
  d$location[i] <- inf$LOKAL
  d$dmyl[i] <- paste(as.character(inf$DAG),as.character(inf$MON),as.character(inf$YEAR),inf$LOKAL)
  
}


d <- d[-which(is.na(d$size)),]
#Fix location names
d$location[d$location=="BR\xc6"] <- rep("BRA",length(d$location[d$location=="BR\xc6"]))
d$location[d$location=="M\xc5\xd8"] <- rep("MAO",length(d$location[d$location=="M\xc5\xd8"]))
write.table(d,file = "/Users/mavast/Desktop/Markus SOS/seal and cormorant size preference/data/corm_diet.csv",sep = ';',dec = '.')

cod_corm <- subset(d,ART=='cod')
setwd("/Users/mavast/Desktop/Markus SOS")
M_age <- read.table('M_sms.csv',sep=';',header=T)
#rm(list=setdiff(ls(),c('cod_corm','M_age','vbgr.fixed','vbgr.sd','N_age')))
#####

# read and prepare seal data
#####
setwd("/Users/mavast/Desktop/Seals")
D <- read.table("scat_otoliths.csv",header=TRUE,
                sep=";",as.is=TRUE)
D$FL.with.SCF <- as.numeric(D$FL.with.SCF)
D$FW.with.SCF <- as.numeric(D$FW.with.SCF)
D$FL.with.SCF[which(D$Species=="cod" | D$Species=="Cod" |
                      D$Species=="cod ")] <- D$FL.with.SCF[which(D$Species=="cod" | D$Species=="Cod" |
                                                                   D$Species=="cod ")]/10
D$FL.with.SCF <- D$FL.with.SCF*10

D$size <- floor(D$FL.with.SCF/10)*10
D <- D[-which(is.na(D$size)),]
D <- D[-which(D$size==0),]
D$Month <- stringr::str_extract(D$Month, "^.{3}")

D$site_abb <- rep("O",length(D$Site))
#Abbreviate locations
D$site_abb[which(D$Site=="Måkläppen")] <- rep("Mk",length(D$Site[which(D$Site=="Måkläppen")]))
D$site_abb[which(D$Site=="Rødsand")] <- rep("Rs",length(D$Site[which(D$Site=="Rødsand")]))
D$site_abb[which(D$Site=="Tat")] <- rep("T",length(D$Site[which(D$Site=="Tat")]))
D$site_abb[which(D$Site=="Utklippan")] <- rep("Uk",length(D$Site[which(D$Site=="Utklippan")]))


table(D$Species)
cod <- D[which(D$Species=="cod" | D$Species=="Cod" |
                 D$Species=="cod "),]
cod$Species <- rep("cod",length(cod$Scat))
flatfish <- D[which(D$Species=="flatfish" | D$Species=="Flatfish"),]
flatfish$Species <- rep("flatfish",length(flatfish$Scat))
goby <- D[which(D$Species=="GOBY" | D$Species=="goby"),]
goby$Species <- rep("goby",length(goby$Scat))
herring <- D[which(D$Species=="herring" | D$Species=="Herring"),]
herring$Species <- rep("herring",length(herring$Scat))
sandeel <- D[which(D$Species=="sandeel"),]
sandeel$Species <- rep("sandeel",length(sandeel$Scat))
sprat <- D[which(D$Species=="sprat" | D$Species=="Sprat"),]
sprat$Species <- rep("sprat",length(sprat$Scat))
whiting <- D[which(D$Species=="whiting" | D$Species=="Whiting"),]
whiting$Species <- rep("whiting",length(whiting$Scat))
other <- D[which(D$Species=="four-bearded rockling" | D$Species=="Unidentified"),]
other$Species <- rep("other",length(other$Scat))
Dall <- rbind(cod,herring,flatfish,goby,sandeel,sprat,whiting)
Dall$yms <- paste(Dall$Year,Dall$Month,Dall$Site)

fish <- c("cod","herring","flatfish","goby",
          "sandeel","sprat","whiting","other")

# plaice NCF used for flatfish - assuming flatfish are largely plaice and flounder, and that NCF's are similar between these two species
# flatfish NCF used for goby (Karin Hüssy). Sandeel NCF used for sprat. NB find a better solutions
NCF <- data.frame(species = fish, ncf = c(1.2,3,1.6,1.6,3.6,3.6,1.3,1)) 

mar <- aggregate(FW.with.SCF~Year+Month+Site+Species,data=Dall,FUN = sum)
mar$ymp <- paste(mar$Year,mar$Month,mar$Site)

mar$cor <- rep(0,length(mar$ymp))
mar$rel_biom <- rep(0,length(mar$ymp))
# with NCF'
for (i in 1:length(mar$ymp)){
  idx <- mar$ymp[i]
  sp <- mar$Species[i]
  mar$cor[i] <- mar$FW.with.SCF[i]*NCF$ncf[NCF$species==sp]
}

mar_tot <- aggregate(cor~Year+Month+Site+ymp,data=mar,FUN = sum)
for (i in 1:length(mar$ymp)){
  idx <- mar$ymp[i]
  mar$rel_biom[i] <- mar$cor[i]/mar_tot$cor[mar_tot$ymp==idx]
  
}
cod_seal <- cod[,c('Scat','Site','Year','Month','size')]
write.table(cod_seal,file = "/Users/mavast/Desktop/Markus SOS/seal and cormorant size preference/data/gSeal_diet_diet.csv",sep = ';',dec = '.')

#rm(list=setdiff(ls(),c('cod_corm','M_age','vbgr.fixed','vbgr.sd','N_age','cod_seal')))
#####

# Seal preference for cod
#####
#Ns <- N_age(M_age,c(0,1,2,3),c(1992,1993,1994))
cod_seal$my <- paste(cod_seal$Month,'-',cod_seal$Year)


samplings <- unique(cod_seal$my)
month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
           "Sep","Oct","Nov","Dec") 

lng.class <- data.frame(l.class = 10*(0:150))

df.seal <- data.frame(l.class = rep(10*(0:150),length(samplings)),
                      prop_eaten = rep(0,151*length(samplings)),
                      eaten = rep(0,151*length(samplings)),
                      not_eaten = rep(0,151*length(samplings)),
                      n_scat = rep(0,151*length(samplings)),
                      year = rep(0,151*length(samplings)),
                      location = rep(0,151*length(samplings)))


for (i in 1:length(samplings)){
  dada <- cod_seal %>% filter(my==samplings[i])
  n_sca <- length(unique(dada$Scat))
  date <- paste(dada$Year[1],'-',which(dada$Month[1]==month)[[1]],'-',15,sep="")
  jday <- yday(date)
  sizes <- c(vbgr.fixed(jday/365),vbgr.fixed(jday/365+1),vbgr.fixed(jday/365+2),vbgr.fixed(jday/365+3),
             vbgr.fixed(jday/365+4),vbgr.fixed(jday/365+5),vbgr.fixed(jday/365+6),vbgr.fixed(jday/365+7))*10
  ages <- c(0,1,2,3,4,5,6,7)
  pop <- N_age(M_age,ages,dada$Year[1])[jday,1:length(ages)]
  sd<- c(vbgr.sd(jday/365,100000),vbgr.sd(jday/365+1,100000),vbgr.sd(jday/365+2,100000),vbgr.sd(jday/365+3,100000),
         vbgr.sd(jday/365+4,100000),vbgr.sd(jday/365+5,100000),vbgr.sd(jday/365+6,100000),vbgr.sd(jday/365+7,100000))*10
  
  if (min(sd)<0){
    sd[sd<0] <-rep(0,length(sd[sd<0]))
  }
  available <- c(floor(rnorm(pop[1],sizes[1],sd[1])/10)*10,floor(rnorm(pop[2],sizes[2],sd[2])/10)*10,
                 floor(rnorm(pop[3],sizes[3],sd[3])/10)*10,floor(rnorm(pop[4],sizes[4],sd[4])/10)*10,
                 floor(rnorm(pop[5],sizes[5],sd[5])/10)*10,floor(rnorm(pop[6],sizes[6],sd[6])/10)*10,
                 floor(rnorm(pop[7],sizes[7],sd[7])/10)*10,floor(rnorm(pop[8],sizes[8],sd[8])/10)*10)
  available[available<0] <- rep(0,length(available[available<0]))
  
  avail <- as.data.frame(table(available))
  colnames(avail) <- c("l.class", "available")
  eat <- as.data.frame(table(dada$size))
  colnames(eat) <- c("l.class", "eaten")
  
  df <- merge(lng.class,avail,by='l.class',all.x = T)
  df <- merge(df,eat,by='l.class',all.x = T)
  df[is.na(df)] <- rep(0,length(df[is.na(df)]))
  
  df$eaten <- df$eaten
  df$not_eaten <- df$available-df$eaten
  df$prop_eaten <- (df$eaten/df$available)
  df[is.na(df)] <- rep(0,length(df[is.na(df)]))
  df$prop_eaten <- df$prop_eaten/sum(df$prop_eaten)
  
  df$n_scat <- rep(n_sca,151)
  df$location <- rep(dada$Site[1],151)
  df$year <- rep(dada$Year[1],151)
  
  df.seal[((i-1)*151+1):(i*151),c('prop_eaten','eaten','not_eaten','n_scat','year','location')] <- df[,c('prop_eaten','eaten','not_eaten','n_scat','year','location')]
}

df.seal.ori <- df.seal
df.seal <- df.seal.ori %>% filter(!is.na(prop_eaten))
df.seal <- df.seal %>% filter(abs(prop_eaten)<1)

# 2015

df.15 <- df.seal %>% filter(year==2015) 
fac.15 <- sum(sort(unique(df.15$n_scat))*table(df.15$n_scat)/151)
mean(df.15$prop_eaten)
df.15$weighted.prop <- df.15$prop_eaten*df.15$n_scat
df.fin.15 <- aggregate(weighted.prop~l.class,data=df.15,FUN=mean)
df.fin.15$weighted.prop <- df.fin.15$weighted.prop*table(df.15$l.class)[1]/fac.15
mean(df.fin.15$weighted.prop)

plot(df.fin.15$l.class,df.fin.15$weighted.prop,xlim=c(0,1000),pch=19,
     xlab = 'cm',ylab = 'preference',col='blue')
# 2016
df.16 <- df.seal %>% filter(year==2016) 
fac.16 <- sum(sort(unique(df.16$n_scat))*table(df.16$n_scat)/151)
mean(df.16$prop_eaten)
df.16$weighted.prop <- df.16$prop_eaten*df.16$n_scat
df.fin.16 <- aggregate(weighted.prop~l.class,data=df.16,FUN=mean)
df.fin.16$weighted.prop <- df.fin.16$weighted.prop*table(df.16$l.class)[1]/fac.16
mean(df.fin.16$weighted.prop)

points(df.fin.16$l.class,df.fin.16$weighted.prop,pch=19,col='red')
# 2017
df.17 <- df.seal %>% filter(year==2017) 
fac.17 <- sum(sort(unique(df.17$n_scat))*table(df.17$n_scat)/151)
mean(df.17$prop_eaten)
df.17$weighted.prop <- df.17$prop_eaten*df.17$n_scat
df.fin.17 <- aggregate(weighted.prop~l.class,data=df.17,FUN=mean)
df.fin.17$weighted.prop <- df.fin.17$weighted.prop*table(df.17$l.class)[1]/fac.17
mean(df.fin.17$weighted.prop)

points(df.fin.17$l.class,df.fin.17$weighted.prop,pch=19,col='green')
# All together - df.fin
fac <- sum(sort(unique(df.seal$n_scat))*table(df.seal$n_scat)/151)
mean(df.seal$prop_eaten)
df.seal$weighted.prop <- df.seal$prop_eaten*df.seal$n_scat
df.seal.fin <- aggregate(weighted.prop~l.class,data=df.seal,FUN=mean)
df.seal.fin$weighted.prop <- df.seal.fin$weighted.prop*table(df.seal$l.class)[1]/fac
mean(df.seal.fin$weighted.prop)
points(df.seal.fin$l.class,df.seal.fin$weighted.prop,pch=19)
legend(0,0.12,legend = c('1992','1993','1994','All'),pch=19,
       col=c('blue','red','green','black'))
#####

#seal fit cod - run multiple times until lowest AIC (~-120). NB improve this method
#####
dataInput <- df.seal.fin %>% filter(l.class<700)
dataInput <- dataInput[-which(dataInput$l.class>100 & dataInput$l.class<560 & dataInput$weighted.prop==0),]

names(dataInput)[names(dataInput)=="l.class"] <- "time"
names(dataInput)[names(dataInput)=="weighted.prop"] <- "intensity"
time <- dataInput$time
for (i in 1:1000){
  normalizedInput <- normalizeData(dataInput)
  parameterVector <- doublesigmoidalFitFunction(normalizedInput,
                                             tryCounter = 200,n_iterations = 1024)
    if (parameterVector$isThisaFit==T & parameterVector$AIC_value<(-76)){
    break
  } 
}  
#parameterVector <- doublesigmoidalFitFunction(normalizedInput,
#                                              tryCounter = 200,n_iterations = 1024)
#parameterVector$AIC_value
#Check the results
if(parameterVector$isThisaFit){
  intensityTheoretical <-
    doublesigmoidalFitFormula(
      time,
      finalAsymptoteIntensityRatio = parameterVector$finalAsymptoteIntensityRatio_Estimate,
      maximum = parameterVector$maximum_Estimate,
      slope1Param = parameterVector$slope1Param_Estimate,
      midPoint1Param = parameterVector$midPoint1Param_Estimate,
      slope2Param = parameterVector$slope2Param_Estimate,
      midPointDistanceParam = parameterVector$midPointDistanceParam_Estimate)
  
  comparisonData <- cbind(dataInput, intensityTheoretical)
  require(ggplot2)
  ggplot(comparisonData) +
    geom_point(aes(x = time, y = intensity)) +
    geom_line(aes(x = time, y = intensityTheoretical), color = "orange") +
    expand_limits(x = 0, y = 0)
}




double_sigmoid <- function(x, maximum, slope1, midPoint1, slope2, midPointDistance) {
  0 * maximum + 
    (maximum / (1 + exp(-slope1 * (x - midPoint1)))) - 
    (maximum / (1 + exp(-slope2 * (x - (midPoint1 + midPointDistance)))))
}
df <- df.seal.fin %>% filter(l.class<700)
df <- df[-which(df$l.class>100 & df$l.class<560 & df$weighted.prop==0),]

fit.seal <- nlsLM(
  weighted.prop ~ double_sigmoid(l.class, maximum, slope1, midPoint1, slope2, midPointDistance),
  data = df,control = nls.lm.control(maxiter = 1024),
  start = list(
    maximum = 0.05,
    slope1 = 0.01,
    midPoint1 = 250,
    slope2 = 0.003,
    midPointDistance = 200
  )
)



prediction <- predict(fit.seal,newdata = df,type='response')
prediction[prediction<0] <- rep(0,length(prediction[prediction<0]))

seal_fit <- cbind(time, intensityTheoretical)
plot(df$l.class, prediction,type='l',lwd=2,xlim=c(0,600))
points(df$l.class, df$weighted.prop)
lines(df$l.class,intensityTheoretical,lwd=2,col='red')
sqrt(mean((prediction-df$weighted.prop)^2))
sqrt(mean((intensityTheoretical-df$weighted.prop)^2))


#rm(list=setdiff(ls(),c('cod_corm','M_age','vbgr.fixed','vbgr.sd','N_age','cod_seal','df.seal.fin','seal_fit')))
#####

# Cormorant preference for cod
#####
cod_corm$dmy <- paste(cod_corm$year,'-',cod_corm$month,'-',cod_corm$day,sep="")


samplings <- unique(cod_corm$dmy)
month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
           "Sep","Oct","Nov","Dec") 

lng.class <- data.frame(l.class = 10*(0:150))

df.corm <- data.frame(l.class = rep(10*(0:150),length(samplings)),
                     prop_eaten = rep(0,151*length(samplings)),
                     eaten = rep(0,151*length(samplings)),
                     not_eaten = rep(0,151*length(samplings)),
                     n_pellet = rep(0,151*length(samplings)),
                     year = rep(0,151*length(samplings)),
                     location = rep(0,151*length(samplings)))


for (i in 1:length(samplings)){
  dada <- cod_corm %>% filter(dmy==samplings[i])
  n_pel <- length(unique(dada$GYLP))
  date <- paste(dada$year[1],'-',which(dada$month[1]==month)[[1]],'-',dada$day[1],sep="")
  jday <- yday(date)
  sizes <- c(vbgr.fixed(jday/365),vbgr.fixed(jday/365+1),vbgr.fixed(jday/365+2),vbgr.fixed(jday/365+3),
             vbgr.fixed(jday/365+4),vbgr.fixed(jday/365+5),vbgr.fixed(jday/365+6),vbgr.fixed(jday/365+7))*10
  ages <- c(0,1,2,3,4,5,6,7)
  pop <- N_age(M_age,ages,dada$year[1])[jday,1:length(ages)]
  sd<- c(vbgr.sd(jday/365,100000),vbgr.sd(jday/365+1,100000),vbgr.sd(jday/365+2,100000),vbgr.sd(jday/365+3,100000),
         vbgr.sd(jday/365+4,100000),vbgr.sd(jday/365+5,100000),vbgr.sd(jday/365+6,100000),vbgr.sd(jday/365+7,100000))*10
  
  if (min(sd)<0){
    sd[sd<0] <-rep(0,length(sd[sd<0]))
  }
  available <- c(floor(rnorm(pop[1],sizes[1],sd[1])/10)*10,floor(rnorm(pop[2],sizes[2],sd[2])/10)*10,
                 floor(rnorm(pop[3],sizes[3],sd[3])/10)*10,floor(rnorm(pop[4],sizes[4],sd[4])/10)*10,
                 floor(rnorm(pop[5],sizes[5],sd[5])/10)*10,floor(rnorm(pop[6],sizes[6],sd[6])/10)*10,
                 floor(rnorm(pop[7],sizes[7],sd[7])/10)*10,floor(rnorm(pop[8],sizes[8],sd[8])/10)*10)
  available[available<0] <- rep(0,length(available[available<0]))
  
  avail <- as.data.frame(table(available))
  colnames(avail) <- c("l.class", "available")
  eat <- as.data.frame(table(dada$size))
  colnames(eat) <- c("l.class", "eaten")
  
  df <- merge(lng.class,avail,by='l.class',all.x = T)
  df <- merge(df,eat,by='l.class',all.x = T)
  df[is.na(df)] <- rep(0,length(df[is.na(df)]))
  
  df$eaten <- df$eaten
  df$not_eaten <- df$available-df$eaten
  df$prop_eaten <- (df$eaten/df$available)
  df[is.na(df)] <- rep(0,length(df[is.na(df)]))
  df$prop_eaten <- df$prop_eaten/sum(df$prop_eaten)
  
  df$n_pellet <- rep(n_pel,151)
  df$location <- rep(dada$location[1],151)
  df$year <- rep(dada$year[1],151)
  
  df.corm[((i-1)*151+1):(i*151),c('prop_eaten','eaten','not_eaten','n_pellet','year','location')] <- df[,c('prop_eaten','eaten','not_eaten','n_pellet','year','location')]
}
#df <- data.frame(x = c('h','h', 'f','f','o','o'),
#                 y = c(4,5,3,3,3.5,3.5))
#barplot(df$x,df$y)

#hist(available,breaks=100,xlim=c(0,600),main = 'Cod population',yaxt='n',xlab='',ylab='Frequency [ind.]')
#hist(dada$size,breaks=30,xlim=c(0,600),main = 'Cod in the diet',yaxt='n',xlab='Length [mm]',ylab='Frequency [ind.]')
# 1992
df.92 <- df.corm %>% filter(year==1992) 
fac.92 <- sum(sort(unique(df.92$n_pellet))*table(df.92$n_pellet)/151)
mean(df.92$prop_eaten)
df.92$weighted.prop <- df.92$prop_eaten*df.92$n_pellet
df.corm.fin.92 <- aggregate(weighted.prop~l.class,data=df.92,FUN=mean)
df.corm.fin.92$weighted.prop <- df.corm.fin.92$weighted.prop*table(df.92$l.class)[1]/fac.92
mean(df.corm.fin.92$weighted.prop)
plot(df.corm.fin.92$l.class,df.corm.fin.92$weighted.prop,xlim=c(0,500),pch=19,
     xlab = 'cm',ylab = 'preference',col='blue')
# 1993
df.93 <- df.corm %>% filter(year==1993) 
fac.93 <- sum(sort(unique(df.93$n_pellet))*table(df.93$n_pellet)/151)
mean(df.93$prop_eaten)
df.93$weighted.prop <- df.93$prop_eaten*df.93$n_pellet
df.corm.fin.93 <- aggregate(weighted.prop~l.class,data=df.93,FUN=mean)
df.corm.fin.93$weighted.prop <- df.corm.fin.93$weighted.prop*table(df.93$l.class)[1]/fac.93
mean(df.corm.fin.93$weighted.prop)
points(df.corm.fin.93$l.class,df.corm.fin.93$weighted.prop,pch=19,col='red')
# 1994
df.94 <- df.corm %>% filter(year==1994) 
fac.94 <- sum(sort(unique(df.94$n_pellet))*table(df.94$n_pellet)/151)
mean(df.94$prop_eaten)
df.94$weighted.prop <- df.94$prop_eaten*df.94$n_pellet
df.corm.fin.94 <- aggregate(weighted.prop~l.class,data=df.94,FUN=mean)
df.corm.fin.94$weighted.prop <- df.corm.fin.94$weighted.prop*table(df.94$l.class)[1]/fac.94
mean(df.corm.fin.94$weighted.prop)
points(df.corm.fin.94$l.class,df.corm.fin.94$weighted.prop,pch=19,col='green')
# All together - df.corm.fin
fac <- sum(sort(unique(df.corm$n_pellet))*table(df.corm$n_pellet)/151)
mean(df.corm$prop_eaten)
df.corm$weighted.prop <- df.corm$prop_eaten*df.corm$n_pellet
df.corm.fin <- aggregate(weighted.prop~l.class,data=df.corm,FUN=mean)
df.corm.fin$weighted.prop <- df.corm.fin$weighted.prop*54/fac
mean(df.corm.fin$weighted.prop)
points(df.corm.fin$l.class,df.corm.fin$weighted.prop,pch=19)
legend(0,0.12,legend = c('1992','1993','1994','All'),pch=19,
       col=c('blue','red','green','black'))
#####

#cormorant fit cod - run multiple times until lowest AIC (~-110). NB improve this method
#####
dataInput <- df.corm.fin %>% filter(l.class<600 & weighted.prop<0.05) #NB removes two outliers
names(dataInput)[names(dataInput)=="l.class"] <- "time"
names(dataInput)[names(dataInput)=="weighted.prop"] <- "intensity"
time <- dataInput$time
for (i in 1:1000){
  normalizedInput <- normalizeData(dataInput)
  parameterVector <- doublesigmoidalFitFunction(normalizedInput,
                                                tryCounter = 200,n_iterations = 1024)
  if (parameterVector$isThisaFit==T & parameterVector$AIC_value<(-33)){
    break
  } 
}
#parameterVector <- doublesigmoidalFitFunction(normalizedInput,
#                                              tryCounter = 200,n_iterations = 1024)
#parameterVector$AIC_value
#Check the results
if(parameterVector$isThisaFit){
  intensityTheoretical <-
    doublesigmoidalFitFormula(
      time,
      finalAsymptoteIntensityRatio = parameterVector$finalAsymptoteIntensityRatio_Estimate,
      maximum = parameterVector$maximum_Estimate,
      slope1Param = parameterVector$slope1Param_Estimate,
      midPoint1Param = parameterVector$midPoint1Param_Estimate,
      slope2Param = parameterVector$slope2Param_Estimate,
      midPointDistanceParam = parameterVector$midPointDistanceParam_Estimate)
  
  comparisonData <- cbind(dataInput, intensityTheoretical)
  require(ggplot2)
  ggplot(comparisonData) +
    geom_point(aes(x = time, y = intensity)) +
    geom_line(aes(x = time, y = intensityTheoretical), color = "orange") +
    expand_limits(x = 0, y = 0)
}
# custom fit
#corm_pref <- function(x, a, b, d, k1, k2) {
#  ifelse(x <= 225, 1 / (d + exp(-k1 * (x - a))),
#         1 / (d + exp(k2 * (x - b))))
#}
#x <- (0:5000)/10
#plot(df.fin$l.class,df.fin$weighted.prop,pch=19,xlim=c(0,500))
#lines(x,corm_pref(x,90,365,40,0.11,0.06),
#      lwd=2,col='orange')
#lines(time,intensityTheoretical,pch=19,col='red',lwd=2,)

cormorant_fit <- cbind(time, intensityTheoretical)
#rm(list=setdiff(ls(),c('cod_corm','M_age','vbgr.fixed',
#                       'vbgr.sd','N_age','cod_seal','df.seal.fin',
#                       'seal_fit','cormorant_fit','df.corm.fin')))

#####

# Seal and cormorant cod preference
#double sigmoidal fit
#####
par(mfrow=c(2,1))
par(mar = c(2.5,4,1,1))
plot(df.corm.fin$l.class,df.corm.fin$weighted.prop,pch=19,
     main = 'Cormorants',xlim=c(0,600),
     ylab='preference',xlab = '', xaxt='n')
lines(cormorant_fit[,1],cormorant_fit[,2]/sum(cormorant_fit[,2]),
      type='l',col='black',lwd=3)
par(mar = c(4,4,1,1))


plot(df.seal.fin$l.class,df.seal.fin$weighted.prop,pch=19,
     main = 'Seals',xlim=c(0,600),
     ylab='preference',xlab = 'cod length [mm]')
lines(seal_fit[,1],seal_fit[,2]/sum(seal_fit[,2]),
      type='l',col='grey60',lwd=3)

par(mfrow=c(1,1))
par(mar = c(4.5,4,1,2))
plot(seal_fit[,1],seal_fit[,2]/sum(seal_fit[,2]),pch=19,col='orange',
     type='l',lwd=2,xlab = 'cod length [mm]',ylab = 'Relative preference')
lines(cormorant_fit[,1],cormorant_fit[,2]/sum(cormorant_fit[,2]),
      type='l',col='red',lwd=2)
legend(500,0.04,legend=c('Seal','Cormorant'),
       pch=19,col=c('orange','red'))
#####


