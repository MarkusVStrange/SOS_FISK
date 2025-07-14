setwd("/Users/mavast/Desktop/Markus SOS")
library(sicegar)
library(ggplot2)
library(minpack.lm)
library(drc)

##############################
# Herring size preference based on field observations
##############################
# Growth model, cohort spread by length, and individuals at by age and time, functions: vbgr.fixed(), vbgr.sd(), N_age
#####
vbgr.fixed <- function(t,n){ # von Bertalanffy growth rate for herring. 
  #Parameters 2005-2010 for Western Baltic spring spawning herring (Gröhsler et al., 2013)
  Linf <- 30.57 # cm
  k <- 0.453 # vbgr growth rate
  hatching <- c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365
  t_hatch <- mean(hatching) # From figure 4 in Polte et al., 2021
  #t_hatch <- rnorm(n,mean=mean(hatching)*365,sd=7.5) # From figure 4 in Polte et al., 2021. Mean is from 4A sd is based on 4B
  L_hatch <- 0.645 # figure 2 in Bauer et al., 2014
  #L_hatch <- rnorm(n,mean=0.645,sd=0.036)# figure 2 in Bauer et al., 2014
  fish.age <- t-t_hatch
  u <- 0.22 # Amplitude of seasonality (same as cod)
  w <- 0.90 # Timing of peak growth
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+fish.age-phi_hatch))+Linf# vbgr estimated length
  return(L_t)
}

vbgr.sd <- function(t,n){ # von Bertalanffy growth rate for herring. 
  #Parameters 2005-2010 for Western Baltic spring spawning herring (Gröhsler et al., 2013)
  
  #Linf <- rnorm(n,mean=30.57,sd=sd(c(31.16,31.81,29.92,29.84,29.34,30.82))) # cm - from table 2 in Gröhsler et al, 2013.
  Linf <- rnorm(n,mean=30.57,sd=2.47254819)
  k.s <- c(0.36,0.365,0.525,0.539,0.552,0.482) # from table 2 in Gröhsler et al, 2013
  #k <- rnorm(n,mean(k.s),sd(k.s)) # vbgr growth rate
  k <- rnorm(n,mean(k.s),0.06009264) # vbgr growth rate
  hatching <- c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365
  #t_hatch <- mean(hatching) # From figure 4 in Polte et al., 2021
  #t_hatch <- rnorm(n,mean=mean(hatching)*365,sd=7.5)/365 # From figure 4 in Polte et al., 2021. Mean is from 4A sd is based on 4B
  t_hatch <- rnorm(n,mean=mean(hatching)*365,sd=sqrt(sd(hatching*365)^2+7.5^2))/365 # From figure 4 in Polte et al., 2021. Mean is from 4A sd is based on 4B
  
  #L_hatch <- 0.645 # figure 2 in Bauer et al., 2014
  L_hatch <- 0.645 # figure 2 in Bauer et al., 2014
  
  sd <- rep(0,length(t))
  for (i in 1:length(sd)){
    L_t <- (L_hatch-Linf)*(exp(-k*(t[i]-t_hatch)))+Linf# vbgr estimated length
    sd[i] <- sd(L_t)
  }
  return(sd)
}


setwd("/Users/mavast/Desktop/Markus SOS")
sam.herring <- read.table("herring in sub.div. 22_24 (SAM).csv",header=T,sep=';')

H <- data.frame(years = rep(sam.herring$Year,9),
                Number = c(sam.herring$age0,sam.herring$age1,
                           sam.herring$age2,sam.herring$age3,
                           sam.herring$age4,sam.herring$age5,
                           sam.herring$age6,sam.herring$age7,
                           sam.herring$age8),
                parameter = rep(sam.herring$Parameter,9),
                ages = rep(0:8,each=length(sam.herring$Year)))
H_age <- cbind(H[which(H$parameter=="N"),],H$Number[which(H$parameter=="F")])
names(H_age)[names(H_age)=="H$Number[which(H$parameter == \"F\")]"] <- "Fishing"
H_age$Z <- rep(0,length(H_age$years))
for (i in 1:length(H_age$years)){
  dada.i <- H_age[i,]
  if (dada.i$years<max(H_age$years) & dada.i$ages<max(H_age$ages)){
    dada.ip1 <- H_age %>% filter(years==(dada.i$years+1) & ages==(dada.i$ages+1))
    H_age$Z[i] <- -log(dada.ip1$Number/dada.i$Number)
  }
  if (dada.i$years==max(H_age$years) | dada.i$ages==max(H_age$ages)){
    H_age$Z[i] <- dada.i$Fishing
  }
  
}



N_age <- function(age,yrs){ # Number of individuals by age - Model N is from the start of the year
  Ns <- matrix(rep(0,(length(age)+1)*length(yrs)*365),ncol = length(age)+1)
  Ns[,length(age)+1] <- rep(yrs,each=365)
  x <- seq(1/365,1,length=365)
  for(i in 1:length(yrs)){
    for(j in 1:length(age)){
      di <- H_age %>% filter(ages==age[j] & years==yrs[i])
      dip1 <- H_age %>% filter(ages==(age[j]+1) & years==(yrs[i]+1))
      
      Ns[(365*(i-1)+1):(365*i),j] <- di$N*exp(-x*di$Z)
    }
  }
  return(Ns)
}
#rm(list=setdiff(ls(),c('vbgr.fixed','vbgr.sd','N_age','H_age')))
#####

#read and prepare cormorant and herring data
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

herring_corm <- subset(d,ART=='herring')
#rm(list=setdiff(ls(),c('herring_corm','H_age','vbgr.fixed','vbgr.sd','N_age')))
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
herring_seal <- herring[,c('Scat','Site','Year','Month','size')]
#rm(list=setdiff(ls(),c('herring_seal','herring_corm','H_age','vbgr.fixed','vbgr.sd','N_age','cod_seal')))
#####

# Seal preference for herring
#####
herring_seal$my <- paste(herring_seal$Month,'-',herring_seal$Year)
herring_seal$size.agg <- floor(herring_seal$size/100*2)/2*100

samplings <- unique(herring_seal$my)
month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
           "Sep","Oct","Nov","Dec") 

lng.class <- data.frame(l.class = 50*(0:10))

df.seal <- data.frame(l.class = rep(50*(0:10),length(samplings)),
                      prop_eaten = rep(0,11*length(samplings)),
                      eaten = rep(0,11*length(samplings)),
                      not_eaten = rep(0,11*length(samplings)),
                      n_scat = rep(0,11*length(samplings)),
                      year = rep(0,11*length(samplings)),
                      location = rep(0,11*length(samplings)),
                      sampling = rep(0,11*length(samplings)))


for (i in 1:length(samplings)){
  dada <- herring_seal %>% filter(my==samplings[i])
  n_sca <- length(unique(dada$Scat))
  date <- paste(dada$Year[1],'-',which(dada$Month[1]==month)[[1]],'-',15,sep="")
  jday <- yday(date)
  sizes <- c(vbgr.fixed(jday/365),vbgr.fixed(jday/365+1),vbgr.fixed(jday/365+2),vbgr.fixed(jday/365+3),
             vbgr.fixed(jday/365+4),vbgr.fixed(jday/365+5),vbgr.fixed(jday/365+6),vbgr.fixed(jday/365+7),vbgr.fixed(jday/365+7))*10
  ages <- c(0,1,2,3,4,5,6,7,8)
  pop <- N_age(ages,dada$Year[1])[jday,1:length(ages)]
  sd<- c(vbgr.sd(jday/365,100000),vbgr.sd(jday/365+1,100000),
         vbgr.sd(jday/365+2,100000),vbgr.sd(jday/365+3,100000),
         vbgr.sd(jday/365+4,100000),vbgr.sd(jday/365+5,100000),
         vbgr.sd(jday/365+6,100000),vbgr.sd(jday/365+7,100000),vbgr.sd(jday/365+8,100000))*10
  
  if (min(sd)<0){
    sd[sd<0] <-rep(0,length(sd[sd<0]))
  }

  available <- c(floor(rnorm(pop[1],sizes[1],sd[1])/10)*10,floor(rnorm(pop[2],sizes[2],sd[2])/10)*10,
                 floor(rnorm(pop[3],sizes[3],sd[3])/10)*10,floor(rnorm(pop[4],sizes[4],sd[4])/10)*10,
                 floor(rnorm(pop[5],sizes[5],sd[5])/10)*10,floor(rnorm(pop[6],sizes[6],sd[6])/10)*10,
                 floor(rnorm(pop[7],sizes[7],sd[7])/10)*10,floor(rnorm(pop[8],sizes[8],sd[8])/10)*10,
                 floor(rnorm(pop[8],sizes[8],sd[8])/10)*10)
  available[available<0] <- rep(0,length(available[available<0]))
  
  avail <- as.data.frame(table(floor(available/100*2)/2*100))
  colnames(avail) <- c("l.class", "available")
  eat <- as.data.frame(table(dada$size.agg))
  colnames(eat) <- c("l.class", "eaten")
  eat$l.class <- as.numeric(as.character(eat$l.class))
  eat$l.class[eat$l.class>300] <- rep(300,length(eat$l.class[eat$l.class>300]))
  
  eat <-  aggregate(eaten~l.class,data=eat,FUN=sum)
  
  df <- merge(lng.class,avail,by='l.class',all.x = T)
  df <- merge(df,eat,by='l.class',all.x = T)
  df[is.na(df)] <- rep(0,length(df[is.na(df)]))
  
  df$eaten <- df$eaten
  df$not_eaten <- df$available-df$eaten
  df$prop_eaten <- (df$eaten/df$available)
  df$n_scat <- rep(n_sca,11)
  df$location <- rep(dada$Site[1],11)
  df$year <- rep(dada$Year[1],11)
  df$sampling <- rep(samplings[i],11)
  
  df.seal[((i-1)*11+1):(i*11),c('prop_eaten','eaten','not_eaten','n_scat','year','location','sampling')] <- 
    df[,c('prop_eaten','eaten','not_eaten','n_scat','year','location','sampling')]
  print(i)
}
df.seal.ori <- df.seal
df.seal <- df.seal.ori %>% filter(!is.na(prop_eaten))
df.seal <- df.seal %>% filter(abs(prop_eaten)<1)
# 2015
df.15 <- df.seal %>% filter(year==2015) 
fac.15 <- aggregate(n_scat~l.class,data=df.15,FUN=sum)
mean(df.15$prop_eaten)
df.15$weighted.prop <- df.15$prop_eaten*df.15$n_scat
df.fin.15 <- aggregate(weighted.prop~l.class,data=df.15,FUN=sum)
df.fin.15$weighted.prop <- df.fin.15$weighted.prop/fac.15$n_scat
mean(df.fin.15$weighted.prop)

plot(df.fin.15$l.class,df.fin.15$weighted.prop,xlim=c(0,500),ylim=c(0,10^-4),
     pch=19,xlab = 'cm',ylab = 'preference',col='blue')
# 2016
df.16 <- df.seal %>% filter(year==2016) 
fac.16 <- aggregate(n_scat~l.class,data=df.16,FUN=sum)
df.16$weighted.prop <- df.16$prop_eaten*df.16$n_scat
df.fin.16 <- aggregate(weighted.prop~l.class,data=df.16,FUN=sum)
df.fin.16$weighted.prop <- df.fin.16$weighted.prop/fac.16$n_scat

points(df.fin.16$l.class,df.fin.16$weighted.prop,pch=19,col='red')
# 2017
df.17 <- df.seal %>% filter(year==2017) 
fac.17 <- aggregate(n_scat~l.class,data=df.17,FUN=sum)
df.17$weighted.prop <- df.17$prop_eaten*df.17$n_scat
df.fin.17 <- aggregate(weighted.prop~l.class,data=df.17,FUN=sum)
df.fin.17$weighted.prop <- df.fin.17$weighted.prop/fac.17$n_scat

points(df.fin.17$l.class,df.fin.17$weighted.prop,pch=19,col='green')
# All together - df.fin
fac <- aggregate(n_scat~l.class,data=df.seal,FUN=sum)
df.seal$weighted.prop <- df.seal$prop_eaten*df.seal$n_scat
df.seal.fin <- aggregate(weighted.prop~l.class,data=df.seal,FUN=sum)
df.seal.fin$weighted.prop <- df.seal.fin$weighted.prop/fac$n_scat
points(df.seal.fin$l.class,df.seal.fin$weighted.prop,pch=19)
legend(0,10^-4,legend = c('2015','2016','2017','All'),pch=19,
       col=c('blue','red','green','black'))

#####

#seal fit cod - run multiple times until lowest AIC (~9.1). NB improve this method
#####
dataInput <- df.seal.fin 
dataInput$weighted.prop[dataInput$l.class>300] <- dataInput$weighted.prop[dataInput$l.class==300] # NB data pruning should be removed when more data available
dataInput <- dataInput %>% filter(dataInput$l.class!=250)
names(dataInput)[names(dataInput)=="l.class"] <- "time"
names(dataInput)[names(dataInput)=="weighted.prop"] <- "intensity"
time <- dataInput$time
normalizedInput <- normalizeData(dataInput)
for (i in 1:1000){
  parameterVector <- sigmoidalFitFunction(normalizedInput,
                                          tryCounter = 200,n_iterations = 1024)
  if (parameterVector$isThisaFit==T & parameterVector$AIC_value<(-43)){
    break
  } 
} 
#parameterVector <- sigmoidalFitFunction(normalizedInput,
#                                        tryCounter = 200,n_iterations = 1024)
#parameterVector$AIC_value
#Check the results
if(parameterVector$isThisaFit){
  time <- 0:400
  intensityTheoretical <- sigmoidalFitFormula(time,
                                              maximum = parameterVector$maximum_Estimate,
                                              slopeParam = parameterVector$slopeParam_Estimate,
                                              midPoint = parameterVector$midPoint_Estimate)
  comparisonData <- cbind(time, intensityTheoretical)
  require(ggplot2)
  ggplot(dataInput) +
    geom_point(aes(x = time, y = intensity)) +
    geom_line(data = comparisonData, aes(x = time, y = intensityTheoretical)) +
    expand_limits(x = 0, y = 0)
}

seal_fit <- cbind(time, intensityTheoretical)

#rm(list=setdiff(ls(),c('herring_corm','H_age','vbgr.fixed','vbgr.sd','N_age','herring_seal','df.seal.fin','seal_fit')))
#####

# Cormorant preference for herring. NB takes a while
#####
herring_corm$dmy <- paste(herring_corm$year,'-',herring_corm$month,'-',herring_corm$day,sep="")


samplings <- unique(herring_corm$dmy)
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
  dada <- herring_corm %>% filter(dmy==samplings[i])
  n_pel <- length(unique(dada$GYLP))
  date <- paste(dada$year[1],'-',which(dada$month[1]==month)[[1]],'-',dada$day[1],sep="")
  jday <- yday(date)
  sizes <- c(vbgr.fixed(jday/365),vbgr.fixed(jday/365+1),vbgr.fixed(jday/365+2),vbgr.fixed(jday/365+3),
             vbgr.fixed(jday/365+4),vbgr.fixed(jday/365+5),vbgr.fixed(jday/365+6),vbgr.fixed(jday/365+7),
             vbgr.fixed(jday/365+8))*10
  ages <- c(0,1,2,3,4,5,6,7,8)
  pop <- N_age(ages,dada$year[1])[jday,1:length(ages)]
  sd<- c(vbgr.sd(jday/365,100000),vbgr.sd(jday/365+1,100000),vbgr.sd(jday/365+2,100000),vbgr.sd(jday/365+3,100000),
         vbgr.sd(jday/365+4,100000),vbgr.sd(jday/365+5,100000),vbgr.sd(jday/365+6,100000),vbgr.sd(jday/365+7,100000),
         vbgr.sd(jday/365+8,100000))*10
  
  if (min(sd)<0){
    sd[sd<0] <-rep(0,length(sd[sd<0]))
  }
  available <- c(floor(rnorm(pop[1],sizes[1],sd[1])/10)*10,floor(rnorm(pop[2],sizes[2],sd[2])/10)*10,
                 floor(rnorm(pop[3],sizes[3],sd[3])/10)*10,floor(rnorm(pop[4],sizes[4],sd[4])/10)*10,
                 floor(rnorm(pop[5],sizes[5],sd[5])/10)*10,floor(rnorm(pop[6],sizes[6],sd[6])/10)*10,
                 floor(rnorm(pop[7],sizes[7],sd[7])/10)*10,floor(rnorm(pop[8],sizes[8],sd[8])/10)*10,
                 floor(rnorm(pop[9],sizes[9],sd[9])/10)*10)
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
  print(paste(round(i/44*100,1),'% done'))
}


df.corm.ori <- df.corm
df.corm <- df.corm.ori %>% filter(!is.na(prop_eaten))
df.corm <- df.corm %>% filter(abs(prop_eaten)<1)
# 1992
df.92 <- df.corm %>% filter(year==1992) 
fac.92 <- aggregate(n_pellet~l.class,data=df.92,FUN=sum)
mean(df.92$prop_eaten)
df.92$weighted.prop <- df.92$prop_eaten*df.92$n_pellet
df.fin.92 <- aggregate(weighted.prop~l.class,data=df.92,FUN=sum)
df.fin.92$weighted.prop <- df.fin.92$weighted.prop/fac.92$n_pellet
mean(df.fin.92$weighted.prop)

plot(df.fin.92$l.class,df.fin.92$weighted.prop,xlim=c(0,500),
     pch=19,xlab = 'cm',ylab = 'preference',col='blue')
# 1993
df.93 <- df.corm %>% filter(year==1993) 
fac.93 <- aggregate(n_pellet~l.class,data=df.93,FUN=sum)
df.93$weighted.prop <- df.93$prop_eaten*df.93$n_pellet
df.fin.93 <- aggregate(weighted.prop~l.class,data=df.93,FUN=sum)
df.fin.93$weighted.prop <- df.fin.93$weighted.prop/fac.93$n_pellet

points(df.fin.93$l.class,df.fin.93$weighted.prop,pch=19,col='red')
# 1994
df.94 <- df.corm %>% filter(year==1994) 
fac.94 <- aggregate(n_pellet~l.class,data=df.94,FUN=sum)
df.94$weighted.prop <- df.94$prop_eaten*df.94$n_pellet
df.fin.94 <- aggregate(weighted.prop~l.class,data=df.94,FUN=sum)
df.fin.94$weighted.prop <- df.fin.94$weighted.prop/fac.94$n_pellet

points(df.fin.94$l.class,df.fin.94$weighted.prop,pch=19,col='green')

# All together - df.fin
fac <- aggregate(n_pellet~l.class,data=df.corm,FUN=sum)
df.corm$weighted.prop <- df.corm$prop_eaten*df.corm$n_pellet
df.corm.fin <- aggregate(weighted.prop~l.class,data=df.corm,FUN=sum)
df.corm.fin$weighted.prop <- df.corm.fin$weighted.prop/fac$n_pellet
points(df.corm.fin$l.class,df.corm.fin$weighted.prop,pch=19)
legend(0,0.25,legend = c('1992','1993','1994','All'),pch=19,
       col=c('blue','red','green','black'))


#####

#cormorant fit cod - run multiple times until lowest AIC (~-49). NB improve this method
#####
dataInput <- df.corm.fin %>% filter(l.class<400)
names(dataInput)[names(dataInput)=="l.class"] <- "time"
names(dataInput)[names(dataInput)=="weighted.prop"] <- "intensity"
time <- dataInput$time
normalizedInput <- normalizeData(dataInput)
for (i in 1:1000){
  parameterVector <- doublesigmoidalFitFunction(normalizedInput,
                                                tryCounter = 200,n_iterations = 1024)
  if (parameterVector$isThisaFit==T & parameterVector$AIC_value<(-45)){
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
#rm(list=setdiff(ls(),c('herring_corm','H_age','vbgr.fixed',
#                       'vbgr.sd','N_age','herring_seal','df.seal.fin',
#                       'seal_fit','cormorant_fit','df.corm.fin')))

#####

# Seal and cormorant cod preference
#double sigmoidal fit
#####
par(mfrow=c(2,1))
par(mar = c(2.5,4,1,1))
plot(df.corm.fin$l.class,df.corm.fin$weighted.prop,pch=19,
     main = 'Cormorants',xlim=c(0,400),
     ylab='preference',xlab = '', xaxt='n')
lines(cormorant_fit[,1],cormorant_fit[,2]/sum(cormorant_fit[,2]),
      type='l',col='grey60',lwd=3)
par(mar = c(4,4,1,1))


plot(df.seal.fin$l.class,df.seal.fin$weighted.prop/max(df.seal.fin$weighted.prop),pch=19,
     main = 'Seals',xlim=c(0,400),
     ylab='preference',xlab = 'cod length [mm]')
lines(seal_fit[,1],seal_fit[,2]/max(seal_fit[,2]),
      type='l',col='grey60',lwd=3)

par(mfrow=c(1,1))
par(mar = c(4.5,4,1,2))
plot(seal_fit[,1],seal_fit[,2]/mean(seal_fit[,2])/max(seal_fit[,2]/mean(seal_fit[,2]))
     ,pch=19,col='orange',
     type='l',lwd=2,xlab = 'cod length [mm]',ylab = 'Relative preference')
lines(cormorant_fit[,1],cormorant_fit[,2]/mean(cormorant_fit[,2])/max(cormorant_fit[,2]/mean(cormorant_fit[,2])),
      type='l',col='red',lwd=2)
legend(0,1,legend=c('Seal','Cormorant'),
       pch=19,col=c('orange','red'))
#####
