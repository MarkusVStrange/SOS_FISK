source('estimate uncertainty on coefficients.R')

data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="")
# Define cod functions and read data
#####
# Growth model, cohort spread by length, and individuals at by age and time, functions: vbgr.fixed(), vbgr.sd(), N_age

vbgrCod <- function(age){ # von Bertalanffy growth rate. Parameters from McQueen et al., 2019 with fixed hatching time
  # define growth parameters
  Linf <- 97.9 # cm. McQueen et al., 2019
  L_hatch <- 0.4 # length at hatching, cm (Pepin et al, 1997)
  k <- 0.22 # vbgr growth rate, from McQueen
  u <- 0.23 # Amplitude of seasonality
  w <- 0.91 # Timing of peak growth
  #t_hatch <- ((46+135)/2)/365 # time of hatching, Julian day / 365 - spawning from Feb. to May (Fiskeatlas), average set to April 1st
  t_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf # vbgr estimated length
  return(L_t)
}


vbgr.sdCod <- function(age,n){ # von Bertalanffy growth rate. Parameters from McQueen et al., 2019 and estimated
  # t in fraction year, Julian day/365
  Linf <- rnorm(n,mean=97.9,sd=cod_sd_coef[3]) # Assymptotic length - with stochasticity
  #Linf <- 154.56 # cm
  L_hatch <- 0.4 # length at hatching, cm (Pepin et al, 1997)
  #t_hatch <- ((46+135)/2)/365
  t_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
  k <- rnorm(n,mean=0.22,sd=cod_sd_coef[1]) # vbgr growth rate
  u <- 0.23 # Amplitude of seasonality
  w <- 0.91 # Timing of peak growth
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  sd <- rep(0,length(t))
  for (i in 1:length(sd)){
    fish.age <- age[i]+rnorm(n,mean=0,sd=cod_sd_coef[2]) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
    L_t <- (L_hatch-Linf)*exp(-k*(phi_t[i]+fish.age-phi_hatch))+Linf # vbgr estimated length
    sd[i] <- sd(L_t)
  }
  return(sd)
}

N_ageCod <- function(dat,age,yrs){ # Number of individuals by age
  dat$Z <- dat$M_var+dat$Fest
  Ns <- matrix(rep(0,(length(age)+1)*length(yrs)*365),ncol = length(age)+1)
  Ns[,length(age)+1] <- rep(yrs,each=365)
  x <- seq(1/365,1,length=365)
  for(i in 1:length(yrs)){
    for(j in 1:length(age)){
      di <- dat %>% filter(ages==age[j] & years==yrs[i])
      dip1 <- dat %>% filter(ages==(age[j]+1) & years==(yrs[i]+1))
      
      Ns[(365*(i-1)+1):(365*i),j] <- di$N*exp(-x*di$Z)
    }
  }
  return(Ns)
}

Cod <- read.table(paste(data_wd,'M_sms.csv',sep=""),sep=';',header=TRUE)
#Cod$N <- readRDS('data/sms_mortality.RDS')$N$N
cormorant <- read.table(paste(data_wd,'corm_diet.csv',sep=""),sep=';',header=TRUE)
gseal <- read.table(paste(data_wd,'gSeal_diet.csv',sep=""),sep=';',header=TRUE)
gseal <- gseal %>% filter(Site!="Utklippan")

#####

# Define herring functions and read data
#####
vbgrHerring <- function(age){ # von Bertalanffy growth rate for herring. 
  #Parameters 2005-2010 for Western Baltic spring spawning herring (Gröhsler et al., 2013)
  Linf <- 30.57 # cm
  k <- 0.453 # from table 2 in Gröhsler et al, 2013
  hatching <- c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365
  t_hatch <- mean(hatching) # From figure 4 in Polte et al., 2021
  #t_hatch <- rnorm(n,mean=mean(hatching)*365,sd=7.5) # From figure 4 in Polte et al., 2021. Mean is from 4A sd is based on 4B
  L_hatch <- 0.645 # figure 2 in Bauer et al., 2014
  u <- 0.23 # Amplitude of seasonality
  w <- 0.91 # Timing of peak growth
  #t_hatch <- ((46+135)/2)/365 # time of hatching, Julian day / 365 - spawning from Feb. to May (Fiskeatlas), average set to April 1st
  t_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf# vbgr estimated length
  return(L_t)
}

vbgr.sdHerring <- function(age,n){ # von Bertalanffy growth rate for herring. 
  #Parameters 2005-2010 for Western Baltic spring spawning herring (Gröhsler et al., 2013)
  Linf <- rnorm(n,mean=30.57,sd=herring_sd_coef[3])
  k <- rnorm(n,0.453,herring_sd_coef[1]) # from table 2 in Gröhsler et al, 2013
  L_hatch <- 0.645 # figure 2 in Bauer et al., 2014
  hatching <- c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365
  t_hatch <- mean(hatching) # From figure 4 in Polte et al., 2021
  
  u <- 0.23 # Amplitude of seasonality
  w <- 0.91 # Timing of peak growth
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  sd <- rep(0,length(t))
  for (i in 1:length(sd)){
    fish.age <- age[i]+rnorm(n,mean=0,sd=herring_sd_coef[2]) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
    L_t <- (L_hatch-Linf)*exp(-k*(phi_t[i]+fish.age-phi_hatch))+Linf # vbgr estimated length
    sd[i] <- sd(L_t)
  }
  return(sd)
}

N_ageHerring <- function(dat,age,yrs){ # Number of individuals by age - Model N is from the start of the year
  Ns <- matrix(rep(0,(length(age)+1)*length(yrs)*365),ncol = length(age)+1)
  Ns[,length(age)+1] <- rep(yrs,each=365)
  x <- seq(1/365,1,length=365)
  for(i in 1:length(yrs)){
    for(j in 1:length(age)){
      di <- Herring %>% filter(ages==age[j] & years==yrs[i])
      dip1 <- Herring %>% filter(ages==(age[j]+1) & years==(yrs[i]+1))
      
      Ns[(365*(i-1)+1):(365*i),j] <- di$N*exp(-x*di$Z)
    }
  }
  return(Ns)
}
Herring <- read.table(paste(data_wd,'Herring.csv',sep=""),sep=';',header=TRUE)
#####

# Define Flatfish functions and read data
#####
# define flounder growth function
vbgrFlounder <- function(age) {
  Linf <- 34.0031653 # estimated in "fit flatfish growths from DATRAS" 
  k <- 0.4537312  # estimated in "fit flatfish growths from DATRAS" 
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  u <- 0.23 # Amplitude of seasonality (same as cod)
  w <- 0.91 # Timing of peak growth
  t_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
  
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
  
  return(L_t)
} #Original

vbgr.sdFlounder <- function(age,n){
  #Parameters 2005-2010 for Western Baltic spring spawning herring (Gröhsler et al., 2013)
  Linf <- rnorm(n,mean=34.0031653,sd=flounder_sd_coef[3])
  k <- rnorm(n,0.4537312,flounder_sd_coef[1]) # from table 2 in Gröhsler et al, 2013
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  u <- 0.23 # Amplitude of seasonality (same as cod)
  w <- 0.91 # Timing of peak growth
  t_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
  
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  sd <- rep(0,length(t))
  for (i in 1:length(sd)){
    fish.age <- age[i]+rnorm(n,mean=0,sd=flounder_sd_coef[2]) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
    L_t <- (L_hatch-Linf)*exp(-k*(phi_t[i]+fish.age-phi_hatch))+Linf # vbgr estimated length
    sd[i] <- sd(L_t)
  }
  return(sd)
}


# define plaice growth function
vbgrPlaice <- function(age) {
  Linf <- 41.637260  # estimated in "fit flatfish growths from DATRAS" 
  k <- 0.240703   # estimated in "fit flatfish growths from DATRAS" 
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  u <- 0.23 # Amplitude of seasonality (same as cod)
  w <- 0.91 # Timing of peak growth
  t_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
  
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
  
  return(L_t)
} #Original

vbgr.sdPlaice <- function(age,n){
  #Parameters 2005-2010 for Western Baltic spring spawning herring (Gröhsler et al., 2013)
  Linf <- rnorm(n,mean=41.637260,sd=plaice_sd_coef[3])
  k <- rnorm(n,0.240703,plaice_sd_coef[1]) # from table 2 in Gröhsler et al, 2013
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  u <- 0.23 # Amplitude of seasonality (same as cod)
  w <- 0.91 # Timing of peak growth
  t_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
  
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  sd <- rep(0,length(t))
  for (i in 1:length(sd)){
    fish.age <- age[i]+rnorm(n,mean=0,sd=plaice_sd_coef[2]) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
    L_t <- (L_hatch-Linf)*exp(-k*(phi_t[i]+fish.age-phi_hatch))+Linf # vbgr estimated length
    sd[i] <- sd(L_t)
  }
  return(sd)
}

# define dab growth function
vbgrDab <- function(age) {
  Linf <- 26.5093927   # estimated in "fit flatfish growths from DATRAS" 
  k <- 0.4846326    # estimated in "fit flatfish growths from DATRAS" 
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  u <- 0.23 # Amplitude of seasonality (same as cod)
  w <- 0.91 # Timing of peak growth
  t_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
  
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
  
  return(L_t)
} #Original

vbgr.sdDab <- function(age,n){
  #Parameters 2005-2010 for Western Baltic spring spawning herring (Gröhsler et al., 2013)
  Linf <- rnorm(n,mean=26.5093927,sd=dab_sd_coef[3]  )
  k <- rnorm(n,0.4846326,dab_sd_coef[1] ) # from table 2 in Gröhsler et al, 2013
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  u <- 0.23 # Amplitude of seasonality (same as cod)
  w <- 0.91 # Timing of peak growth
  t_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching

  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  sd <- rep(0,length(t))
  for (i in 1:length(sd)){
    fish.age <- age[i]+rnorm(n,mean=0,sd=dab_sd_coef[2] ) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
    L_t <- (L_hatch-Linf)*exp(-k*(phi_t[i]+fish.age-phi_hatch))+Linf # vbgr estimated length
    sd[i] <- sd(L_t)
  }
  return(sd)
}

Flatfish <- read.table(paste(data_wd,"N_est.flatfish.csv",sep=""),sep=';',header=TRUE)
Flatfish$yday <- round(Flatfish$jd*365)
names(Flatfish)[names(Flatfish)=="Species"] <- "species"
#####


# TIP: Run cod together
# Seal preference for cod
#####
gseal.cod <- gseal %>% filter(Species %in% c("cod"))
gseal.cod$my <- paste(gseal.cod$Month,'-',gseal.cod$Year)
samplings <- unique(gseal.cod$my)
month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
           "Sep","Oct","Nov","Dec") 
width <- 10 # mm. width of the size classes
vec.lengths <- seq(width/2,1500+width/2,by=width)
n.lengths <- length(vec.lengths)
gseal.cod$size <- floor(gseal.cod$FL.with.SCF/width)*width+width/2


df.seal <- data.frame(l.class = rep(vec.lengths,length(samplings)),
                      n_eaten = rep(0,n.lengths*length(samplings)),
                      n =  rep(0,n.lengths*length(samplings)),
                      sampling = rep(samplings,each=n.lengths),
                      n_scat =  rep(0,n.lengths*length(samplings)),
                      year =  rep(0,n.lengths*length(samplings)))

cod_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning

for (i in 1:length(samplings)){
  dada <- gseal.cod %>% filter(my==samplings[i])
  n_sca <- length(unique(dada$Scat))
  date <- paste(dada$Year[1],'-',which(dada$Month[1]==month)[[1]],'-',15,sep="")
  jday <- yday(date)/365
  sizes <- c(vbgrCod(jday-cod_hatch),vbgrCod(jday-cod_hatch+1),vbgrCod(jday-cod_hatch+2),vbgrCod(jday-cod_hatch+3),
             vbgrCod(jday-cod_hatch+4),vbgrCod(jday-cod_hatch+5),vbgrCod(jday-cod_hatch+6),vbgrCod(jday-cod_hatch+7))*10
  ages <- c(0,1,2,3,4,5,6,7)
  cod.i <- N_ageCod(Cod,ages,dada$Year[1])[jday*365,1:length(ages)]
  sd<- c(vbgr.sdCod(jday-cod_hatch,100000),vbgr.sdCod(jday-cod_hatch+1,100000),vbgr.sdCod(jday-cod_hatch+2,100000),vbgr.sdCod(jday-cod_hatch+3,100000),
         vbgr.sdCod(jday-cod_hatch+4,100000),vbgr.sdCod(jday-cod_hatch+5,100000),vbgr.sdCod(jday-cod_hatch+6,100000),vbgr.sdCod(jday-cod_hatch+7,100000))*10
  
  if (min(sd)<0){
    sd[sd<0] <-rep(0,length(sd[sd<0]))
  }
  available <- c(floor(rnorm(cod.i[1],sizes[1],sd[1])/width)*width+width/2,floor(rnorm(cod.i[2],sizes[2],sd[2])/width)*width+width/2,
                 floor(rnorm(cod.i[3],sizes[3],sd[3])/width)*width+width/2,floor(rnorm(cod.i[4],sizes[4],sd[4])/width)*width+width/2,
                 floor(rnorm(cod.i[5],sizes[5],sd[5])/width)*width+width/2,floor(rnorm(cod.i[6],sizes[6],sd[6])/width)*width+width/2,
                 floor(rnorm(cod.i[7],sizes[7],sd[7])/width)*width+width/2,floor(rnorm(cod.i[8],sizes[8],sd[8])/width)*width+width/2)
  available[available<0] <- rep(width/2,length(available[available<0]))

  avail <- as.data.frame(table(available))
  colnames(avail) <- c("l.class", "available")
  eat <- as.data.frame(table(dada$size))
  colnames(eat) <- c("l.class", "eaten")
  
  ava.idx <- unique(df.seal$l.class) %in% avail$l.class
  eat.idx <- unique(df.seal$l.class) %in% eat$l.class
  i.idx <- (((i-1)*n.lengths+1):(i*n.lengths))
  df.seal$n[i.idx[ava.idx]] <- avail$available
  df.seal$n_eaten[i.idx[eat.idx]] <- eat$eaten
  df.seal$n_scat[i.idx] <- rep(length(unique(dada$Scat)),n.lengths)
#  df$location <- rep(dada$Site[1],151)
  df.seal$year[i.idx] <- rep(dada$Year[1],n.lengths)
}
df.seal$odds <- df.seal$n_eaten/df.seal$n
df.seal <- df.seal %>% filter(!is.na(odds))
plot(df.seal$l.class,df.seal$odds)
# Remove outliers - NB put in again when more data
df.seal <- df.seal %>% filter(odds<0.03)

# 2015
df.15 <- df.seal %>% filter(year==2015) 
tot.scat.15 <- sum(aggregate(n_scat~sampling,data=df.15,FUN=mean)$n_scat)
mean(df.15$odds)
df.15$weighted.odds <- (df.15$n_eaten/df.15$n)*df.15$n_scat
df.fin.15 <- aggregate(weighted.odds~l.class,data=df.15,FUN=sum)
df.fin.15$weighted.odds <- df.fin.15$weighted.odds/tot.scat.15
mean(df.fin.15$weighted.odds)
# 2016
df.16 <- df.seal %>% filter(year==2016) 
tot.scat.16 <- sum(aggregate(n_scat~sampling,data=df.16,FUN=mean)$n_scat)
mean(df.16$odds)
df.16$weighted.odds <- (df.16$n_eaten/df.16$n)*df.16$n_scat
df.fin.16 <- aggregate(weighted.odds~l.class,data=df.16,FUN=sum)
df.fin.16$weighted.odds <- df.fin.16$weighted.odds/tot.scat.16
mean(df.fin.16$weighted.odds)
# 2017
df.17 <- df.seal %>% filter(year==2017) 
tot.scat.17 <- sum(aggregate(n_scat~sampling,data=df.17,FUN=mean)$n_scat)
mean(df.17$odds)
df.17$weighted.odds <- (df.17$n_eaten/df.17$n)*df.17$n_scat
df.fin.17 <- aggregate(weighted.odds~l.class,data=df.17,FUN=sum)
df.fin.17$weighted.odds <- df.fin.17$weighted.odds/tot.scat.17
mean(df.fin.17$weighted.odds)
# All together - df.fin
tot.scat <- sum(aggregate(n_scat~sampling,data=df.seal,FUN=mean)$n_scat)
mean(df.seal$odds)
df.seal$weighted.odds <- (df.seal$n_eaten/df.seal$n)*df.seal$n_scat
df.fin.seal <- aggregate(weighted.odds~l.class,data=df.seal,FUN=mean)
df.fin.seal$weighted.odds <- df.fin.seal$weighted.odds*table(df.seal$l.class)[1]/tot.scat
mean(df.fin.seal$weighted.odds)

plot(df.fin.15$l.class,df.fin.15$weighted.odds,xlim=c(0,600),ylim=c(0,0.015),pch=19,
     xlab = 'cm',ylab = 'preference',col='blue')
points(df.fin.16$l.class,df.fin.16$weighted.odds,pch=19,col='red')
points(df.fin.17$l.class,df.fin.17$weighted.odds,pch=19,col='green')
points(df.fin.seal$l.class,df.fin.seal$weighted.odds,pch=19)
legend(450,0.015,legend = c('2015','2016','2017','All'),pch=19,
       col=c('blue','red','green','black'))

scats <- data.frame(l.class = sort(unique(gseal.cod$size)),
                    n = as.numeric(table(gseal.cod$size)))

df.fin.seal <- left_join(df.fin.seal,scats,by="l.class")
df.fin.seal$n[is.na(df.fin.seal$n)] <- 0
df.fin.seal$weight <- df.fin.seal$n+1


df.fin.seal <- df.fin.seal[-which(df.fin.seal$l.class>400 & df.fin.seal$weighted.odds==0),]
double_sigmoid <- function(x, maximum, slope1, midPoint1, slope2, midPointDistance) {
  0 * maximum + 
    (maximum / (1 + exp(-slope1 * (x - midPoint1)))) - 
    (maximum / (1 + exp(-slope2 * (x - (midPoint1 + midPointDistance)))))
}

par <- list(logMax = -5,logS1=-2,logMid=6,logS2=-2,logDis = 4)

#df.fin.seal <- df.fin.seal %>% filter(l.class<700)
#df.fin.seal$l.class <- df.fin.seal$l.class+5
#df.fin.seal <- df.fin.seal[-which(df.fin.seal$l.class<150 & df.fin.seal$weighted.odds>0.0038),]
obj_fn <- function(par) {
  max <- exp(par[["logMax"]])
  s1 <- exp(par[["logS1"]])
  mid <- exp(par[["logMid"]])
  s2 <- exp(par[["logS2"]])
  dis <- exp(par[["logDis"]])
  pred <- double_sigmoid(df.fin.seal$l.class,max,s1,mid,s2,dis)
  sum((df.fin.seal$weighted.odds - pred)^2)
}
# df.fin.seal$weight*
opt <- nlminb(par, objective = obj_fn)

seal_param <- exp(opt$par)

seal_fit <- double_sigmoid((0:7000)/10,seal_param[1],seal_param[2],seal_param[3],seal_param[4],seal_param[5])
plot(df.fin.seal$l.class,df.fin.seal$weighted.odds,ylab = "preference",
     xlab="length [mm]",main="Grey Seal")
lines((0:7000)/10,seal_fit,
      lwd=2,col="red")


make_my_function <- function(seal_param) {
  maximum <- seal_param[1]   
  slope1 <- seal_param[2]    
  midPoint <- seal_param[3]      
  slope2 <- seal_param[4]    
  midPointDistance <- seal_param[5]   
  
  function(x) { # x is length in mm
    pref <- (maximum / (1 + exp(-slope1 * (x - midPoint)))) - 
      (maximum / (1 + exp(-slope2 * (x - (midPoint + midPointDistance)))))
    pref <- pref/max(pref)
    return(pref)
  }
}
seal_cod.pref <- make_my_function(seal_param)

#####
# cormorant preference for cod
#####
cod_corm <- cormorant %>% filter(ART=="cod")
cod_corm$dmy <- paste(cod_corm$year,'-',cod_corm$month,'-',cod_corm$day,sep="")
width <- 30 # mm. width of the size classes
cod_corm$size <- floor(cod_corm$FISKLGD/width)*width+width/2
samplings <- unique(cod_corm$dmy)
month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
           "Sep","Oct","Nov","Dec")
vec.lengths <- seq(width/2,1500+width/2,by=width)
n.lengths <- length(vec.lengths)

df.corm <- data.frame(l.class = rep(vec.lengths,length(samplings)),
                      n_eaten = rep(0,n.lengths*length(samplings)),
                      n =  rep(0,n.lengths*length(samplings)),
                      sampling = rep(samplings,each=n.lengths),
                      n_scat =  rep(0,n.lengths*length(samplings)),
                      year =  rep(0,n.lengths*length(samplings)))

cod_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
start_time <- Sys.time() #NB takes 5 min
for (i in 1:length(samplings)){
  dada <- cod_corm %>% filter(dmy==samplings[i])
  n_sca <- length(unique(dada$GYLP))
  date <- paste(dada$year[1],'-',which(dada$month[1]==month)[[1]],'-',dada$day[1],sep="")
  jday <- yday(date)/365
  sizes <- c(vbgrCod(jday-cod_hatch),vbgrCod(jday-cod_hatch+1),vbgrCod(jday-cod_hatch+2),vbgrCod(jday-cod_hatch+3),
             vbgrCod(jday-cod_hatch+4),vbgrCod(jday-cod_hatch+5),vbgrCod(jday-cod_hatch+6),vbgrCod(jday-cod_hatch+7))*10
  ages <- 0:7
  cod.i <- N_ageCod(Cod,ages,dada$year[1])[jday*365,1:length(ages)]
  sd<- c(vbgr.sdCod(jday-cod_hatch,100000),vbgr.sdCod(jday-cod_hatch+1,100000),vbgr.sdCod(jday-cod_hatch+2,100000),vbgr.sdCod(jday-cod_hatch+3,100000),
         vbgr.sdCod(jday-cod_hatch+4,100000),vbgr.sdCod(jday-cod_hatch+5,100000),vbgr.sdCod(jday-cod_hatch+6,100000),vbgr.sdCod(jday-cod_hatch+7,100000))*10
  
  if (min(sd)<0){
    sd[sd<0] <-rep(0,length(sd[sd<0]))
  }
  available <- c(floor(rnorm(cod.i[1],sizes[1],sd[1])/width)*width+width/2,floor(rnorm(cod.i[2],sizes[2],sd[2])/width)*width+width/2,
                 floor(rnorm(cod.i[3],sizes[3],sd[3])/width)*width+width/2,floor(rnorm(cod.i[4],sizes[4],sd[4])/width)*width+width/2,
                 floor(rnorm(cod.i[5],sizes[5],sd[5])/width)*width+width/2,floor(rnorm(cod.i[6],sizes[6],sd[6])/width)*width+width/2,
                 floor(rnorm(cod.i[7],sizes[7],sd[7])/width)*width+width/2,floor(rnorm(cod.i[8],sizes[8],sd[8])/width)*width+width/2)
  available[available<0] <- rep(NA,length(available[available<0]))
  
  
  avail <- as.data.frame(table(available))
  colnames(avail) <- c("l.class", "available")
  eat <- as.data.frame(table(dada$size))
  colnames(eat) <- c("l.class", "eaten")
  
  ava.idx <- unique(df.corm$l.class) %in% avail$l.class
  eat.idx <- unique(df.corm$l.class) %in% eat$l.class
  i.idx <- (((i-1)*n.lengths+1):(i*n.lengths))
  df.corm$n[i.idx[ava.idx]] <- avail$available
  df.corm$n_eaten[i.idx[eat.idx]] <- eat$eaten
  df.corm$n_scat[i.idx] <- rep(n_sca,n.lengths)
  #  df$location <- rep(dada$Site[1],n.lengths)
  df.corm$year[i.idx] <- rep(dada$year[1],n.lengths)
}
end_time <- Sys.time()
elapsed <- end_time - start_time
print(elapsed)
df.corm$odds <- df.corm$n_eaten/df.corm$n
df.corm <- df.corm %>% filter(!is.na(odds))
plot(df.corm$l.class,df.corm$odds)
# Remove outliers - NB put in again when more data
#df.corm <- df.corm %>% filter(odds<0.0065)
# 2015
df.15 <- df.corm %>% filter(year==1992) 
tot.scat.15 <- sum(aggregate(n_scat~sampling,data=df.15,FUN=mean)$n_scat)
mean(df.15$odds)
df.15$weighted.odds <- (df.15$n_eaten/df.15$n)*df.15$n_scat
df.fin.15 <- aggregate(weighted.odds~l.class,data=df.15,FUN=sum)
df.fin.15$weighted.odds <- df.fin.15$weighted.odds/tot.scat.15
mean(df.fin.15$weighted.odds)
# 2016
df.16 <- df.corm %>% filter(year==1993) 
tot.scat.16 <- sum(aggregate(n_scat~sampling,data=df.16,FUN=mean)$n_scat)
mean(df.16$odds)
df.16$weighted.odds <- (df.16$n_eaten/df.16$n)*df.16$n_scat
df.fin.16 <- aggregate(weighted.odds~l.class,data=df.16,FUN=sum)
df.fin.16$weighted.odds <- df.fin.16$weighted.odds/tot.scat.16
mean(df.fin.16$weighted.odds)
# 2017
df.17 <- df.corm %>% filter(year==1994) 
tot.scat.17 <- sum(aggregate(n_scat~sampling,data=df.17,FUN=mean)$n_scat)
mean(df.17$odds)
df.17$weighted.odds <- (df.17$n_eaten/df.17$n)*df.17$n_scat
df.fin.17 <- aggregate(weighted.odds~l.class,data=df.17,FUN=sum)
df.fin.17$weighted.odds <- df.fin.17$weighted.odds/tot.scat.17
mean(df.fin.17$weighted.odds)
# All together - df.fin
tot.scat <- sum(aggregate(n_scat~sampling,data=df.corm,FUN=mean)$n_scat)
mean(df.corm$odds)
df.corm$weighted.odds <- (df.corm$n_eaten/df.corm$n)*df.corm$n_scat
df.fin.corm <- aggregate(weighted.odds~l.class,data=df.corm,FUN=mean)
df.fin.corm$weighted.odds <- df.fin.corm$weighted.odds*table(df.corm$l.class)[1]/tot.scat
mean(df.fin.corm$weighted.odds)

plot(df.fin.15$l.class,df.fin.15$weighted.odds,xlim=c(0,500),
     ylim=c(0,max(c(df.fin.15$weighted.odds,df.fin.16$weighted.odds,df.fin.17$weighted.odds))),pch=19,
     xlab = 'cm',ylab = 'preference',col='blue')
points(df.fin.16$l.class,df.fin.16$weighted.odds,pch=19,col='red')
points(df.fin.17$l.class,df.fin.17$weighted.odds,pch=19,col='green')
points(df.fin.corm$l.class,df.fin.corm$weighted.odds,pch=19)
legend(350,max(c(df.fin.15$weighted.odds,df.fin.16$weighted.odds,df.fin.17$weighted.odds))
       ,legend = c('2015','2016','2017','All'),pch=19,
       col=c('blue','red','green','black'))

scats <- data.frame(l.class = sort(unique(cormorant$size)),
                    n = as.numeric(table(cormorant$size)))

df.fin.corm <- left_join(df.fin.corm,scats,by="l.class")
df.fin.corm$n[is.na(df.fin.corm$n)] <- 0
df.fin.corm$weight <- df.fin.corm$n+1


df.fin.corm <- df.fin.corm %>% filter(weighted.odds<0.0014)
double_sigmoid <- function(x, maximum, slope1, midPoint1, slope2, midPointDistance) {
  (maximum / (1 + exp(-slope1 * (x - midPoint1)))) - 
    (maximum / (1 + exp(-slope2 * (x - (midPoint1 + midPointDistance)))))
}

par <- list(logMax = -5,logS1=-2,logMid=6,logS2=-2,logDis = 4)


obj_fn <- function(par) {
  max <- exp(par[["logMax"]])
  s1 <- exp(par[["logS1"]])
  mid <- exp(par[["logMid"]])
  s2 <- exp(par[["logS2"]])
  dis <- exp(par[["logDis"]])
  pred <- double_sigmoid(df.fin.corm$l.class,max,s1,mid,s2,dis)
  sum((df.fin.corm$weighted.odds - pred)^2)
}
# df.fin.seal$weight*
opt <- nlminb(par, objective = obj_fn)
corm_param <- exp(opt$par)

corm_fit <- double_sigmoid((0:7000)/10,corm_param[1],corm_param[2],corm_param[3],corm_param[4],corm_param[5])
plot(df.fin.corm$l.class,df.fin.corm$weighted.odds,ylab = "preference",
     xlab="length [mm]",main="Cormorant",xlim=c(0,500))
lines((0:7000)/10,corm_fit,
      lwd=2,col="orange")




plot((0:7000)/10,corm_fit/max(corm_fit),
     lwd=3,col="darkblue",ylab = "preference",
     xlab="length [mm]",main="Cod",type='l',ylim=c(0,1.53))
lines((0:7000)/10,seal_fit/max(seal_fit),
      lwd=3,col='red')
points(df.fin.corm$l.class,df.fin.corm$weighted.odds/max(corm_fit),pc=19,cex=1,col="darkblue")
points(df.fin.seal$l.class,df.fin.seal$weighted.odds/max(seal_fit),pc=19,cex=1,col="red")
legend(500,1,legend=c('Seal','Cormorant'),
       pch=19,col=c('red','darkblue'))

make_my_function <- function(corm_param) {
  maximum <- corm_param[1]   
  slope1 <- corm_param[2]    
  midPoint <- corm_param[3]      
  slope2 <- corm_param[4]    
  midPointDistance <- corm_param[5]   
  
  function(x) { # x is length in mm
    pref <- (maximum / (1 + exp(-slope1 * (x - midPoint)))) - 
      (maximum / (1 + exp(-slope2 * (x - (midPoint + midPointDistance)))))
    pref <- pref/max(pref)
    return(pref)
  }
}

# Create the function
corm_cod.pref <- make_my_function(corm_param)

#####

# TIP: Run herring together
# Seal preference for herring
#####
gseal.herring <- gseal %>% filter(Species %in% c("herring"))
gseal.herring$my <- paste(gseal.herring$Month,'-',gseal.herring$Year)
width <- 20 # mm. width of the size classes
gseal.herring$size <- floor(gseal.herring$FL.with.SCF/width)*width+width/2
samplings <- unique(gseal.herring$my)
month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
           "Sep","Oct","Nov","Dec")
vec.lengths <- seq(width/2,500+width/2,by=width)
n.lengths <- length(vec.lengths)
df.seal <- data.frame(l.class = rep(vec.lengths,length(samplings)),
                      n_eaten = rep(0,n.lengths*length(samplings)),
                      n =  rep(0,n.lengths*length(samplings)),
                      sampling = rep(samplings,each=n.lengths),
                      n_scat =  rep(0,n.lengths*length(samplings)),
                      year =  rep(0,n.lengths*length(samplings)))

herring_hatch <- mean(c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365)

for (i in 1:length(samplings)){
  dada <- gseal.herring %>% filter(my==samplings[i])
  n_sca <- length(unique(dada$Scat))
  date <- paste(dada$Year[1],'-',which(dada$Month[1]==month)[[1]],'-',15,sep="")
  jday <- yday(date)/365
  sizes <- c(vbgrHerring(jday-herring_hatch),vbgrHerring(jday-herring_hatch+1),
             vbgrHerring(jday-herring_hatch+2),vbgrHerring(jday-herring_hatch+3),
             vbgrHerring(jday-herring_hatch+4),vbgrHerring(jday-herring_hatch+5),
             vbgrHerring(jday-herring_hatch+6),vbgrHerring(jday-herring_hatch+7),
             vbgrHerring(jday-herring_hatch+8))*10
  ages <- 0:8
  herring.i <- N_ageHerring(Herring,ages,dada$Year[1])[jday*365,1:length(ages)]
  
  sd<- c(vbgr.sdHerring(jday-herring_hatch,100000),vbgr.sdHerring(jday-herring_hatch+1,100000),
         vbgr.sdHerring(jday-herring_hatch+2,100000),vbgr.sdHerring(jday-herring_hatch+3,100000),
         vbgr.sdHerring(jday-herring_hatch+4,100000),vbgr.sdHerring(jday-herring_hatch+5,100000),
         vbgr.sdHerring(jday-herring_hatch+6,100000),vbgr.sdHerring(jday-herring_hatch+7,100000),
         vbgr.sdHerring(jday-herring_hatch+7,100000))*10
  if (min(sd)<0){
    sd[sd] <-rep(0,length(sd[sd<0]))
  }
  available <- c(floor(rnorm(herring.i[1],sizes[1],sd[1])/width)*width+width/2,floor(rnorm(herring.i[2],sizes[2],sd[2])/width)*width+width/2,
                 floor(rnorm(herring.i[3],sizes[3],sd[3])/width)*width+width/2,floor(rnorm(herring.i[4],sizes[4],sd[4])/width)*width+width/2,
                 floor(rnorm(herring.i[5],sizes[5],sd[5])/width)*width+width/2,floor(rnorm(herring.i[6],sizes[6],sd[6])/width)*width+width/2,
                 floor(rnorm(herring.i[7],sizes[7],sd[7])/width)*width+width/2,floor(rnorm(herring.i[8],sizes[8],sd[8])/width)*width+width/2,
                 floor(rnorm(herring.i[9],sizes[9],sd[9])/width)*width+width/2)
  available[available<0] <- rep(width/2,length(available[available<0]))
  
  
  avail <- as.data.frame(table(available))
  colnames(avail) <- c("l.class", "available")
  eat <- as.data.frame(table(dada$size))
  colnames(eat) <- c("l.class", "eaten")
  
  ava.idx <- unique(df.seal$l.class) %in% avail$l.class
  eat.idx <- unique(df.seal$l.class) %in% eat$l.class
  i.idx <- (((i-1)*n.lengths+1):(i*n.lengths))
  df.seal$n[i.idx[ava.idx]] <- avail$available
  df.seal$n_eaten[i.idx[eat.idx]] <- eat$eaten
  df.seal$n_scat[i.idx] <- rep(length(unique(dada$Scat)),n.lengths)
  #  df$location <- rep(dada$Site[1],n.lengths)
  df.seal$year[i.idx] <- rep(dada$Year[1],n.lengths)
  
 
}
df.seal$odds <- df.seal$n_eaten/df.seal$n
#df.seal$n[df.seal$n==0] <- 1
df.seal <- df.seal %>% filter(n>0)
plot(df.seal$l.class,df.seal$odds)


# 2015
df.15 <- df.seal %>% filter(year==2015) 
tot.scat.15 <- sum(aggregate(n_scat~sampling,data=df.15,FUN=mean)$n_scat)
mean(df.15$odds)
df.15$weighted.odds <- (df.15$n_eaten/df.15$n)*df.15$n_scat
df.fin.15 <- aggregate(weighted.odds~l.class,data=df.15,FUN=sum)
df.fin.15$weighted.odds <- df.fin.15$weighted.odds/tot.scat.15
mean(df.fin.15$weighted.odds)
# 2016
df.16 <- df.seal %>% filter(year==2016) 
tot.scat.16 <- sum(aggregate(n_scat~sampling,data=df.16,FUN=mean)$n_scat)
mean(df.16$odds)
df.16$weighted.odds <- (df.16$n_eaten/df.16$n)*df.16$n_scat
df.fin.16 <- aggregate(weighted.odds~l.class,data=df.16,FUN=sum)
df.fin.16$weighted.odds <- df.fin.16$weighted.odds/tot.scat.16
mean(df.fin.16$weighted.odds)
# 2017
df.17 <- df.seal %>% filter(year==2017) 
tot.scat.17 <- sum(aggregate(n_scat~sampling,data=df.17,FUN=mean)$n_scat)
mean(df.17$odds)
df.17$weighted.odds <- (df.17$n_eaten/df.17$n)*df.17$n_scat
df.fin.17 <- aggregate(weighted.odds~l.class,data=df.17,FUN=sum)
df.fin.17$weighted.odds <- df.fin.17$weighted.odds/tot.scat.17
mean(df.fin.17$weighted.odds)
# All together - df.fin
tot.scat <- sum(aggregate(n_scat~sampling,data=df.seal,FUN=mean)$n_scat)
mean(df.seal$odds)
df.seal$weighted.odds <- (df.seal$n_eaten/df.seal$n)*df.seal$n_scat
df.fin.seal <- aggregate(weighted.odds~l.class,data=df.seal,FUN=mean)
df.fin.seal$weighted.odds <- df.fin.seal$weighted.odds*table(df.seal$l.class)[1]/tot.scat
mean(df.fin.seal$weighted.odds)

plot(df.fin.15$l.class,df.fin.15$weighted.odds,xlim=c(0,600),ylim=c(0,0.0002),pch=19,
     xlab = 'cm',ylab = 'preference',col='blue')
points(df.fin.16$l.class,df.fin.16$weighted.odds,pch=19,col='red')
points(df.fin.17$l.class,df.fin.17$weighted.odds,pch=19,col='green')
points(df.fin.seal$l.class,df.fin.seal$weighted.odds,pch=19)
legend(450,0.015,legend = c('2015','2016','2017','All'),pch=19,
       col=c('blue','red','green','black'))

scats <- data.frame(l.class = sort(unique(gseal.herring$size)),
                    n = as.numeric(table(gseal.herring$size)))

df.fin.seal <- left_join(df.fin.seal,scats,by="l.class")
df.fin.seal$n[is.na(df.fin.seal$n)] <- 0
df.fin.seal$weight <- df.fin.seal$n+1





df.fin.seal <- df.fin.seal %>% filter(l.class<350 & weighted.odds>0)
#df.fin.seal$weighted.odds[df.fin.seal$l.class>355] <- max(df.fin.seal$weighted.odds)
#df.fin.seal <- df.fin.seal %>% filter(l.class<360 & l.class>150)
sigmoid <- function(x, maximum,slope, midpoint) {
  #maximum <- max(df.fin.seal$weighted.odds) 
  maximum / (1 + exp(-slope * (x - midpoint))) # set maximum to the highest remaining observation
}


obj_fn <- function(par) {
  s <- exp(par[["logS"]])
  mid <- exp(par[["logMid"]])
  maximum <- exp(par[["logMax"]])
  pred <- sigmoid(df.fin.seal$l.class,maximum,s,mid)
  sum((df.fin.seal$weighted.odds - pred)^2)
}
starts <- expand.grid(logMax =log(max(df.fin.seal$weighted.odds)),logS = seq(log(0.02), log(0.1), 0.01), logMid = seq(log(250), log(350), 0.005))
results <- lapply(1:nrow(starts), function(i) {
  par <- as.list(starts[i,])
  opt <- try(nlminb(par, objective = obj_fn), silent = TRUE)
  if (!inherits(opt, "try-error")) {
    return(list(par = exp(opt$par), obj = opt$objective))
  } else {
    return(NULL)
  }
})
best <- results[[which.min(sapply(results, function(x) if (!is.null(x)) x$obj else Inf))]]
best


sigmoid(1:400,max(df.fin.seal$weighted.odds),0.03,300)
seal_param <- best$par

seal_fit <- sigmoid((0:5000)/10,seal_param[1],seal_param[2],seal_param[3])
plot(df.fin.seal$l.class,df.fin.seal$weighted.odds,ylab = "preference",
     xlab="length [mm]",main="Grey Seal",)
lines((0:5000)/10,seal_fit,
      lwd=2,col="red")
#lines((0:5000)/10,sigmoid((0:5000)/10,0.04,300))


make_my_function <- function(seal_param) {
  maximum <- seal_param[1]   
  slope <- seal_param[2]    
  midPoint <- seal_param[3]      

  function(x) { # x is length in mm
    pref <-   maximum / (1 + exp(-slope * (x - midPoint))) # set maximum to the highest remaining observation
    pref <- pref/max(pref)
    return(pref)
  }
}

# Create the function
seal_herring.pref <- make_my_function(seal_param)


#####
# cormorant preference for herring
#####
herring_corm <- cormorant %>% filter(ART=="herring")
herring_corm$dmy <- paste(herring_corm$year,'-',herring_corm$month,'-',herring_corm$day,sep="")
width <- 20 # mm. width of the size classes
#gseal$size <- floor(gseal$FL.with.SCF/50)*50
herring_corm$size <- floor(herring_corm$FISKLGD/width)*width+width/2
samplings <- unique(herring_corm$dmy)
month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
           "Sep","Oct","Nov","Dec")
vec.lengths <- seq(width/2,500+width/2,by=width)
n.lengths <- length(vec.lengths)

df.corm <- data.frame(l.class = rep(vec.lengths,length(samplings)),
                      n_eaten = rep(0,n.lengths*length(samplings)),
                      n =  rep(0,n.lengths*length(samplings)),
                      sampling = rep(samplings,each=n.lengths),
                      n_scat =  rep(0,n.lengths*length(samplings)),
                      year =  rep(0,n.lengths*length(samplings)))

herring_hatch <- mean(c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365)
start_time <- Sys.time() #NB takes 5 min
for (i in 1:length(samplings)){
  dada <- herring_corm %>% filter(dmy==samplings[i])
  n_sca <- length(unique(dada$GYLP))
  date <- paste(dada$year[1],'-',which(dada$month[1]==month)[[1]],'-',dada$day[1],sep="")
  jday <- yday(date)/365
  sizes <- c(vbgrHerring(jday-herring_hatch),vbgrHerring(jday-herring_hatch+1),
             vbgrHerring(jday-herring_hatch+2),vbgrHerring(jday-herring_hatch+3),
             vbgrHerring(jday-herring_hatch+4),vbgrHerring(jday-herring_hatch+5),
             vbgrHerring(jday-herring_hatch+6),vbgrHerring(jday-herring_hatch+7),
             vbgrHerring(jday-herring_hatch+8))*10
  ages <- 0:8
  herring.i <- N_ageHerring(Herring,ages,dada$year[1])[jday*365,1:length(ages)]
  sd<- c(vbgr.sdHerring(jday-herring_hatch,100000),vbgr.sdHerring(jday-herring_hatch+1,100000),
         vbgr.sdHerring(jday-herring_hatch+2,100000),vbgr.sdHerring(jday-herring_hatch+3,100000),
         vbgr.sdHerring(jday-herring_hatch+4,100000),vbgr.sdHerring(jday-herring_hatch+5,100000),
         vbgr.sdHerring(jday-herring_hatch+6,100000),vbgr.sdHerring(jday-herring_hatch+7,100000),
         vbgr.sdHerring(jday-herring_hatch+8,100000))*10
  
  if (min(sd)<0){
    sd[sd<0] <-rep(0,length(sd[sd<0]))
  }
  available <- c(floor(rnorm(herring.i[1],sizes[1],sd[1])/width)*width+width/2,floor(rnorm(herring.i[2],sizes[2],sd[2])/width)*width+width/2,
                 floor(rnorm(herring.i[3],sizes[3],sd[3])/width)*width+width/2,floor(rnorm(herring.i[4],sizes[4],sd[4])/width)*width+width/2,
                 floor(rnorm(herring.i[5],sizes[5],sd[5])/width)*width+width/2,floor(rnorm(herring.i[6],sizes[6],sd[6])/width)*width+width/2,
                 floor(rnorm(herring.i[7],sizes[7],sd[7])/width)*width+width/2,floor(rnorm(herring.i[8],sizes[8],sd[8])/width)*width+width/2,
                 floor(rnorm(herring.i[9],sizes[9],sd[9])/width)*width+width/2)
  available[available<0] <- rep(width/2,length(available[available<0]))
  
  avail <- as.data.frame(table(available))
  colnames(avail) <- c("l.class", "available")
  eat <- as.data.frame(table(dada$size))
  colnames(eat) <- c("l.class", "eaten")
  
  ava.idx <- unique(df.corm$l.class) %in% avail$l.class
  eat.idx <- unique(df.corm$l.class) %in% eat$l.class
  i.idx <- (((i-1)*n.lengths+1):(i*n.lengths))
  df.corm$n[i.idx[ava.idx]] <- avail$available
  df.corm$n_eaten[i.idx[eat.idx]] <- eat$eaten
  df.corm$n_scat[i.idx] <- rep(n_sca,n.lengths)
  #  df$location <- rep(dada$Site[1],n.lengths)
  df.corm$year[i.idx] <- rep(dada$year[1],n.lengths)
}
end_time <- Sys.time()
elapsed <- end_time - start_time
print(elapsed)
df.corm$odds <- df.corm$n_eaten/df.corm$n
df.corm <- df.corm %>% filter(n>0)
plot(df.corm$l.class,df.corm$odds)
# Remove outliers - NB put in again when more data
#df.corm <- df.corm %>% filter(odds<0.0065)
# 2015
df.15 <- df.corm %>% filter(year==1992) 
tot.scat.15 <- sum(aggregate(n_scat~sampling,data=df.15,FUN=mean)$n_scat)
mean(df.15$odds)
df.15$weighted.odds <- (df.15$n_eaten/df.15$n)*df.15$n_scat
df.fin.15 <- aggregate(weighted.odds~l.class,data=df.15,FUN=sum)
df.fin.15$weighted.odds <- df.fin.15$weighted.odds/tot.scat.15
mean(df.fin.15$weighted.odds)
# 2016
df.16 <- df.corm %>% filter(year==1993) 
tot.scat.16 <- sum(aggregate(n_scat~sampling,data=df.16,FUN=mean)$n_scat)
mean(df.16$odds)
df.16$weighted.odds <- (df.16$n_eaten/df.16$n)*df.16$n_scat
df.fin.16 <- aggregate(weighted.odds~l.class,data=df.16,FUN=sum)
df.fin.16$weighted.odds <- df.fin.16$weighted.odds/tot.scat.16
mean(df.fin.16$weighted.odds)
# 2017
df.17 <- df.corm %>% filter(year==1994) 
tot.scat.17 <- sum(aggregate(n_scat~sampling,data=df.17,FUN=mean)$n_scat)
mean(df.17$odds)
df.17$weighted.odds <- (df.17$n_eaten/df.17$n)*df.17$n_scat
df.fin.17 <- aggregate(weighted.odds~l.class,data=df.17,FUN=sum)
df.fin.17$weighted.odds <- df.fin.17$weighted.odds/tot.scat.17
mean(df.fin.17$weighted.odds)
# All together - df.fin
tot.scat <- sum(aggregate(n_scat~sampling,data=df.corm,FUN=mean)$n_scat)
mean(df.corm$odds)
df.corm$weighted.odds <- (df.corm$n_eaten/df.corm$n)*df.corm$n_scat
df.fin.corm <- aggregate(weighted.odds~l.class,data=df.corm,FUN=mean)
df.fin.corm$weighted.odds <- df.fin.corm$weighted.odds*table(df.corm$l.class)[1]/tot.scat
mean(df.fin.corm$weighted.odds)

plot(df.fin.15$l.class,df.fin.15$weighted.odds,xlim=c(0,500),
     ylim=c(0,max(c(df.fin.15$weighted.odds,df.fin.16$weighted.odds,df.fin.17$weighted.odds))),pch=19,
     xlab = 'cm',ylab = 'preference',col='blue')
points(df.fin.16$l.class,df.fin.16$weighted.odds,pch=19,col='red')
points(df.fin.17$l.class,df.fin.17$weighted.odds,pch=19,col='green')
points(df.fin.corm$l.class,df.fin.corm$weighted.odds,pch=19)
legend(350,max(c(df.fin.15$weighted.odds,df.fin.16$weighted.odds,df.fin.17$weighted.odds))
       ,legend = c('2015','2016','2017','All'),pch=19,
       col=c('blue','red','green','black'))

scats <- data.frame(l.class = sort(unique(cormorant$size)),
                    n = as.numeric(table(cormorant$size)))

df.fin.corm <- left_join(df.fin.corm,scats,by="l.class")
df.fin.corm$n[is.na(df.fin.corm$n)] <- 0
df.fin.corm$weight <- df.fin.corm$n+1

bell_curve <- function(x, A, mu, sigma) {
  A * exp(- (x - mu)^2 / (2 * sigma^2))
}

df.fin.corm <- df.fin.corm %>% filter(weighted.odds<(2*10^-5))
fit <- nls(weighted.odds ~ bell_curve(l.class,A,mu,sigma),
           data = df.fin.corm,
           start = list(A = max(df.fin.corm$weighted.odds), mu = 200, sigma = 50))

corm_param <- coef(fit)

corm_fit <- bell_curve((0:5000)/10,corm_param[1],corm_param[2],corm_param[3])
plot(df.fin.corm$l.class,df.fin.corm$weighted.odds,ylab = "preference",
     xlab="length [mm]",main="Cormorant",xlim=c(0,500))
lines((0:5000)/10,corm_fit,
      lwd=2,col="orange")
#lines((0:7000)/10,bell_curve((0:7000)/10,max(df.fin.corm$weighted.odds),180,65),
#      lwd=2,col="darkblue")

sum(abs(df.fin.corm$weighted.odds-bell_curve(df.fin.corm$l.class,corm_param[1],corm_param[2],corm_param[3])))
sum(abs(df.fin.corm$weighted.odds-bell_curve(df.fin.corm$l.class,max(df.fin.corm$weighted.odds),180,65)))


plot((0:5000)/10,corm_fit/max(corm_fit),
     lwd=3,col="darkblue",ylab = "preference",
     xlab="length [mm]",main="herring",type='l',ylim=c(0,1.5))
lines((0:5000)/10,seal_fit/max(seal_fit),
      lwd=3,col='red')
points(df.fin.corm$l.class,df.fin.corm$weighted.odds/max(corm_fit),pc=19,cex=1,col="darkblue")
points(df.fin.seal$l.class,df.fin.seal$weighted.odds/max(seal_fit),pc=19,cex=1,col="red")
legend(355,0.5,legend=c('Seal','Cormorant'),
       pch=19,col=c('red','darkblue'))

make_my_function <- function(corm_param) {
  A <- corm_param[1]   
  mu <- corm_param[2]    
  sigma <- corm_param[3]      

  function(x) { # x is length in mm
    pref <- A * exp(- (x - mu)^2 / (2 * sigma^2))
    pref <- pref/max(pref)
    return(pref)
  }
}

# Create the function
corm_herring.pref <- make_my_function(corm_param)

#####

# TIP: Run Flatfish together
# Seal preference for flatfish (flounder, plaice, and dab)
#####
gseal.flatfish <- gseal %>% filter(Species %in% c("flatfish"))
gseal.flatfish$my <- paste(gseal.flatfish$Month,'-',gseal.flatfish$Year)
width <- 30 # mm. width of the size classes
gseal.flatfish$size <- floor(gseal.flatfish$FL.with.SCF/width)*width+width/2
samplings <- unique(gseal.flatfish$my)
month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
           "Sep","Oct","Nov","Dec")
vec.lengths <- seq(width/2,700+width/2,by=width)
n.lengths <- length(vec.lengths)
df.seal <- data.frame(l.class = rep(vec.lengths,length(samplings)),
                      n_eaten = rep(0,n.lengths*length(samplings)),
                      n =  rep(0,n.lengths*length(samplings)),
                      sampling = rep(samplings,each=n.lengths),
                      n_scat =  rep(0,n.lengths*length(samplings)),
                      year =  rep(0,n.lengths*length(samplings)))

Flounder_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Plaice_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Dab_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

for (i in 1:length(samplings)){
  dada <- gseal.flatfish %>% filter(my==samplings[i])
  n_sca <- length(unique(dada$Scat))
  date <- paste(dada$Year[1],'-',which(dada$Month[1]==month)[[1]],'-',15,sep="")
  jday <- yday(date)/365
  sizesFlounder <- c(vbgrFlounder(jday-Flounder_hatch),vbgrFlounder(jday-Flounder_hatch+1),
                     vbgrFlounder(jday-Flounder_hatch+2),vbgrFlounder(jday-Flounder_hatch+3),
                     vbgrFlounder(jday-Flounder_hatch+4),vbgrFlounder(jday-Flounder_hatch+5),
                     vbgrFlounder(jday-Flounder_hatch+6),vbgrFlounder(jday-Flounder_hatch+7),
                     vbgrFlounder(jday-Flounder_hatch+8),vbgrFlounder(jday-Flounder_hatch+9),
                     vbgrFlounder(jday-Flounder_hatch+10))*10
  sizesPlaice <- c(vbgrPlaice(jday-Plaice_hatch),vbgrPlaice(jday-Plaice_hatch+1),
                     vbgrPlaice(jday-Plaice_hatch+2),vbgrPlaice(jday-Plaice_hatch+3),
                     vbgrPlaice(jday-Plaice_hatch+4),vbgrPlaice(jday-Plaice_hatch+5),
                     vbgrPlaice(jday-Plaice_hatch+6),vbgrPlaice(jday-Plaice_hatch+7),
                     vbgrPlaice(jday-Plaice_hatch+8),vbgrPlaice(jday-Plaice_hatch+9),
                     vbgrPlaice(jday-Plaice_hatch+10))*10
  sizesDab <- c(vbgrDab(jday-Dab_hatch),vbgrDab(jday-Dab_hatch+1),
                     vbgrDab(jday-Dab_hatch+2),vbgrDab(jday-Dab_hatch+3),
                     vbgrDab(jday-Dab_hatch+4),vbgrDab(jday-Dab_hatch+5),
                     vbgrDab(jday-Dab_hatch+6),vbgrDab(jday-Dab_hatch+7),
                     vbgrDab(jday-Dab_hatch+8),vbgrDab(jday-Dab_hatch+9),
                     vbgrDab(jday-Dab_hatch+10))*10
  fl.i <- Flatfish %>% filter(yday==round(jday*365))
  fl.i$Age[fl.i$Age>=11] <- 10+jday
  flatfish.i <- aggregate(N~Age+species,data=fl.i,FUN=sum)
  
  sdFlounder<- c(vbgr.sdFlounder(jday-Flounder_hatch,100000),vbgr.sdFlounder(jday-Flounder_hatch+1,100000),
         vbgr.sdFlounder(jday-Flounder_hatch+2,100000),vbgr.sdFlounder(jday-Flounder_hatch+3,100000),
         vbgr.sdFlounder(jday-Flounder_hatch+4,100000),vbgr.sdFlounder(jday-Flounder_hatch+5,100000),
         vbgr.sdFlounder(jday-Flounder_hatch+6,100000),vbgr.sdFlounder(jday-Flounder_hatch+7,100000),
         vbgr.sdFlounder(jday-Flounder_hatch+8,100000),vbgr.sdFlounder(jday-Flounder_hatch+9,100000),
         vbgr.sdFlounder(jday-Flounder_hatch+10,100000))*10
  sdPlaice<- c(vbgr.sdPlaice(jday-Plaice_hatch,100000),vbgr.sdPlaice(jday-Plaice_hatch+1,100000),
                 vbgr.sdPlaice(jday-Plaice_hatch+2,100000),vbgr.sdPlaice(jday-Plaice_hatch+3,100000),
                 vbgr.sdPlaice(jday-Plaice_hatch+4,100000),vbgr.sdPlaice(jday-Plaice_hatch+5,100000),
                 vbgr.sdPlaice(jday-Plaice_hatch+6,100000),vbgr.sdPlaice(jday-Plaice_hatch+7,100000),
                 vbgr.sdPlaice(jday-Plaice_hatch+8,100000),vbgr.sdPlaice(jday-Plaice_hatch+9,100000),
                 vbgr.sdPlaice(jday-Plaice_hatch+10,100000))*10
  sdDab<- c(vbgr.sdDab(jday-Dab_hatch,100000),vbgr.sdDab(jday-Dab_hatch+1,100000),
                 vbgr.sdDab(jday-Dab_hatch+2,100000),vbgr.sdDab(jday-Dab_hatch+3,100000),
                 vbgr.sdDab(jday-Dab_hatch+4,100000),vbgr.sdDab(jday-Dab_hatch+5,100000),
                 vbgr.sdDab(jday-Dab_hatch+6,100000),vbgr.sdDab(jday-Dab_hatch+7,100000),
                 vbgr.sdDab(jday-Dab_hatch+8,100000),vbgr.sdDab(jday-Dab_hatch+9,100000),
                 vbgr.sdDab(jday-Dab_hatch+10,100000))*10
  
  
  if (min(sdFlounder)<0){
    sdFlounder[sdFlounder] <-rep(0,length(sdFlounder[sdFlounder<0]))
  }
  if (min(sdPlaice)<0){
    sdPlaice[sdPlaice] <-rep(0,length(sdPlaice[sdPlaice<0]))
  }
  if (min(sdDab)<0){
    sdDab[sdDab] <-rep(0,length(sdDab[sdDab<0]))
  }
  
  
  
  
  availableFlounder <- c(floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][1],sizesFlounder[1],sdFlounder[1])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][2],sizesFlounder[2],sdFlounder[2])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][3],sizesFlounder[3],sdFlounder[3])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][4],sizesFlounder[4],sdFlounder[4])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][5],sizesFlounder[5],sdFlounder[5])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][6],sizesFlounder[6],sdFlounder[6])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][7],sizesFlounder[7],sdFlounder[7])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][8],sizesFlounder[8],sdFlounder[8])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][9],sizesFlounder[9],sdFlounder[9])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][10],sizesFlounder[10],sdFlounder[10])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][11],sizesFlounder[11],sdFlounder[11])/width)*width+width/2)
  availableFlounder[availableFlounder<0] <- rep(width/2,length(availableFlounder[availableFlounder<0]))
  availablePlaice <- c(floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][1],sizesPlaice[1],sdPlaice[1])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][2],sizesPlaice[2],sdPlaice[2])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][3],sizesPlaice[3],sdPlaice[3])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][4],sizesPlaice[4],sdPlaice[4])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][5],sizesPlaice[5],sdPlaice[5])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][6],sizesPlaice[6],sdPlaice[6])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][7],sizesPlaice[7],sdPlaice[7])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][8],sizesPlaice[8],sdPlaice[8])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][9],sizesPlaice[9],sdPlaice[9])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][10],sizesPlaice[10],sdPlaice[10])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][11],sizesPlaice[11],sdPlaice[11])/width)*width+width/2)
  availablePlaice[availablePlaice<0] <- rep(width/2,length(availablePlaice[availablePlaice<0]))
  availableDab <- c(floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][1],sizesDab[1],sdDab[1])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][2],sizesDab[2],sdDab[2])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][3],sizesDab[3],sdDab[3])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][4],sizesDab[4],sdDab[4])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][5],sizesDab[5],sdDab[5])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][6],sizesDab[6],sdDab[6])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][7],sizesDab[7],sdDab[7])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][8],sizesDab[8],sdDab[8])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][9],sizesDab[9],sdDab[9])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][10],sizesDab[10],sdDab[10])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][11],sizesDab[11],sdDab[11])/width)*width+width/2)
  availableDab[availableDab<0] <- rep(width/2,length(availableDab[availableDab<0]))
  available <- c(availableFlounder,availablePlaice,availableDab)
  
  avail <- as.data.frame(table(available))
  colnames(avail) <- c("l.class", "available")
  eat <- as.data.frame(table(dada$size))
  colnames(eat) <- c("l.class", "eaten")
  
  ava.idx <- unique(df.seal$l.class) %in% avail$l.class
  eat.idx <- unique(df.seal$l.class) %in% eat$l.class
  i.idx <- (((i-1)*n.lengths+1):(i*n.lengths))
  df.seal$n[i.idx[ava.idx]] <- avail$available
  df.seal$n_eaten[i.idx[eat.idx]] <- eat$eaten
  df.seal$n_scat[i.idx] <- rep(length(unique(dada$Scat)),n.lengths)
  #  df$location <- rep(dada$Site[1],n.lengths)
  df.seal$year[i.idx] <- rep(dada$Year[1],n.lengths)
  
  print(i)
}
df.seal$odds <- df.seal$n_eaten/df.seal$n
#df.seal$n[df.seal$n==0] <- 1
df.seal <- df.seal %>% filter(n>0)
plot(df.seal$l.class,df.seal$odds)


# 2015
df.15 <- df.seal %>% filter(year==2015) 
tot.scat.15 <- sum(aggregate(n_scat~sampling,data=df.15,FUN=mean)$n_scat)
mean(df.15$odds)
df.15$weighted.odds <- (df.15$n_eaten/df.15$n)*df.15$n_scat
df.fin.15 <- aggregate(weighted.odds~l.class,data=df.15,FUN=sum)
df.fin.15$weighted.odds <- df.fin.15$weighted.odds/tot.scat.15
mean(df.fin.15$weighted.odds)
# 2016
df.16 <- df.seal %>% filter(year==2016) 
tot.scat.16 <- sum(aggregate(n_scat~sampling,data=df.16,FUN=mean)$n_scat)
mean(df.16$odds)
df.16$weighted.odds <- (df.16$n_eaten/df.16$n)*df.16$n_scat
df.fin.16 <- aggregate(weighted.odds~l.class,data=df.16,FUN=sum)
df.fin.16$weighted.odds <- df.fin.16$weighted.odds/tot.scat.16
mean(df.fin.16$weighted.odds)
# 2017
df.17 <- df.seal %>% filter(year==2017) 
tot.scat.17 <- sum(aggregate(n_scat~sampling,data=df.17,FUN=mean)$n_scat)
mean(df.17$odds)
df.17$weighted.odds <- (df.17$n_eaten/df.17$n)*df.17$n_scat
df.fin.17 <- aggregate(weighted.odds~l.class,data=df.17,FUN=sum)
df.fin.17$weighted.odds <- df.fin.17$weighted.odds/tot.scat.17
mean(df.fin.17$weighted.odds)
# All together - df.fin
tot.scat <- sum(aggregate(n_scat~sampling,data=df.seal,FUN=mean)$n_scat)
mean(df.seal$odds)
df.seal$weighted.odds <- (df.seal$n_eaten/df.seal$n)*df.seal$n_scat
df.fin.seal <- aggregate(weighted.odds~l.class,data=df.seal,FUN=mean)
df.fin.seal$weighted.odds <- df.fin.seal$weighted.odds*table(df.seal$l.class)[1]/tot.scat
mean(df.fin.seal$weighted.odds)

plot(df.fin.15$l.class,df.fin.15$weighted.odds,xlim=c(0,600),ylim=c(0,0.4),pch=19,
     xlab = 'cm',ylab = 'preference',col='blue')
points(df.fin.16$l.class,df.fin.16$weighted.odds,pch=19,col='red')
points(df.fin.17$l.class,df.fin.17$weighted.odds,pch=19,col='green')
points(df.fin.seal$l.class,df.fin.seal$weighted.odds,pch=19)
legend(0,0.4,legend = c('2015','2016','2017','All'),pch=19,
       col=c('blue','red','green','black'))







df.fin.seal <- df.fin.seal %>% filter(l.class<400)
#df.fin.seal$weighted.odds[df.fin.seal$l.class>355] <- max(df.fin.seal$weighted.odds)
#df.fin.seal <- df.fin.seal %>% filter(l.class<360 & l.class>150)
sigmoid <- function(x, maximum,slope, midpoint) {
  #maximum <- max(df.fin.seal$weighted.odds) 
  maximum / (1 + exp(-slope * (x - midpoint))) # set maximum to the highest remaining observation
}


obj_fn <- function(par) {
  s <- exp(par[["logS"]])
  mid <- exp(par[["logMid"]])
  maximum <- exp(par[["logMax"]])
  pred <- sigmoid(df.fin.seal$l.class,maximum,s,mid)
  sum((df.fin.seal$weighted.odds - pred)^2)
}
starts <- expand.grid(logMax =log(max(df.fin.seal$weighted.odds)),logS = seq(log(0.02), log(0.1), 0.01), logMid = seq(log(250), log(350), 0.005))
results <- lapply(1:nrow(starts), function(i) {
  par <- as.list(starts[i,])
  opt <- try(nlminb(par, objective = obj_fn), silent = TRUE)
  if (!inherits(opt, "try-error")) {
    return(list(par = exp(opt$par), obj = opt$objective))
  } else {
    return(NULL)
  }
})
best <- results[[which.min(sapply(results, function(x) if (!is.null(x)) x$obj else Inf))]]
best


seal_param <- best$par

seal_fit <- sigmoid((0:5000)/10,seal_param[1],seal_param[2],seal_param[3])
plot(df.fin.seal$l.class,df.fin.seal$weighted.odds,ylab = "preference",
     xlab="length [mm]",main="Grey Seal",xlim=c(0,700),ylim=c(0,0.04))
lines((0:5000)/10,seal_fit,
      lwd=2,col="red")
lines((0:5000)/10,sigmoid((0:5000)/10,0.03745708,0.04,323))
seal_fit <- sigmoid((0:5000)/10,0.03745708,0.04,323)
seal_param <- c(0.001,0.04,310)

make_my_function <- function(seal_param) {
  maximum <- seal_param[1]   
  slope <- seal_param[2]    
  midPoint <- seal_param[3]      
  
  function(x) { # x is length in mm
    pref <-   maximum / (1 + exp(-slope * (x - midPoint))) # set maximum to the highest remaining observation
    pref <- pref/max(pref)
    return(pref)
  }
}

# Create the function
seal_flatfish.pref <- make_my_function(seal_param)
#####
# cormorant preference for flatfish
#####
flatfish_corm <- cormorant %>% filter(ART=="flatfish")
flatfish_corm$dmy <- paste(flatfish_corm$year,'-',flatfish_corm$month,'-',flatfish_corm$day,sep="")
width <- 20 # mm. width of the size classes
#gseal$size <- floor(gseal$FL.with.SCF/50)*50
flatfish_corm$size <- floor(flatfish_corm$FISKLGD/width)*width+width/2
samplings <- unique(flatfish_corm$dmy)
month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
           "Sep","Oct","Nov","Dec")
vec.lengths <- seq(width/2,700+width/2,by=width)
n.lengths <- length(vec.lengths)

df.corm <- data.frame(l.class = rep(vec.lengths,length(samplings)),
                      n_eaten = rep(0,n.lengths*length(samplings)),
                      n =  rep(0,n.lengths*length(samplings)),
                      sampling = rep(samplings,each=n.lengths),
                      n_scat =  rep(0,n.lengths*length(samplings)),
                      year =  rep(0,n.lengths*length(samplings)))

Flounder_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Plaice_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Dab_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

start_time <- Sys.time() #NB takes 5 min
for (i in 1:length(samplings)){
  dada <- flatfish_corm %>% filter(dmy==samplings[i])
  n_sca <- length(unique(dada$GYLP))
  date <- paste(dada$year[1],'-',which(dada$month[1]==month)[[1]],'-',dada$day[1],sep="")
  jday <- yday(date)/365
  sizesFlounder <- c(vbgrFlounder(jday-Flounder_hatch),vbgrFlounder(jday-Flounder_hatch+1),
                     vbgrFlounder(jday-Flounder_hatch+2),vbgrFlounder(jday-Flounder_hatch+3),
                     vbgrFlounder(jday-Flounder_hatch+4),vbgrFlounder(jday-Flounder_hatch+5),
                     vbgrFlounder(jday-Flounder_hatch+6),vbgrFlounder(jday-Flounder_hatch+7),
                     vbgrFlounder(jday-Flounder_hatch+8),vbgrFlounder(jday-Flounder_hatch+9),
                     vbgrFlounder(jday-Flounder_hatch+10))*10
  sizesPlaice <- c(vbgrPlaice(jday-Plaice_hatch),vbgrPlaice(jday-Plaice_hatch+1),
                   vbgrPlaice(jday-Plaice_hatch+2),vbgrPlaice(jday-Plaice_hatch+3),
                   vbgrPlaice(jday-Plaice_hatch+4),vbgrPlaice(jday-Plaice_hatch+5),
                   vbgrPlaice(jday-Plaice_hatch+6),vbgrPlaice(jday-Plaice_hatch+7),
                   vbgrPlaice(jday-Plaice_hatch+8),vbgrPlaice(jday-Plaice_hatch+9),
                   vbgrPlaice(jday-Plaice_hatch+10))*10
  sizesDab <- c(vbgrDab(jday-Dab_hatch),vbgrDab(jday-Dab_hatch+1),
                vbgrDab(jday-Dab_hatch+2),vbgrDab(jday-Dab_hatch+3),
                vbgrDab(jday-Dab_hatch+4),vbgrDab(jday-Dab_hatch+5),
                vbgrDab(jday-Dab_hatch+6),vbgrDab(jday-Dab_hatch+7),
                vbgrDab(jday-Dab_hatch+8),vbgrDab(jday-Dab_hatch+9),
                vbgrDab(jday-Dab_hatch+10))*10
  fl.i <- Flatfish %>% filter(yday==round(jday*365))
  fl.i$Age[fl.i$Age>=11] <- 10+jday
  flatfish.i <- aggregate(N~Age+species,data=fl.i,FUN=sum)
  
  sdFlounder<- c(vbgr.sdFlounder(jday-Flounder_hatch,100000),vbgr.sdFlounder(jday-Flounder_hatch+1,100000),
                 vbgr.sdFlounder(jday-Flounder_hatch+2,100000),vbgr.sdFlounder(jday-Flounder_hatch+3,100000),
                 vbgr.sdFlounder(jday-Flounder_hatch+4,100000),vbgr.sdFlounder(jday-Flounder_hatch+5,100000),
                 vbgr.sdFlounder(jday-Flounder_hatch+6,100000),vbgr.sdFlounder(jday-Flounder_hatch+7,100000),
                 vbgr.sdFlounder(jday-Flounder_hatch+8,100000),vbgr.sdFlounder(jday-Flounder_hatch+9,100000),
                 vbgr.sdFlounder(jday-Flounder_hatch+10,100000))*10
  sdPlaice<- c(vbgr.sdPlaice(jday-Plaice_hatch,100000),vbgr.sdPlaice(jday-Plaice_hatch+1,100000),
               vbgr.sdPlaice(jday-Plaice_hatch+2,100000),vbgr.sdPlaice(jday-Plaice_hatch+3,100000),
               vbgr.sdPlaice(jday-Plaice_hatch+4,100000),vbgr.sdPlaice(jday-Plaice_hatch+5,100000),
               vbgr.sdPlaice(jday-Plaice_hatch+6,100000),vbgr.sdPlaice(jday-Plaice_hatch+7,100000),
               vbgr.sdPlaice(jday-Plaice_hatch+8,100000),vbgr.sdPlaice(jday-Plaice_hatch+9,100000),
               vbgr.sdPlaice(jday-Plaice_hatch+10,100000))*10
  sdDab<- c(vbgr.sdDab(jday-Dab_hatch,100000),vbgr.sdDab(jday-Dab_hatch+1,100000),
            vbgr.sdDab(jday-Dab_hatch+2,100000),vbgr.sdDab(jday-Dab_hatch+3,100000),
            vbgr.sdDab(jday-Dab_hatch+4,100000),vbgr.sdDab(jday-Dab_hatch+5,100000),
            vbgr.sdDab(jday-Dab_hatch+6,100000),vbgr.sdDab(jday-Dab_hatch+7,100000),
            vbgr.sdDab(jday-Dab_hatch+8,100000),vbgr.sdDab(jday-Dab_hatch+9,100000),
            vbgr.sdDab(jday-Dab_hatch+10,100000))*10
  
  
  if (min(sdFlounder)<0){
    sdFlounder[sdFlounder] <-rep(0,length(sdFlounder[sdFlounder<0]))
  }
  if (min(sdPlaice)<0){
    sdPlaice[sdPlaice] <-rep(0,length(sdPlaice[sdPlaice<0]))
  }
  if (min(sdDab)<0){
    sdDab[sdDab] <-rep(0,length(sdDab[sdDab<0]))
  }
  
  
  
  
  availableFlounder <- c(floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][1],sizesFlounder[1],sdFlounder[1])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][2],sizesFlounder[2],sdFlounder[2])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][3],sizesFlounder[3],sdFlounder[3])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][4],sizesFlounder[4],sdFlounder[4])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][5],sizesFlounder[5],sdFlounder[5])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][6],sizesFlounder[6],sdFlounder[6])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][7],sizesFlounder[7],sdFlounder[7])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][8],sizesFlounder[8],sdFlounder[8])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][9],sizesFlounder[9],sdFlounder[9])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][10],sizesFlounder[10],sdFlounder[10])/width)*width+width/2,
                         floor(rnorm(flatfish.i$N[flatfish.i$species=="flounder"][11],sizesFlounder[11],sdFlounder[11])/width)*width+width/2)
  availableFlounder[availableFlounder<0] <- rep(width/2,length(availableFlounder[availableFlounder<0]))
  availablePlaice <- c(floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][1],sizesPlaice[1],sdPlaice[1])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][2],sizesPlaice[2],sdPlaice[2])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][3],sizesPlaice[3],sdPlaice[3])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][4],sizesPlaice[4],sdPlaice[4])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][5],sizesPlaice[5],sdPlaice[5])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][6],sizesPlaice[6],sdPlaice[6])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][7],sizesPlaice[7],sdPlaice[7])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][8],sizesPlaice[8],sdPlaice[8])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][9],sizesPlaice[9],sdPlaice[9])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][10],sizesPlaice[10],sdPlaice[10])/width)*width+width/2,
                       floor(rnorm(flatfish.i$N[flatfish.i$species=="plaice"][11],sizesPlaice[11],sdPlaice[11])/width)*width+width/2)
  availablePlaice[availablePlaice<0] <- rep(width/2,length(availablePlaice[availablePlaice<0]))
  availableDab <- c(floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][1],sizesDab[1],sdDab[1])/width)*width+width/2,
                    floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][2],sizesDab[2],sdDab[2])/width)*width+width/2,
                    floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][3],sizesDab[3],sdDab[3])/width)*width+width/2,
                    floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][4],sizesDab[4],sdDab[4])/width)*width+width/2,
                    floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][5],sizesDab[5],sdDab[5])/width)*width+width/2,
                    floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][6],sizesDab[6],sdDab[6])/width)*width+width/2,
                    floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][7],sizesDab[7],sdDab[7])/width)*width+width/2,
                    floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][8],sizesDab[8],sdDab[8])/width)*width+width/2,
                    floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][9],sizesDab[9],sdDab[9])/width)*width+width/2,
                    floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][10],sizesDab[10],sdDab[10])/width)*width+width/2,
                    floor(rnorm(flatfish.i$N[flatfish.i$species=="dab"][11],sizesDab[11],sdDab[11])/width)*width+width/2)
  availableDab[availableDab<0] <- rep(width/2,length(availableDab[availableDab<0]))
  available <- c(availableFlounder,availablePlaice,availableDab)
  
  avail <- as.data.frame(table(available))
  colnames(avail) <- c("l.class", "available")
  eat <- as.data.frame(table(dada$size))
  colnames(eat) <- c("l.class", "eaten")
  
  ava.idx <- unique(df.corm$l.class) %in% avail$l.class
  eat.idx <- unique(df.corm$l.class) %in% eat$l.class
  i.idx <- (((i-1)*n.lengths+1):(i*n.lengths))
  df.corm$n[i.idx[ava.idx]] <- avail$available
  df.corm$n_eaten[i.idx[eat.idx]] <- eat$eaten
  df.corm$n_scat[i.idx] <- rep(n_sca,n.lengths)
  #  df$location <- rep(dada$Site[1],n.lengths)
  df.corm$year[i.idx] <- rep(dada$year[1],n.lengths)
  
  print(i)
}
end_time <- Sys.time()
elapsed <- end_time - start_time
print(elapsed)
df.corm$odds <- df.corm$n_eaten/df.corm$n
df.corm <- df.corm %>% filter(n>0)
plot(df.corm$l.class,df.corm$odds)
# Remove outliers - NB put in again when more data
#df.corm <- df.corm %>% filter(odds<0.0065)
# 2015
df.15 <- df.corm %>% filter(year==1992) 
tot.scat.15 <- sum(aggregate(n_scat~sampling,data=df.15,FUN=mean)$n_scat)
mean(df.15$odds)
df.15$weighted.odds <- (df.15$n_eaten/df.15$n)*df.15$n_scat
df.fin.15 <- aggregate(weighted.odds~l.class,data=df.15,FUN=sum)
df.fin.15$weighted.odds <- df.fin.15$weighted.odds/tot.scat.15
mean(df.fin.15$weighted.odds)
# 2016
df.16 <- df.corm %>% filter(year==1993) 
tot.scat.16 <- sum(aggregate(n_scat~sampling,data=df.16,FUN=mean)$n_scat)
mean(df.16$odds)
df.16$weighted.odds <- (df.16$n_eaten/df.16$n)*df.16$n_scat
df.fin.16 <- aggregate(weighted.odds~l.class,data=df.16,FUN=sum)
df.fin.16$weighted.odds <- df.fin.16$weighted.odds/tot.scat.16
mean(df.fin.16$weighted.odds)
# 2017
df.17 <- df.corm %>% filter(year==1994) 
tot.scat.17 <- sum(aggregate(n_scat~sampling,data=df.17,FUN=mean)$n_scat)
mean(df.17$odds)
df.17$weighted.odds <- (df.17$n_eaten/df.17$n)*df.17$n_scat
df.fin.17 <- aggregate(weighted.odds~l.class,data=df.17,FUN=sum)
df.fin.17$weighted.odds <- df.fin.17$weighted.odds/tot.scat.17
mean(df.fin.17$weighted.odds)
# All together - df.fin
tot.scat <- sum(aggregate(n_scat~sampling,data=df.corm,FUN=mean)$n_scat)
mean(df.corm$odds)
df.corm$weighted.odds <- (df.corm$n_eaten/df.corm$n)*df.corm$n_scat
df.fin.corm <- aggregate(weighted.odds~l.class,data=df.corm,FUN=mean)
df.fin.corm$weighted.odds <- df.fin.corm$weighted.odds*table(df.corm$l.class)[1]/tot.scat
mean(df.fin.corm$weighted.odds)

plot(df.fin.15$l.class,df.fin.15$weighted.odds,xlim=c(0,500),
     ylim=c(0,max(c(df.fin.15$weighted.odds,df.fin.16$weighted.odds,df.fin.17$weighted.odds))),pch=19,
     xlab = 'cm',ylab = 'preference',col='blue')
points(df.fin.16$l.class,df.fin.16$weighted.odds,pch=19,col='red')
points(df.fin.17$l.class,df.fin.17$weighted.odds,pch=19,col='green')
points(df.fin.corm$l.class,df.fin.corm$weighted.odds,pch=19)
legend(350,max(c(df.fin.15$weighted.odds,df.fin.16$weighted.odds,df.fin.17$weighted.odds))
       ,legend = c('2015','2016','2017','All'),pch=19,
       col=c('blue','red','green','black'))


df.fin.corm <- df.fin.corm %>% filter(weighted.odds<0.02)
double_sigmoid <- function(x, maximum, slope1, midPoint1, slope2, midPointDistance) {
  (maximum / (1 + exp(-slope1 * (x - midPoint1)))) - 
    (maximum / (1 + exp(-slope2 * (x - (midPoint1 + midPointDistance)))))
}

par <- list(logMax = log(0.013),logS1=log(0.12),logMid=log(70),logS2=log(0.09),logDis = log(150))


obj_fn <- function(par) {
  max <- exp(par[["logMax"]])
  s1 <- exp(par[["logS1"]])
  mid <- exp(par[["logMid"]])
  s2 <- exp(par[["logS2"]])
  dis <- exp(par[["logDis"]])
  pred <- double_sigmoid(df.fin.corm$l.class,max,s1,mid,s2,dis)
  sum((df.fin.corm$weighted.odds - pred)^2)
}
# df.fin.seal$weight*
opt <- nlminb(par, objective = obj_fn)
corm_param <- exp(opt$par)
corm_fit <- double_sigmoid((0:5000)/10,corm_param[1],corm_param[2],corm_param[3],corm_param[4],corm_param[5])

plot(df.fin.corm$l.class,df.fin.corm$weighted.odds,ylab = "preference",
     xlab="length [mm]",main="Cormorant",xlim=c(0,500))
lines((0:5000)/10,corm_fit,
      lwd=2,col="orange")
lines((0:5000)/10,double_sigmoid((0:5000)/10,0.013,0.12,70,0.09,150),
      lwd=2,col="darkblue")
lines((0:5000)/10,double_sigmoid((0:5000)/10,0.01165488,0.27474958,61.59843126,0.1,165),
      lwd=2,col="darkblue")

corm_param <- exp(opt$par)/exp(opt$par)*c(0.01165488,0.27474958,61.59843126,0.1,165)
  
corm_fit <- double_sigmoid((0:5000)/10,corm_param[1],corm_param[2],corm_param[3],corm_param[4],corm_param[5])


plot((0:5000)/10,corm_fit/max(corm_fit),
     lwd=3,col="darkblue",ylab = "Preference",
     xlab="length [mm]",main="Flatfish",type='l',ylim=c(0,1.5))
lines((0:5000)/10,seal_fit/max(seal_fit),
      lwd=3,col='red')
points(df.fin.corm$l.class,df.fin.corm$weighted.odds/max(corm_fit),pc=19,cex=1,col="darkblue")
points(df.fin.seal$l.class,df.fin.seal$weighted.odds/max(seal_fit),pc=19,cex=1,col="red")
legend(355,0.5,legend=c('Seal','Cormorant'),
       pch=19,col=c('darkblue','red'))

make_my_function <- function(corm_param) {
  maximum <- corm_param[1]   
  slope1 <- corm_param[2]    
  midPoint <- corm_param[3]      
  slope2 <- corm_param[4]    
  midPointDistance <- corm_param[5]   
  
  function(x) { # x is length in mm
    pref <- (maximum / (1 + exp(-slope1 * (x - midPoint)))) - 
      (maximum / (1 + exp(-slope2 * (x - (midPoint + midPointDistance)))))
    pref <- pref/max(pref)
    return(pref)
  }
}

# Create the function
corm_flatfish.pref <- make_my_function(corm_param)
#####
rm(list=setdiff(ls(),c('seal_cod.pref','seal_herring.pref','seal_flatfish.pref',
                       'corm_cod.pref','corm_herring.pref','corm_flatfish.pref')))


# Save and plot functions - NB save does not work
#####
x <- 0:700
plot(x,seal_cod.pref(x),col="red",lty=1,type = 'l',lwd=2)
lines(x,seal_herring.pref(x),col="red",lty=2,lwd=2)
lines(x,seal_flatfish.pref(x),col="red",lty=3,lwd=2)
lines(x,corm_cod.pref(x),col="darkblue",lty=1,lwd=2)
lines(x,corm_herring.pref(x),col="darkblue",lty=2,lwd=2)
lines(x,corm_flatfish.pref(x),col="darkblue",lty=3,lwd=2)

save(corm_cod.pref, file = "seal_cod_pref.RData")
save(corm_cod.pref, file = "seal_herring_pref.RData")
save(corm_cod.pref, file = "seal_flatfish_pref.RData")
save(corm_cod.pref, file = "corm_cod_pref.RData")
save(corm_cod.pref, file = "corm_herring_pref.RData")
save(corm_cod.pref, file = "corm_flatfish_pref.RData")
#####

