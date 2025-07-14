#run "size preference"

rm(list=setdiff(ls(),c('seal_cod.pref','seal_herring.pref','seal_flatfish.pref',
                       'corm_cod.pref','corm_herring.pref','corm_flatfish.pref')))

setwd("C:/Users/mavast/Desktop/Markus SOS/seal and cormorant size preference")
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
  Linf <- rnorm(n,mean=97.9,sd=19.99284) # Assymptotic length - with stochasticity
  #Linf <- 154.56 # cm
  L_hatch <- 0.4 # length at hatching, cm (Pepin et al, 1997)
  #t_hatch <- ((46+135)/2)/365
  t_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
  k <- rnorm(n,mean=0.22,sd=4.188023e-07) # vbgr growth rate
  u <- 0.23 # Amplitude of seasonality
  w <- 0.91 # Timing of peak growth
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  sd <- rep(0,length(t))
  for (i in 1:length(sd)){
    fish.age <- age[i]+rnorm(n,mean=0,sd=8.582605e-02) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
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

Cod <- read.table('data/M_sms.csv',sep=';',header=TRUE)
cormorant <- read.table('data/corm_diet.csv',sep=';',header=TRUE)
gseal <- read.table('data/gSeal_diet.csv',sep=';',header=TRUE)
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
  Linf <- rnorm(n,mean=30.57,sd=3.346919)
  k <- rnorm(n,0.453,5.649682e-07) # from table 2 in Gröhsler et al, 2013
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
    fish.age <- age[i]+rnorm(n,mean=0,sd=3.739677e-02) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
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
      if(yrs>=1991){
        di <- Herring %>% filter(ages==age[j] & years==yrs[i])
        dip1 <- Herring %>% filter(ages==(age[j]+1) & years==(yrs[i]+1))
        Ns[(365*(i-1)+1):(365*i),j] <- di$N*exp(-x*di$Z)
      }
      if(yrs<1991){
        di <- Herring %>% filter(ages==age[j] & years==1991)
        dip1 <- Herring %>% filter(ages==(age[j]+1) & years==(1991+1))
        Ns[(365*(i-1)+1):(365*i),j] <- di$N*exp(-x*di$Z)
      }
    }
  }
  return(Ns)
}
Herring <- read.table('data/Herring.csv',sep=';',header=TRUE)
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
  Linf <- rnorm(n,mean=34.0031653,sd=5.326541)
  k <- rnorm(n,0.4537312,1.575315e-06) # from table 2 in Gröhsler et al, 2013
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
    fish.age <- age[i]+rnorm(n,mean=0,sd=8.611682e-02) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
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
  Linf <- rnorm(n,mean=41.637260,sd=5.88133054 )
  k <- rnorm(n,0.240703,0.02414717) # from table 2 in Gröhsler et al, 2013
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
    fish.age <- age[i]+rnorm(n,mean=0,sd=0.06868120) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
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
  Linf <- rnorm(n,mean=26.5093927,sd=3.62925898  )
  k <- rnorm(n,0.4846326,0.15026563 ) # from table 2 in Gröhsler et al, 2013
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
    fish.age <- age[i]+rnorm(n,mean=0,sd=0.08027334 ) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
    L_t <- (L_hatch-Linf)*exp(-k*(phi_t[i]+fish.age-phi_hatch))+Linf # vbgr estimated length
    sd[i] <- sd(L_t)
  }
  return(sd)
}

Flatfish <- read.table("data/N_est.flatfish.csv",sep=';',header=TRUE)
Flatfish$yday <- round(Flatfish$jd*365)
names(Flatfish)[names(Flatfish)=="Species"] <- "species"

#Flatfish <- read.table("data/N_est.flatfish.csv",sep=';',header=TRUE)
#Flatfish$jd <- (Flatfish$jd*365+1)/365
#Flatfish$yday <- round(Flatfish$jd*365)
#####
setwd("C:/Users/mavast/Desktop/Markus SOS/Dirichlet diet estimation")
LW <- read.table("data/length-weight.csv",header=TRUE,sep=';')
month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
           "Sep","Oct","Nov","Dec") 
width <- 10
# Seal food index
#####
prey <- c("cod","herring","flounder","plaice","dab")
samplings <- unique(gseal$yms)
n_samplings <- length(samplings)
n_prey <- length(prey)


sealFood <- data.frame(sampling=rep("O",n_samplings*n_prey),
                       biomass=rep(0,n_samplings*n_prey),
                       species=rep("O",n_samplings*n_prey))

cod_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
herring_hatch <- mean(c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365)
Flounder_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Plaice_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Dab_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

for (i in 1:length(samplings)){
  # sampling specific
  dada <- gseal %>% filter(yms==samplings[i])
  n_sca <- length(unique(dada$Scat))
  date <- paste(dada$Year[1],'-',which(dada$Month[1]==month)[[1]],'-',15,sep="")
  jday <- yday(date)/365
  # cod
  ####
  cod.s <- c(vbgrCod(jday-cod_hatch),vbgrCod(jday-cod_hatch+1),vbgrCod(jday-cod_hatch+2),vbgrCod(jday-cod_hatch+3),
             vbgrCod(jday-cod_hatch+4),vbgrCod(jday-cod_hatch+5),vbgrCod(jday-cod_hatch+6),vbgrCod(jday-cod_hatch+7))*10
  cod.age <- c(0,1,2,3,4,5,6,7)
  cod.i <- N_ageCod(Cod,cod.age,dada$Year[1])[jday*365,1:length(cod.age)]
  cod.sd<- c(vbgr.sdCod(jday-cod_hatch,100000),vbgr.sdCod(jday-cod_hatch+1,100000),vbgr.sdCod(jday-cod_hatch+2,100000),vbgr.sdCod(jday-cod_hatch+3,100000),
         vbgr.sdCod(jday-cod_hatch+4,100000),vbgr.sdCod(jday-cod_hatch+5,100000),vbgr.sdCod(jday-cod_hatch+6,100000),vbgr.sdCod(jday-cod_hatch+7,100000))*10
  
  if (min(cod.sd)<0){
    cod.sd[cod.sd<0] <-rep(0,length(cod.sd[cod.sd<0]))
  }
  available.cod <- c(floor(rnorm(cod.i[1],cod.s[1],cod.sd[1])/width)*width+width/2,floor(rnorm(cod.i[2],cod.s[2],cod.sd[2])/width)*width+width/2,
                 floor(rnorm(cod.i[3],cod.s[3],cod.sd[3])/width)*width+width/2,floor(rnorm(cod.i[4],cod.s[4],cod.sd[4])/width)*width+width/2,
                 floor(rnorm(cod.i[5],cod.s[5],cod.sd[5])/width)*width+width/2,floor(rnorm(cod.i[6],cod.s[6],cod.sd[6])/width)*width+width/2,
                 floor(rnorm(cod.i[7],cod.s[7],cod.sd[7])/width)*width+width/2,floor(rnorm(cod.i[8],cod.s[8],cod.sd[8])/width)*width+width/2)
  available.cod[available.cod<0] <- rep(width/2,length(available.cod[available.cod<0]))
  
  cod.av <- as.data.frame(table(available.cod))
  cod.av$N <- as.numeric(cod.av$Freq)
  cod.av$l <- as.numeric(as.character(cod.av$available.cod))/10 # cm
  cod.coef <- LW %>% filter(species=="cod" & year==dada$Year[1])
  if(length(cod.coef$year)==0){
    cod.coef <- LW %>% filter(species=="cod" & year==min(year))
  } 
  cod.av$B <- cod.av$N*(cod.coef$a*cod.av$l^cod.coef$b) # g
  cod.av$avB <- pmax(0,seal_cod.pref(cod.av$l*10))*cod.av$B
  ####
  # herring
  ####
  herring.s <- c(vbgrHerring(jday-herring_hatch),vbgrHerring(jday-herring_hatch+1),
                 vbgrHerring(jday-herring_hatch+2),vbgrHerring(jday-herring_hatch+3),
                 vbgrHerring(jday-herring_hatch+4),vbgrHerring(jday-herring_hatch+5),
                 vbgrHerring(jday-herring_hatch+6),vbgrHerring(jday-herring_hatch+7),
                 vbgrHerring(jday-herring_hatch+8))*10
  herring.age <- 0:8
  herring.i <- N_ageHerring(herring,herring.age,dada$Year[1])[jday*365,1:length(herring.age)]
  herring.sd<- c(vbgr.sdHerring(jday-herring_hatch,100000),vbgr.sdHerring(jday-herring_hatch+1,100000),
                 vbgr.sdHerring(jday-herring_hatch+2,100000),vbgr.sdHerring(jday-herring_hatch+3,100000),
                 vbgr.sdHerring(jday-herring_hatch+4,100000),vbgr.sdHerring(jday-herring_hatch+5,100000),
                 vbgr.sdHerring(jday-herring_hatch+6,100000),vbgr.sdHerring(jday-herring_hatch+7,100000),
                 vbgr.sdHerring(jday-herring_hatch+8,100000))*10
  
  if (min(herring.sd)<0){
    herring.sd[herring.sd<0] <-rep(0,length(herring.sd[herring.sd<0]))
  }
  available.herring <- c(floor(rnorm(herring.i[1],herring.s[1],herring.sd[1])/width)*width+width/2,floor(rnorm(herring.i[2],herring.s[2],herring.sd[2])/width)*width+width/2,
                     floor(rnorm(herring.i[3],herring.s[3],herring.sd[3])/width)*width+width/2,floor(rnorm(herring.i[4],herring.s[4],herring.sd[4])/width)*width+width/2,
                     floor(rnorm(herring.i[5],herring.s[5],herring.sd[5])/width)*width+width/2,floor(rnorm(herring.i[6],herring.s[6],herring.sd[6])/width)*width+width/2,
                     floor(rnorm(herring.i[7],herring.s[7],herring.sd[7])/width)*width+width/2,floor(rnorm(herring.i[8],herring.s[8],herring.sd[8])/width)*width+width/2,
                     floor(rnorm(herring.i[9],herring.s[9],herring.sd[9])/width)*width+width/2)
  available.herring[available.herring<0] <- rep(width/2,length(available.herring[available.herring<0]))
  
  herring.av <- as.data.frame(table(available.herring))
  herring.av$N <- as.numeric(herring.av$Freq)
  herring.av$l <- as.numeric(as.character(herring.av$available.herring))/10 # cm
  herring.coef <- LW %>% filter(species=="herring")
  herring.av$B <- herring.av$N*(herring.coef$a*herring.av$l^herring.coef$b) # g
  herring.av$avB <- pmax(0,seal_herring.pref(herring.av$l*10))*herring.av$B
  ####
  
  # flatfish
  ####
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
  fl.i <- Flatfish %>% filter(year==dada$Year[1] & yday==round(jday*365))
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
  
  ####
  #flounder
  flounder.av <- as.data.frame(table(availableFlounder))
  flounder.av$N <- as.numeric(flounder.av$Freq)
  flounder.av$l <- as.numeric(as.character(flounder.av$availableFlounder))/10 # cm
  flounder.coef <- LW %>% filter(species=="flounder" & year==dada$Year[1])
  if(length(flounder.coef$year)==0){
    flounder.coef <- LW %>% filter(species=="flounder" & year==min(year))
  }  
  flounder.av$B <- flounder.av$N*(flounder.coef$a*flounder.av$l^flounder.coef$b) # g
  flounder.av$avB <- pmax(0,seal_flatfish.pref(flounder.av$l*10))*flounder.av$B
  
  #plaice
  plaice.av <- as.data.frame(table(availablePlaice))
  plaice.av$N <- as.numeric(plaice.av$Freq)
  plaice.av$l <- as.numeric(as.character(plaice.av$availablePlaice))/10 # cm
  plaice.coef <- LW %>% filter(species=="plaice" & year==dada$Year[1])
  if(length(plaice.coef$year)==0){
    plaice.coef <- LW %>% filter(species=="plaice" & year==min(year))
  }  
  plaice.av$B <- plaice.av$N*(plaice.coef$a*plaice.av$l^plaice.coef$b) # g
  plaice.av$avB <- pmax(0,seal_flatfish.pref(plaice.av$l*10))*plaice.av$B
  
  #dab
  dab.av <- as.data.frame(table(availableDab))
  dab.av$N <- as.numeric(dab.av$Freq)
  dab.av$l <- as.numeric(as.character(dab.av$availableDab))/10 # cm
  dab.coef <- LW %>% filter(species=="dab" & year==dada$Year[1])
  if(length(dab.coef$year)==0){
    dab.coef <- LW %>% filter(species=="dab" & year==min(year))
  }  
  dab.av$B <- dab.av$N*(dab.coef$a*dab.av$l^dab.coef$b) # g
  dab.av$avB <- pmax(0,seal_flatfish.pref(dab.av$l*10))*dab.av$B
  sealFood$sampling[((i-1)*n_prey+1):(i*n_prey)] <- samplings[i]
  sealFood$species[((i-1)*n_prey+1):(i*n_prey)] <- prey
  sealFood$biomass[(i-1)*n_prey+1] <- sum(cod.av$avB)
  sealFood$biomass[(i-1)*n_prey+2] <- sum(herring.av$avB)
  sealFood$biomass[(i-1)*n_prey+3] <- sum(flounder.av$avB)
  sealFood$biomass[(i-1)*n_prey+4] <- sum(plaice.av$avB)
  sealFood$biomass[(i-1)*n_prey+5] <- sum(dab.av$avB)
  
  print(i)
}

#####

# cormorant food index - quarter
#####
prey <- c("cod","herring","flounder","plaice","dab")
samplings <- unique(cormorant$dmyl)
n_samplings <- length(samplings)
n_prey <- length(prey)


cormorantFood <- data.frame(sampling=rep("O",n_samplings*n_prey),
                       biomass=rep(0,n_samplings*n_prey),
                       species=rep("O",n_samplings*n_prey))

cod_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
herring_hatch <- mean(c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365)
Flounder_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Plaice_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Dab_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

for (i in 1:length(samplings)){
  # sampling specific
  dada <- cormorant %>% filter(dmyl==samplings[i])
  date <- paste(dada$year[1],'-',which(dada$month[1]==month)[[1]],'-',dada$day[1],sep="")
  jday <- yday(date)/365
  # cod
  ####
  cod.s <- c(vbgrCod(jday-cod_hatch),vbgrCod(jday-cod_hatch+1),vbgrCod(jday-cod_hatch+2),vbgrCod(jday-cod_hatch+3),
             vbgrCod(jday-cod_hatch+4),vbgrCod(jday-cod_hatch+5),vbgrCod(jday-cod_hatch+6),vbgrCod(jday-cod_hatch+7))*10
  cod.age <- c(0,1,2,3,4,5,6,7)
  cod.i <- N_ageCod(Cod,cod.age,dada$year[1])[jday*365,1:length(cod.age)]
  cod.sd<- c(vbgr.sdCod(jday-cod_hatch,100000),vbgr.sdCod(jday-cod_hatch+1,100000),vbgr.sdCod(jday-cod_hatch+2,100000),vbgr.sdCod(jday-cod_hatch+3,100000),
             vbgr.sdCod(jday-cod_hatch+4,100000),vbgr.sdCod(jday-cod_hatch+5,100000),vbgr.sdCod(jday-cod_hatch+6,100000),vbgr.sdCod(jday-cod_hatch+7,100000))*10
  
  if (min(cod.sd)<0){
    cod.sd[cod.sd<0] <-rep(0,length(cod.sd[cod.sd<0]))
  }
  available.cod <- c(floor(rnorm(cod.i[1],cod.s[1],cod.sd[1])/width)*width+width/2,floor(rnorm(cod.i[2],cod.s[2],cod.sd[2])/width)*width+width/2,
                     floor(rnorm(cod.i[3],cod.s[3],cod.sd[3])/width)*width+width/2,floor(rnorm(cod.i[4],cod.s[4],cod.sd[4])/width)*width+width/2,
                     floor(rnorm(cod.i[5],cod.s[5],cod.sd[5])/width)*width+width/2,floor(rnorm(cod.i[6],cod.s[6],cod.sd[6])/width)*width+width/2,
                     floor(rnorm(cod.i[7],cod.s[7],cod.sd[7])/width)*width+width/2,floor(rnorm(cod.i[8],cod.s[8],cod.sd[8])/width)*width+width/2)
  available.cod[available.cod<0] <- rep(width/2,length(available.cod[available.cod<0]))
  
  cod.av <- as.data.frame(table(available.cod))
  cod.av$N <- as.numeric(cod.av$Freq)
  cod.av$l <- as.numeric(as.character(cod.av$available.cod))/10 # cm
  cod.coef <- LW %>% filter(species=="cod" & year==dada$year[1])
  if(length(cod.coef$year)==0){
    cod.coef <- LW %>% filter(species=="cod" & year==min(year))
  } 
  cod.av$B <- cod.av$N*(cod.coef$a*cod.av$l^cod.coef$b) # g
  cod.av$avB <- pmax(0,corm_cod.pref(cod.av$l*10))*cod.av$B
  ####
  # herring
  ####
  herring.s <- c(vbgrHerring(jday-herring_hatch),vbgrHerring(jday-herring_hatch+1),
                 vbgrHerring(jday-herring_hatch+2),vbgrHerring(jday-herring_hatch+3),
                 vbgrHerring(jday-herring_hatch+4),vbgrHerring(jday-herring_hatch+5),
                 vbgrHerring(jday-herring_hatch+6),vbgrHerring(jday-herring_hatch+7),
                 vbgrHerring(jday-herring_hatch+8))*10
  herring.age <- 0:8
  herring.i <- N_ageHerring(herring,herring.age,dada$year[1])[jday*365,1:length(herring.age)]
  herring.sd<- c(vbgr.sdHerring(jday-herring_hatch,100000),vbgr.sdHerring(jday-herring_hatch+1,100000),
                 vbgr.sdHerring(jday-herring_hatch+2,100000),vbgr.sdHerring(jday-herring_hatch+3,100000),
                 vbgr.sdHerring(jday-herring_hatch+4,100000),vbgr.sdHerring(jday-herring_hatch+5,100000),
                 vbgr.sdHerring(jday-herring_hatch+6,100000),vbgr.sdHerring(jday-herring_hatch+7,100000),
                 vbgr.sdHerring(jday-herring_hatch+8,100000))*10
  
  if (min(herring.sd)<0){
    herring.sd[herring.sd<0] <-rep(0,length(herring.sd[herring.sd<0]))
  }
  available.herring <- c(floor(rnorm(herring.i[1],herring.s[1],herring.sd[1])/width)*width+width/2,floor(rnorm(herring.i[2],herring.s[2],herring.sd[2])/width)*width+width/2,
                         floor(rnorm(herring.i[3],herring.s[3],herring.sd[3])/width)*width+width/2,floor(rnorm(herring.i[4],herring.s[4],herring.sd[4])/width)*width+width/2,
                         floor(rnorm(herring.i[5],herring.s[5],herring.sd[5])/width)*width+width/2,floor(rnorm(herring.i[6],herring.s[6],herring.sd[6])/width)*width+width/2,
                         floor(rnorm(herring.i[7],herring.s[7],herring.sd[7])/width)*width+width/2,floor(rnorm(herring.i[8],herring.s[8],herring.sd[8])/width)*width+width/2,
                         floor(rnorm(herring.i[9],herring.s[9],herring.sd[9])/width)*width+width/2)
  available.herring[available.herring<0] <- rep(width/2,length(available.herring[available.herring<0]))
  
  herring.av <- as.data.frame(table(available.herring))
  herring.av$N <- as.numeric(herring.av$Freq)
  herring.av$l <- as.numeric(as.character(herring.av$available.herring))/10 # cm
  herring.coef <- LW %>% filter(species=="herring")
  herring.av$B <- herring.av$N*(herring.coef$a*herring.av$l^herring.coef$b) # g
  herring.av$avB <- pmax(0,corm_herring.pref(herring.av$l*10))*herring.av$B
  ####
  
  # flatfish
  ####
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
  fl.i <- Flatfish %>% filter(year==dada$year[1] & yday==round(jday*365))
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
  
  ####
  #flounder
  flounder.av <- as.data.frame(table(availableFlounder))
  flounder.av$N <- as.numeric(flounder.av$Freq)
  flounder.av$l <- as.numeric(as.character(flounder.av$availableFlounder))/10 # cm
  flounder.coef <- LW %>% filter(species=="flounder" & year==dada$year[1])
  if(length(flounder.coef$year)==0){
    flounder.coef <- LW %>% filter(species=="flounder" & year==min(year))
  }  
  flounder.av$B <- flounder.av$N*(flounder.coef$a*flounder.av$l^flounder.coef$b) # g
  flounder.av$avB <- pmax(0,corm_flatfish.pref(flounder.av$l*10))*flounder.av$B
  
  #plaice
  plaice.av <- as.data.frame(table(availablePlaice))
  plaice.av$N <- as.numeric(plaice.av$Freq)
  plaice.av$l <- as.numeric(as.character(plaice.av$availablePlaice))/10 # cm
  plaice.coef <- LW %>% filter(species=="plaice" & year==dada$year[1])
  if(length(plaice.coef$year)==0){
    LW.plaice <- LW %>% filter(species=="plaice")
    plaice.coef <- LW.plaice %>% filter(year==min(year))
  }  
  plaice.av$B <- plaice.av$N*(plaice.coef$a*plaice.av$l^plaice.coef$b) # g
  plaice.av$avB <- pmax(0,corm_flatfish.pref(plaice.av$l*10))*plaice.av$B
  
  #dab
  dab.av <- as.data.frame(table(availableDab))
  dab.av$N <- as.numeric(dab.av$Freq)
  dab.av$l <- as.numeric(as.character(dab.av$availableDab))/10 # cm
  dab.coef <- LW %>% filter(species=="dab" & year==dada$year[1])
  if(length(dab.coef$year)==0){
    LW.dab <- LW %>% filter(species=="dab")
    dab.coef <- LW.dab %>% filter(species=="dab" & year==min(year))
  }  
  dab.av$B <- dab.av$N*(dab.coef$a*dab.av$l^dab.coef$b) # g
  dab.av$avB <- pmax(0,corm_flatfish.pref(dab.av$l*10))*dab.av$B
  cormorantFood$sampling[((i-1)*n_prey+1):(i*n_prey)] <- samplings[i]
  cormorantFood$species[((i-1)*n_prey+1):(i*n_prey)] <- prey
  cormorantFood$biomass[(i-1)*n_prey+1] <- sum(cod.av$avB)
  cormorantFood$biomass[(i-1)*n_prey+2] <- sum(herring.av$avB)
  cormorantFood$biomass[(i-1)*n_prey+3] <- sum(flounder.av$avB)
  cormorantFood$biomass[(i-1)*n_prey+4] <- sum(plaice.av$avB)
  cormorantFood$biomass[(i-1)*n_prey+5] <- sum(dab.av$avB)
  
  print(i)
}

#####

# cormorant food index - month
#####
prey <- c("cod","herring","flounder","plaice","dab")
months <- c(Jan=1,Feb=2,Mar=3,Apr=4,May=5,Jun=6,Jul=7,Aug=8,Sep=9,Oct=10,Nov=11,Dec=12)
samplings <- unique(cd$ym)  # get cd from DirichletDiet "Prepare cormorant data - month"
n_samplings <- length(samplings)
n_prey <- length(prey)


cormorantFood <- data.frame(sampling=rep("O",n_samplings*n_prey),
                            biomass=rep(0,n_samplings*n_prey),
                            species=rep("O",n_samplings*n_prey))

cod_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
herring_hatch <- mean(c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365)
Flounder_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Plaice_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Dab_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

for (i in 1:length(samplings)){
  # sampling specific
  year.i <- as.numeric(substr(samplings[i],1,4))
  month.i <- substr(samplings[i],6,8)
  date <- paste(year.i,'-',months[month.i],'-',15,sep="")
  jday <- yday(date)/365
  # cod
  ####
  cod.s <- c(vbgrCod(jday-cod_hatch),vbgrCod(jday-cod_hatch+1),vbgrCod(jday-cod_hatch+2),vbgrCod(jday-cod_hatch+3),
             vbgrCod(jday-cod_hatch+4),vbgrCod(jday-cod_hatch+5),vbgrCod(jday-cod_hatch+6),vbgrCod(jday-cod_hatch+7))*10
  cod.age <- c(0,1,2,3,4,5,6,7)
  cod.i <- N_ageCod(Cod,cod.age,year.i)[jday*365,1:length(cod.age)]
  cod.sd<- c(vbgr.sdCod(jday-cod_hatch,100000),vbgr.sdCod(jday-cod_hatch+1,100000),vbgr.sdCod(jday-cod_hatch+2,100000),vbgr.sdCod(jday-cod_hatch+3,100000),
             vbgr.sdCod(jday-cod_hatch+4,100000),vbgr.sdCod(jday-cod_hatch+5,100000),vbgr.sdCod(jday-cod_hatch+6,100000),vbgr.sdCod(jday-cod_hatch+7,100000))*10
  
  if (min(cod.sd)<0){
    cod.sd[cod.sd<0] <-rep(0,length(cod.sd[cod.sd<0]))
  }
  available.cod <- c(floor(rnorm(cod.i[1],cod.s[1],cod.sd[1])/width)*width+width/2,floor(rnorm(cod.i[2],cod.s[2],cod.sd[2])/width)*width+width/2,
                     floor(rnorm(cod.i[3],cod.s[3],cod.sd[3])/width)*width+width/2,floor(rnorm(cod.i[4],cod.s[4],cod.sd[4])/width)*width+width/2,
                     floor(rnorm(cod.i[5],cod.s[5],cod.sd[5])/width)*width+width/2,floor(rnorm(cod.i[6],cod.s[6],cod.sd[6])/width)*width+width/2,
                     floor(rnorm(cod.i[7],cod.s[7],cod.sd[7])/width)*width+width/2,floor(rnorm(cod.i[8],cod.s[8],cod.sd[8])/width)*width+width/2)
  available.cod[available.cod<0] <- rep(width/2,length(available.cod[available.cod<0]))
  
  cod.av <- as.data.frame(table(available.cod))
  cod.av$N <- as.numeric(cod.av$Freq)
  cod.av$l <- as.numeric(as.character(cod.av$available.cod))/10 # cm
  cod.coef <- LW %>% filter(species=="cod" & year==year.i)
  if(length(cod.coef$year)==0){
    cod.coef <- LW %>% filter(species=="cod" & year==min(year))
  } 
  cod.av$B <- cod.av$N*(cod.coef$a*cod.av$l^cod.coef$b) # g
  cod.av$avB <- pmax(0,corm_cod.pref(cod.av$l*10))*cod.av$B
  ####
  # herring
  ####
  herring.s <- c(vbgrHerring(jday-herring_hatch),vbgrHerring(jday-herring_hatch+1),
                 vbgrHerring(jday-herring_hatch+2),vbgrHerring(jday-herring_hatch+3),
                 vbgrHerring(jday-herring_hatch+4),vbgrHerring(jday-herring_hatch+5),
                 vbgrHerring(jday-herring_hatch+6),vbgrHerring(jday-herring_hatch+7),
                 vbgrHerring(jday-herring_hatch+8))*10
  herring.age <- 0:8
  herring.i <- N_ageHerring(herring,herring.age,year.i)[jday*365,1:length(herring.age)]
  herring.sd<- c(vbgr.sdHerring(jday-herring_hatch,100000),vbgr.sdHerring(jday-herring_hatch+1,100000),
                 vbgr.sdHerring(jday-herring_hatch+2,100000),vbgr.sdHerring(jday-herring_hatch+3,100000),
                 vbgr.sdHerring(jday-herring_hatch+4,100000),vbgr.sdHerring(jday-herring_hatch+5,100000),
                 vbgr.sdHerring(jday-herring_hatch+6,100000),vbgr.sdHerring(jday-herring_hatch+7,100000),
                 vbgr.sdHerring(jday-herring_hatch+8,100000))*10
  
  if (min(herring.sd)<0){
    herring.sd[herring.sd<0] <-rep(0,length(herring.sd[herring.sd<0]))
  }
  available.herring <- c(floor(rnorm(herring.i[1],herring.s[1],herring.sd[1])/width)*width+width/2,floor(rnorm(herring.i[2],herring.s[2],herring.sd[2])/width)*width+width/2,
                         floor(rnorm(herring.i[3],herring.s[3],herring.sd[3])/width)*width+width/2,floor(rnorm(herring.i[4],herring.s[4],herring.sd[4])/width)*width+width/2,
                         floor(rnorm(herring.i[5],herring.s[5],herring.sd[5])/width)*width+width/2,floor(rnorm(herring.i[6],herring.s[6],herring.sd[6])/width)*width+width/2,
                         floor(rnorm(herring.i[7],herring.s[7],herring.sd[7])/width)*width+width/2,floor(rnorm(herring.i[8],herring.s[8],herring.sd[8])/width)*width+width/2,
                         floor(rnorm(herring.i[9],herring.s[9],herring.sd[9])/width)*width+width/2)
  available.herring[available.herring<0] <- rep(width/2,length(available.herring[available.herring<0]))
  
  herring.av <- as.data.frame(table(available.herring))
  herring.av$N <- as.numeric(herring.av$Freq)
  herring.av$l <- as.numeric(as.character(herring.av$available.herring))/10 # cm
  herring.coef <- LW %>% filter(species=="herring")
  herring.av$B <- herring.av$N*(herring.coef$a*herring.av$l^herring.coef$b) # g
  herring.av$avB <- pmax(0,corm_herring.pref(herring.av$l*10))*herring.av$B
  ####
  
  # flatfish
  ####
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
  fl.i <- Flatfish %>% filter(year==year.i & yday==round(jday*365))
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
  
  ####
  #flounder
  flounder.av <- as.data.frame(table(availableFlounder))
  flounder.av$N <- as.numeric(flounder.av$Freq)
  flounder.av$l <- as.numeric(as.character(flounder.av$availableFlounder))/10 # cm
  flounder.coef <- LW %>% filter(species=="flounder" & year==year.i)
  if(length(flounder.coef$year)==0){
    flounder.coef <- LW %>% filter(species=="flounder" & year==min(year))
  }  
  flounder.av$B <- flounder.av$N*(flounder.coef$a*flounder.av$l^flounder.coef$b) # g
  flounder.av$avB <- pmax(0,corm_flatfish.pref(flounder.av$l*10))*flounder.av$B
  
  #plaice
  plaice.av <- as.data.frame(table(availablePlaice))
  plaice.av$N <- as.numeric(plaice.av$Freq)
  plaice.av$l <- as.numeric(as.character(plaice.av$availablePlaice))/10 # cm
  plaice.coef <- LW %>% filter(species=="plaice" & year==year.i)
  if(length(plaice.coef$year)==0){
    LW.plaice <- LW %>% filter(species=="plaice")
    plaice.coef <- LW.plaice %>% filter(year==min(year))
  }  
  plaice.av$B <- plaice.av$N*(plaice.coef$a*plaice.av$l^plaice.coef$b) # g
  plaice.av$avB <- pmax(0,corm_flatfish.pref(plaice.av$l*10))*plaice.av$B
  
  #dab
  dab.av <- as.data.frame(table(availableDab))
  dab.av$N <- as.numeric(dab.av$Freq)
  dab.av$l <- as.numeric(as.character(dab.av$availableDab))/10 # cm
  dab.coef <- LW %>% filter(species=="dab" & year==year.i)
  if(length(dab.coef$year)==0){
    LW.dab <- LW %>% filter(species=="dab")
    dab.coef <- LW.dab %>% filter(species=="dab" & year==min(year))
  }  
  dab.av$B <- dab.av$N*(dab.coef$a*dab.av$l^dab.coef$b) # g
  dab.av$avB <- pmax(0,corm_flatfish.pref(dab.av$l*10))*dab.av$B
  cormorantFood$sampling[((i-1)*n_prey+1):(i*n_prey)] <- samplings[i]
  cormorantFood$species[((i-1)*n_prey+1):(i*n_prey)] <- prey
  cormorantFood$biomass[(i-1)*n_prey+1] <- sum(cod.av$avB)
  cormorantFood$biomass[(i-1)*n_prey+2] <- sum(herring.av$avB)
  cormorantFood$biomass[(i-1)*n_prey+3] <- sum(flounder.av$avB)
  cormorantFood$biomass[(i-1)*n_prey+4] <- sum(plaice.av$avB)
  cormorantFood$biomass[(i-1)*n_prey+5] <- sum(dab.av$avB)
  
  print(i)
}
write.table(cormorantFood,"cormorantFood_month.csv",row.names = FALSE,sep=';')

#####

# seal food index - month
#####
prey <- c("cod","herring","flounder","plaice","dab")
months <- c(Jan=1,Feb=2,Mar=3,Apr=4,May=5,Jun=6,Jul=7,Aug=8,Sep=9,Oct=10,Nov=11,Dec=12)
samplings <- unique(sd$ym)  # get cd from DirichletDiet "Prepare seal data - month"
n_samplings <- length(samplings)
n_prey <- length(prey)


sealFood <- data.frame(sampling=rep("O",n_samplings*n_prey),
                            biomass=rep(0,n_samplings*n_prey),
                            species=rep("O",n_samplings*n_prey))

cod_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
herring_hatch <- mean(c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365)
Flounder_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Plaice_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Dab_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

for (i in 1:length(samplings)){
  # sampling specific
  year.i <- as.numeric(substr(samplings[i],1,4))
  month.i <- substr(samplings[i],6,8)
  date <- paste(year.i,'-',months[month.i],'-',15,sep="")
  jday <- yday(date)/365
  # cod
  ####
  cod.s <- c(vbgrCod(jday-cod_hatch),vbgrCod(jday-cod_hatch+1),vbgrCod(jday-cod_hatch+2),vbgrCod(jday-cod_hatch+3),
             vbgrCod(jday-cod_hatch+4),vbgrCod(jday-cod_hatch+5),vbgrCod(jday-cod_hatch+6),vbgrCod(jday-cod_hatch+7))*10
  cod.age <- c(0,1,2,3,4,5,6,7)
  cod.i <- N_ageCod(Cod,cod.age,year.i)[jday*365,1:length(cod.age)]
  cod.sd<- c(vbgr.sdCod(jday-cod_hatch,100000),vbgr.sdCod(jday-cod_hatch+1,100000),vbgr.sdCod(jday-cod_hatch+2,100000),vbgr.sdCod(jday-cod_hatch+3,100000),
             vbgr.sdCod(jday-cod_hatch+4,100000),vbgr.sdCod(jday-cod_hatch+5,100000),vbgr.sdCod(jday-cod_hatch+6,100000),vbgr.sdCod(jday-cod_hatch+7,100000))*10
  
  if (min(cod.sd)<0){
    cod.sd[cod.sd<0] <-rep(0,length(cod.sd[cod.sd<0]))
  }
  available.cod <- c(floor(rnorm(cod.i[1],cod.s[1],cod.sd[1])/width)*width+width/2,floor(rnorm(cod.i[2],cod.s[2],cod.sd[2])/width)*width+width/2,
                     floor(rnorm(cod.i[3],cod.s[3],cod.sd[3])/width)*width+width/2,floor(rnorm(cod.i[4],cod.s[4],cod.sd[4])/width)*width+width/2,
                     floor(rnorm(cod.i[5],cod.s[5],cod.sd[5])/width)*width+width/2,floor(rnorm(cod.i[6],cod.s[6],cod.sd[6])/width)*width+width/2,
                     floor(rnorm(cod.i[7],cod.s[7],cod.sd[7])/width)*width+width/2,floor(rnorm(cod.i[8],cod.s[8],cod.sd[8])/width)*width+width/2)
  available.cod[available.cod<0] <- rep(width/2,length(available.cod[available.cod<0]))
  
  cod.av <- as.data.frame(table(available.cod))
  cod.av$N <- as.numeric(cod.av$Freq)
  cod.av$l <- as.numeric(as.character(cod.av$available.cod))/10 # cm
  cod.coef <- LW %>% filter(species=="cod" & year==year.i)
  if(length(cod.coef$year)==0){
    cod.coef <- LW %>% filter(species=="cod" & year==min(year))
  } 
  cod.av$B <- cod.av$N*(cod.coef$a*cod.av$l^cod.coef$b) # g
  cod.av$avB <- pmax(0,seal_cod.pref(cod.av$l*10))*cod.av$B
  ####
  # herring
  ####
  herring.s <- c(vbgrHerring(jday-herring_hatch),vbgrHerring(jday-herring_hatch+1),
                 vbgrHerring(jday-herring_hatch+2),vbgrHerring(jday-herring_hatch+3),
                 vbgrHerring(jday-herring_hatch+4),vbgrHerring(jday-herring_hatch+5),
                 vbgrHerring(jday-herring_hatch+6),vbgrHerring(jday-herring_hatch+7),
                 vbgrHerring(jday-herring_hatch+8))*10
  herring.age <- 0:8
  herring.i <- N_ageHerring(herring,herring.age,year.i)[jday*365,1:length(herring.age)]
  herring.sd<- c(vbgr.sdHerring(jday-herring_hatch,100000),vbgr.sdHerring(jday-herring_hatch+1,100000),
                 vbgr.sdHerring(jday-herring_hatch+2,100000),vbgr.sdHerring(jday-herring_hatch+3,100000),
                 vbgr.sdHerring(jday-herring_hatch+4,100000),vbgr.sdHerring(jday-herring_hatch+5,100000),
                 vbgr.sdHerring(jday-herring_hatch+6,100000),vbgr.sdHerring(jday-herring_hatch+7,100000),
                 vbgr.sdHerring(jday-herring_hatch+8,100000))*10
  
  if (min(herring.sd)<0){
    herring.sd[herring.sd<0] <-rep(0,length(herring.sd[herring.sd<0]))
  }
  available.herring <- c(floor(rnorm(herring.i[1],herring.s[1],herring.sd[1])/width)*width+width/2,floor(rnorm(herring.i[2],herring.s[2],herring.sd[2])/width)*width+width/2,
                         floor(rnorm(herring.i[3],herring.s[3],herring.sd[3])/width)*width+width/2,floor(rnorm(herring.i[4],herring.s[4],herring.sd[4])/width)*width+width/2,
                         floor(rnorm(herring.i[5],herring.s[5],herring.sd[5])/width)*width+width/2,floor(rnorm(herring.i[6],herring.s[6],herring.sd[6])/width)*width+width/2,
                         floor(rnorm(herring.i[7],herring.s[7],herring.sd[7])/width)*width+width/2,floor(rnorm(herring.i[8],herring.s[8],herring.sd[8])/width)*width+width/2,
                         floor(rnorm(herring.i[9],herring.s[9],herring.sd[9])/width)*width+width/2)
  available.herring[available.herring<0] <- rep(width/2,length(available.herring[available.herring<0]))
  
  herring.av <- as.data.frame(table(available.herring))
  herring.av$N <- as.numeric(herring.av$Freq)
  herring.av$l <- as.numeric(as.character(herring.av$available.herring))/10 # cm
  herring.coef <- LW %>% filter(species=="herring")
  herring.av$B <- herring.av$N*(herring.coef$a*herring.av$l^herring.coef$b) # g
  herring.av$avB <- pmax(0,seal_herring.pref(herring.av$l*10))*herring.av$B
  ####
  
  # flatfish
  ####
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
  fl.i <- Flatfish %>% filter(year==year.i & yday==round(jday*365))
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
  
  ####
  #flounder
  flounder.av <- as.data.frame(table(availableFlounder))
  flounder.av$N <- as.numeric(flounder.av$Freq)
  flounder.av$l <- as.numeric(as.character(flounder.av$availableFlounder))/10 # cm
  flounder.coef <- LW %>% filter(species=="flounder" & year==year.i)
  if(length(flounder.coef$year)==0){
    flounder.coef <- LW %>% filter(species=="flounder" & year==min(year))
  }  
  flounder.av$B <- flounder.av$N*(flounder.coef$a*flounder.av$l^flounder.coef$b) # g
  flounder.av$avB <- pmax(0,seal_flatfish.pref(flounder.av$l*10))*flounder.av$B
  
  #plaice
  plaice.av <- as.data.frame(table(availablePlaice))
  plaice.av$N <- as.numeric(plaice.av$Freq)
  plaice.av$l <- as.numeric(as.character(plaice.av$availablePlaice))/10 # cm
  plaice.coef <- LW %>% filter(species=="plaice" & year==year.i)
  if(length(plaice.coef$year)==0){
    LW.plaice <- LW %>% filter(species=="plaice")
    plaice.coef <- LW.plaice %>% filter(year==min(year))
  }  
  plaice.av$B <- plaice.av$N*(plaice.coef$a*plaice.av$l^plaice.coef$b) # g
  plaice.av$avB <- pmax(0,seal_flatfish.pref(plaice.av$l*10))*plaice.av$B
  
  #dab
  dab.av <- as.data.frame(table(availableDab))
  dab.av$N <- as.numeric(dab.av$Freq)
  dab.av$l <- as.numeric(as.character(dab.av$availableDab))/10 # cm
  dab.coef <- LW %>% filter(species=="dab" & year==year.i)
  if(length(dab.coef$year)==0){
    LW.dab <- LW %>% filter(species=="dab")
    dab.coef <- LW.dab %>% filter(species=="dab" & year==min(year))
  }  
  dab.av$B <- dab.av$N*(dab.coef$a*dab.av$l^dab.coef$b) # g
  dab.av$avB <- pmax(0,seal_flatfish.pref(dab.av$l*10))*dab.av$B
  sealFood$sampling[((i-1)*n_prey+1):(i*n_prey)] <- samplings[i]
  sealFood$species[((i-1)*n_prey+1):(i*n_prey)] <- prey
  sealFood$biomass[(i-1)*n_prey+1] <- sum(cod.av$avB)
  sealFood$biomass[(i-1)*n_prey+2] <- sum(herring.av$avB)
  sealFood$biomass[(i-1)*n_prey+3] <- sum(flounder.av$avB)
  sealFood$biomass[(i-1)*n_prey+4] <- sum(plaice.av$avB)
  sealFood$biomass[(i-1)*n_prey+5] <- sum(dab.av$avB)
  
  print(i)
}
write.table(sealFood,"sealFood_month.csv",row.names = FALSE,sep=';')

#####

# normalize prey biomass
#####
cormorantFood$normB <- 0
sealFood$normB <- 0
for(i in 1:length(sealFood$normB)){
  ms <- mean(sealFood$biomass[sealFood$species==sealFood$species[i]])
  sealFood$normB[i] <- sealFood$biomass[i]/ms
}
for(i in 1:length(cormorantFood$normB)){
  mc <- mean(cormorantFood$biomass[cormorantFood$species==cormorantFood$species[i]])
  cormorantFood$normB[i] <- cormorantFood$biomass[i]/mc
}
#####
write.table(cormorantFood,"cormorantFood.csv",row.names = FALSE,sep=';')
write.table(sealFood,"sealFood.csv",row.names = FALSE,sep=';')


# simulate seal food index quarter 2 1991-2023
#####
prey <- c("cod","herring","flounder","plaice","dab")
samplings <- paste(rep(1991:2023),2)
n_samplings <- length(samplings)
n_prey <- length(prey)


sealFood <- data.frame(sampling=rep("O",n_samplings*n_prey),
                       biomass=rep(0,n_samplings*n_prey),
                       species=rep("O",n_samplings*n_prey))

cod_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
herring_hatch <- mean(c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365)
Flounder_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Plaice_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Dab_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

for (i in 1:length(samplings)){
  # sampling specific
  quarter <- as.numeric(substr(samplings[i],6,7))
  jday <- (c(1/4,2/4,3/4,4/4)-1/8)[quarter]
  year.i <- as.numeric(substr(samplings[i],1,4))
  # cod
  ####
  cod.s <- c(vbgrCod(jday-cod_hatch),vbgrCod(jday-cod_hatch+1),vbgrCod(jday-cod_hatch+2),vbgrCod(jday-cod_hatch+3),
             vbgrCod(jday-cod_hatch+4),vbgrCod(jday-cod_hatch+5),vbgrCod(jday-cod_hatch+6),vbgrCod(jday-cod_hatch+7))*10
  cod.age <- 0:7
  cod.i <- N_ageCod(Cod,cod.age,year.i)[jday*365,1:length(cod.age)]
  cod.sd<- c(vbgr.sdCod(jday-cod_hatch,100000),vbgr.sdCod(jday-cod_hatch+1,100000),vbgr.sdCod(jday-cod_hatch+2,100000),vbgr.sdCod(jday-cod_hatch+3,100000),
             vbgr.sdCod(jday-cod_hatch+4,100000),vbgr.sdCod(jday-cod_hatch+5,100000),vbgr.sdCod(jday-cod_hatch+6,100000),vbgr.sdCod(jday-cod_hatch+7,100000))*10
  
  if (min(cod.sd)<0){
    cod.sd[cod.sd<0] <-rep(0,length(cod.sd[cod.sd<0]))
  }
  available.cod <- c(floor(rnorm(cod.i[1],cod.s[1],cod.sd[1])/width)*width+width/2,floor(rnorm(cod.i[2],cod.s[2],cod.sd[2])/width)*width+width/2,
                     floor(rnorm(cod.i[3],cod.s[3],cod.sd[3])/width)*width+width/2,floor(rnorm(cod.i[4],cod.s[4],cod.sd[4])/width)*width+width/2,
                     floor(rnorm(cod.i[5],cod.s[5],cod.sd[5])/width)*width+width/2,floor(rnorm(cod.i[6],cod.s[6],cod.sd[6])/width)*width+width/2,
                     floor(rnorm(cod.i[7],cod.s[7],cod.sd[7])/width)*width+width/2,floor(rnorm(cod.i[8],cod.s[8],cod.sd[8])/width)*width+width/2)
  available.cod[available.cod<0] <- rep(width/2,length(available.cod[available.cod<0]))
  
  cod.av <- as.data.frame(table(available.cod))
  cod.av$N <- as.numeric(cod.av$Freq)
  cod.av$l <- as.numeric(as.character(cod.av$available.cod))/10 # cm
  cod.coef <- LW %>% filter(species=="cod" & year==year.i)
  if(length(cod.coef$year)==0){
    cod.coef <- LW %>% filter(species=="cod" & year==min(year))
  } 
  cod.av$B <- cod.av$N*(cod.coef$a*cod.av$l^cod.coef$b) # g
  cod.av$avB <- pmax(0,seal_cod.pref(cod.av$l*10))*cod.av$B
  ####
  # herring
  ####
  herring.s <- c(vbgrHerring(jday-herring_hatch),vbgrHerring(jday-herring_hatch+1),
                 vbgrHerring(jday-herring_hatch+2),vbgrHerring(jday-herring_hatch+3),
                 vbgrHerring(jday-herring_hatch+4),vbgrHerring(jday-herring_hatch+5),
                 vbgrHerring(jday-herring_hatch+6),vbgrHerring(jday-herring_hatch+7),
                 vbgrHerring(jday-herring_hatch+8))*10
  herring.age <- 0:8
  herring.i <- N_ageHerring(herring,herring.age,year.i)[jday*365,1:length(herring.age)]
  herring.sd<- c(vbgr.sdHerring(jday-herring_hatch,100000),vbgr.sdHerring(jday-herring_hatch+1,100000),
                 vbgr.sdHerring(jday-herring_hatch+2,100000),vbgr.sdHerring(jday-herring_hatch+3,100000),
                 vbgr.sdHerring(jday-herring_hatch+4,100000),vbgr.sdHerring(jday-herring_hatch+5,100000),
                 vbgr.sdHerring(jday-herring_hatch+6,100000),vbgr.sdHerring(jday-herring_hatch+7,100000),
                 vbgr.sdHerring(jday-herring_hatch+8,100000))*10
  
  if (min(herring.sd)<0){
    herring.sd[herring.sd<0] <-rep(0,length(herring.sd[herring.sd<0]))
  }
  available.herring <- c(floor(rnorm(herring.i[1],herring.s[1],herring.sd[1])/width)*width+width/2,floor(rnorm(herring.i[2],herring.s[2],herring.sd[2])/width)*width+width/2,
                         floor(rnorm(herring.i[3],herring.s[3],herring.sd[3])/width)*width+width/2,floor(rnorm(herring.i[4],herring.s[4],herring.sd[4])/width)*width+width/2,
                         floor(rnorm(herring.i[5],herring.s[5],herring.sd[5])/width)*width+width/2,floor(rnorm(herring.i[6],herring.s[6],herring.sd[6])/width)*width+width/2,
                         floor(rnorm(herring.i[7],herring.s[7],herring.sd[7])/width)*width+width/2,floor(rnorm(herring.i[8],herring.s[8],herring.sd[8])/width)*width+width/2,
                         floor(rnorm(herring.i[9],herring.s[9],herring.sd[9])/width)*width+width/2)
  available.herring[available.herring<0] <- rep(width/2,length(available.herring[available.herring<0]))
  
  herring.av <- as.data.frame(table(available.herring))
  herring.av$N <- as.numeric(herring.av$Freq)
  herring.av$l <- as.numeric(as.character(herring.av$available.herring))/10 # cm
  herring.coef <- LW %>% filter(species=="herring")
  herring.av$B <- herring.av$N*(herring.coef$a*herring.av$l^herring.coef$b) # g
  herring.av$avB <- pmax(0,seal_herring.pref(herring.av$l*10))*herring.av$B
  ####
  
  # flatfish
  ####
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
  fl.i <- Flatfish %>% filter(year==year.i & yday==round(jday*365))
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
  
  ####
  #flounder
  flounder.av <- as.data.frame(table(availableFlounder))
  flounder.av$N <- as.numeric(flounder.av$Freq)
  flounder.av$l <- as.numeric(as.character(flounder.av$availableFlounder))/10 # cm
  flounder.coef <- LW %>% filter(species=="flounder" & year==year.i)
  if(length(flounder.coef$year)==0){
    flounder.coef <- LW %>% filter(species=="flounder" & year==1992)
  }  
  flounder.av$B <- flounder.av$N*(flounder.coef$a*flounder.av$l^flounder.coef$b) # g
  flounder.av$avB <- pmax(0,seal_flatfish.pref(flounder.av$l*10))*flounder.av$B
  
  #plaice
  plaice.av <- as.data.frame(table(availablePlaice))
  plaice.av$N <- as.numeric(plaice.av$Freq)
  plaice.av$l <- as.numeric(as.character(plaice.av$availablePlaice))/10 # cm
  plaice.coef <- LW %>% filter(species=="plaice" & year==year.i)
  if(length(plaice.coef$year)==0){
    plaice.coef <- LW %>% filter(species=="plaice" & year==1994)
  }  
  plaice.av$B <- plaice.av$N*(plaice.coef$a*plaice.av$l^plaice.coef$b) # g
  plaice.av$avB <- pmax(0,seal_flatfish.pref(plaice.av$l*10))*plaice.av$B
  
  #dab
  dab.av <- as.data.frame(table(availableDab))
  dab.av$N <- as.numeric(dab.av$Freq)
  dab.av$l <- as.numeric(as.character(dab.av$availableDab))/10 # cm
  dab.coef <- LW %>% filter(species=="dab" & year==year.i)
  if(length(dab.coef$year)==0){
    dab.coef <- LW %>% filter(species=="dab" & year==1994)
  }  
  dab.av$B <- dab.av$N*(dab.coef$a*dab.av$l^dab.coef$b) # g
  dab.av$avB <- pmax(0,seal_flatfish.pref(dab.av$l*10))*dab.av$B
  sealFood$sampling[((i-1)*n_prey+1):(i*n_prey)] <- samplings[i]
  sealFood$species[((i-1)*n_prey+1):(i*n_prey)] <- prey
  sealFood$biomass[(i-1)*n_prey+1] <- sum(cod.av$avB)
  sealFood$biomass[(i-1)*n_prey+2] <- sum(herring.av$avB)
  sealFood$biomass[(i-1)*n_prey+3] <- sum(flounder.av$avB)
  sealFood$biomass[(i-1)*n_prey+4] <- sum(plaice.av$avB)
  sealFood$biomass[(i-1)*n_prey+5] <- sum(dab.av$avB)
  
  print(i)
}

write.table(sealFood,"Sealfood_sim.csv",row.names = FALSE,sep=';')
#####

# simulate cormorant food index quarter 2 1991-2023
#####
prey <- c("cod","herring","flounder","plaice","dab")
samplings <- paste(rep(1991:2023),2)
n_samplings <- length(samplings)
n_prey <- length(prey)


cormFood <- data.frame(sampling=rep("O",n_samplings*n_prey),
                       biomass=rep(0,n_samplings*n_prey),
                       species=rep("O",n_samplings*n_prey))

cod_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
herring_hatch <- mean(c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365)
Flounder_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Plaice_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Dab_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

for (i in 1:length(samplings)){
  # sampling specific
  quarter <- as.numeric(substr(samplings[i],6,7))
  jday <- (c(1/4,2/4,3/4,4/4)-1/8)[quarter]
  year.i <- as.numeric(substr(samplings[i],1,4))
  # cod
  ####
  cod.s <- c(vbgrCod(jday-cod_hatch),vbgrCod(jday-cod_hatch+1),vbgrCod(jday-cod_hatch+2),vbgrCod(jday-cod_hatch+3),
             vbgrCod(jday-cod_hatch+4),vbgrCod(jday-cod_hatch+5),vbgrCod(jday-cod_hatch+6),vbgrCod(jday-cod_hatch+7))*10
  cod.age <- 0:7
  cod.i <- N_ageCod(Cod,cod.age,year.i)[jday*365,1:length(cod.age)]
  cod.sd<- c(vbgr.sdCod(jday-cod_hatch,100000),vbgr.sdCod(jday-cod_hatch+1,100000),vbgr.sdCod(jday-cod_hatch+2,100000),vbgr.sdCod(jday-cod_hatch+3,100000),
             vbgr.sdCod(jday-cod_hatch+4,100000),vbgr.sdCod(jday-cod_hatch+5,100000),vbgr.sdCod(jday-cod_hatch+6,100000),vbgr.sdCod(jday-cod_hatch+7,100000))*10
  
  if (min(cod.sd)<0){
    cod.sd[cod.sd<0] <-rep(0,length(cod.sd[cod.sd<0]))
  }
  available.cod <- c(floor(rnorm(cod.i[1],cod.s[1],cod.sd[1])/width)*width+width/2,floor(rnorm(cod.i[2],cod.s[2],cod.sd[2])/width)*width+width/2,
                     floor(rnorm(cod.i[3],cod.s[3],cod.sd[3])/width)*width+width/2,floor(rnorm(cod.i[4],cod.s[4],cod.sd[4])/width)*width+width/2,
                     floor(rnorm(cod.i[5],cod.s[5],cod.sd[5])/width)*width+width/2,floor(rnorm(cod.i[6],cod.s[6],cod.sd[6])/width)*width+width/2,
                     floor(rnorm(cod.i[7],cod.s[7],cod.sd[7])/width)*width+width/2,floor(rnorm(cod.i[8],cod.s[8],cod.sd[8])/width)*width+width/2)
  available.cod[available.cod<0] <- rep(width/2,length(available.cod[available.cod<0]))
  
  cod.av <- as.data.frame(table(available.cod))
  cod.av$N <- as.numeric(cod.av$Freq)
  cod.av$l <- as.numeric(as.character(cod.av$available.cod))/10 # cm
  cod.coef <- LW %>% filter(species=="cod" & year==year.i)
  if(length(cod.coef$year)==0){
    cod.coef <- LW %>% filter(species=="cod" & year==min(year))
  } 
  cod.av$B <- cod.av$N*(cod.coef$a*cod.av$l^cod.coef$b) # g
  cod.av$avB <- pmax(0,corm_cod.pref(cod.av$l*10))*cod.av$B
  ####
  # herring
  ####
  herring.s <- c(vbgrHerring(jday-herring_hatch),vbgrHerring(jday-herring_hatch+1),
                 vbgrHerring(jday-herring_hatch+2),vbgrHerring(jday-herring_hatch+3),
                 vbgrHerring(jday-herring_hatch+4),vbgrHerring(jday-herring_hatch+5),
                 vbgrHerring(jday-herring_hatch+6),vbgrHerring(jday-herring_hatch+7),
                 vbgrHerring(jday-herring_hatch+8))*10
  herring.age <- 0:8
  herring.i <- N_ageHerring(herring,herring.age,year.i)[jday*365,1:length(herring.age)]
  herring.sd<- c(vbgr.sdHerring(jday-herring_hatch,100000),vbgr.sdHerring(jday-herring_hatch+1,100000),
                 vbgr.sdHerring(jday-herring_hatch+2,100000),vbgr.sdHerring(jday-herring_hatch+3,100000),
                 vbgr.sdHerring(jday-herring_hatch+4,100000),vbgr.sdHerring(jday-herring_hatch+5,100000),
                 vbgr.sdHerring(jday-herring_hatch+6,100000),vbgr.sdHerring(jday-herring_hatch+7,100000),
                 vbgr.sdHerring(jday-herring_hatch+8,100000))*10
  
  if (min(herring.sd)<0){
    herring.sd[herring.sd<0] <-rep(0,length(herring.sd[herring.sd<0]))
  }
  available.herring <- c(floor(rnorm(herring.i[1],herring.s[1],herring.sd[1])/width)*width+width/2,floor(rnorm(herring.i[2],herring.s[2],herring.sd[2])/width)*width+width/2,
                         floor(rnorm(herring.i[3],herring.s[3],herring.sd[3])/width)*width+width/2,floor(rnorm(herring.i[4],herring.s[4],herring.sd[4])/width)*width+width/2,
                         floor(rnorm(herring.i[5],herring.s[5],herring.sd[5])/width)*width+width/2,floor(rnorm(herring.i[6],herring.s[6],herring.sd[6])/width)*width+width/2,
                         floor(rnorm(herring.i[7],herring.s[7],herring.sd[7])/width)*width+width/2,floor(rnorm(herring.i[8],herring.s[8],herring.sd[8])/width)*width+width/2,
                         floor(rnorm(herring.i[9],herring.s[9],herring.sd[9])/width)*width+width/2)
  available.herring[available.herring<0] <- rep(width/2,length(available.herring[available.herring<0]))
  
  herring.av <- as.data.frame(table(available.herring))
  herring.av$N <- as.numeric(herring.av$Freq)
  herring.av$l <- as.numeric(as.character(herring.av$available.herring))/10 # cm
  herring.coef <- LW %>% filter(species=="herring")
  herring.av$B <- herring.av$N*(herring.coef$a*herring.av$l^herring.coef$b) # g
  herring.av$avB <- pmax(0,corm_herring.pref(herring.av$l*10))*herring.av$B
  ####
  
  # flatfish
  ####
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
  fl.i <- Flatfish %>% filter(year==year.i & yday==round(jday*365))
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
  
  ####
  #flounder
  flounder.av <- as.data.frame(table(availableFlounder))
  flounder.av$N <- as.numeric(flounder.av$Freq)
  flounder.av$l <- as.numeric(as.character(flounder.av$availableFlounder))/10 # cm
  flounder.coef <- LW %>% filter(species=="flounder" & year==year.i)
  if(length(flounder.coef$year)==0){
    flounder.coef <- LW %>% filter(species=="flounder" & year==1992)
  }  
  flounder.av$B <- flounder.av$N*(flounder.coef$a*flounder.av$l^flounder.coef$b) # g
  flounder.av$avB <- pmax(0,corm_flatfish.pref(flounder.av$l*10))*flounder.av$B
  
  #plaice
  plaice.av <- as.data.frame(table(availablePlaice))
  plaice.av$N <- as.numeric(plaice.av$Freq)
  plaice.av$l <- as.numeric(as.character(plaice.av$availablePlaice))/10 # cm
  plaice.coef <- LW %>% filter(species=="plaice" & year==year.i)
  if(length(plaice.coef$year)==0){
    plaice.coef <- LW %>% filter(species=="plaice" & year==1994)
  }  
  plaice.av$B <- plaice.av$N*(plaice.coef$a*plaice.av$l^plaice.coef$b) # g
  plaice.av$avB <- pmax(0,corm_flatfish.pref(plaice.av$l*10))*plaice.av$B
  
  #dab
  dab.av <- as.data.frame(table(availableDab))
  dab.av$N <- as.numeric(dab.av$Freq)
  dab.av$l <- as.numeric(as.character(dab.av$availableDab))/10 # cm
  dab.coef <- LW %>% filter(species=="dab" & year==year.i)
  if(length(dab.coef$year)==0){
    dab.coef <- LW %>% filter(species=="dab" & year==1994)
  }  
  dab.av$B <- dab.av$N*(dab.coef$a*dab.av$l^dab.coef$b) # g
  dab.av$avB <- pmax(0,corm_flatfish.pref(dab.av$l*10))*dab.av$B
  cormFood$sampling[((i-1)*n_prey+1):(i*n_prey)] <- samplings[i]
  cormFood$species[((i-1)*n_prey+1):(i*n_prey)] <- prey
  cormFood$biomass[(i-1)*n_prey+1] <- sum(cod.av$avB)
  cormFood$biomass[(i-1)*n_prey+2] <- sum(herring.av$avB)
  cormFood$biomass[(i-1)*n_prey+3] <- sum(flounder.av$avB)
  cormFood$biomass[(i-1)*n_prey+4] <- sum(plaice.av$avB)
  cormFood$biomass[(i-1)*n_prey+5] <- sum(dab.av$avB)
  
  print(i)
}

write.table(cormFood,"Cormfood_sim.csv",row.names = FALSE,sep=';')
#####






# simulate seal food index each day 1991-2023
#####
prey <- c("cod","herring","flounder","plaice","dab")
samplings <- paste(rep(1991:2023,each=365),rep(1:365,length(1991:2023)))
n_samplings <- length(samplings)
n_prey <- length(prey)


sealFood <- data.frame(sampling=rep("O",n_samplings*n_prey),
                       biomass=rep(0,n_samplings*n_prey),
                       species=rep("O",n_samplings*n_prey))

cod_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
herring_hatch <- mean(c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365)
Flounder_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Plaice_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Dab_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

for (i in 1:length(samplings)){
  # sampling specific
  jday <- as.numeric(substr(samplings[i],6,nchar(samplings[i])))/365
  year.i <- as.numeric(substr(samplings[i],1,4))
  # cod
  ####
  cod.s <- c(vbgrCod(jday-cod_hatch),vbgrCod(jday-cod_hatch+1),vbgrCod(jday-cod_hatch+2),vbgrCod(jday-cod_hatch+3),
             vbgrCod(jday-cod_hatch+4),vbgrCod(jday-cod_hatch+5),vbgrCod(jday-cod_hatch+6),vbgrCod(jday-cod_hatch+7))*10
  cod.age <- 0:7
  cod.i <- N_ageCod(Cod,cod.age,year.i)[jday*365,1:length(cod.age)]
  cod.sd<- c(vbgr.sdCod(jday-cod_hatch,100000),vbgr.sdCod(jday-cod_hatch+1,100000),vbgr.sdCod(jday-cod_hatch+2,100000),vbgr.sdCod(jday-cod_hatch+3,100000),
             vbgr.sdCod(jday-cod_hatch+4,100000),vbgr.sdCod(jday-cod_hatch+5,100000),vbgr.sdCod(jday-cod_hatch+6,100000),vbgr.sdCod(jday-cod_hatch+7,100000))*10
  
  if (min(cod.sd)<0){
    cod.sd[cod.sd<0] <-rep(0,length(cod.sd[cod.sd<0]))
  }
  available.cod <- c(floor(rnorm(cod.i[1],cod.s[1],cod.sd[1])/width)*width+width/2,floor(rnorm(cod.i[2],cod.s[2],cod.sd[2])/width)*width+width/2,
                     floor(rnorm(cod.i[3],cod.s[3],cod.sd[3])/width)*width+width/2,floor(rnorm(cod.i[4],cod.s[4],cod.sd[4])/width)*width+width/2,
                     floor(rnorm(cod.i[5],cod.s[5],cod.sd[5])/width)*width+width/2,floor(rnorm(cod.i[6],cod.s[6],cod.sd[6])/width)*width+width/2,
                     floor(rnorm(cod.i[7],cod.s[7],cod.sd[7])/width)*width+width/2,floor(rnorm(cod.i[8],cod.s[8],cod.sd[8])/width)*width+width/2)
  available.cod[available.cod<0] <- rep(width/2,length(available.cod[available.cod<0]))
  
  cod.av <- as.data.frame(table(available.cod))
  cod.av$N <- as.numeric(cod.av$Freq)
  cod.av$l <- as.numeric(as.character(cod.av$available.cod))/10 # cm
  cod.coef <- LW %>% filter(species=="cod" & year==year.i)
  if(length(cod.coef$year)==0){
    cod.coef <- LW %>% filter(species=="cod" & year==min(year))
  } 
  cod.av$B <- cod.av$N*(cod.coef$a*cod.av$l^cod.coef$b) # g
  cod.av$avB <- pmax(0,seal_cod.pref(cod.av$l*10))*cod.av$B
  ####
  # herring
  ####
  herring.s <- c(vbgrHerring(jday-herring_hatch),vbgrHerring(jday-herring_hatch+1),
                 vbgrHerring(jday-herring_hatch+2),vbgrHerring(jday-herring_hatch+3),
                 vbgrHerring(jday-herring_hatch+4),vbgrHerring(jday-herring_hatch+5),
                 vbgrHerring(jday-herring_hatch+6),vbgrHerring(jday-herring_hatch+7),
                 vbgrHerring(jday-herring_hatch+8))*10
  herring.age <- 0:8
  herring.i <- N_ageHerring(herring,herring.age,year.i)[jday*365,1:length(herring.age)]
  herring.sd<- c(vbgr.sdHerring(jday-herring_hatch,100000),vbgr.sdHerring(jday-herring_hatch+1,100000),
                 vbgr.sdHerring(jday-herring_hatch+2,100000),vbgr.sdHerring(jday-herring_hatch+3,100000),
                 vbgr.sdHerring(jday-herring_hatch+4,100000),vbgr.sdHerring(jday-herring_hatch+5,100000),
                 vbgr.sdHerring(jday-herring_hatch+6,100000),vbgr.sdHerring(jday-herring_hatch+7,100000),
                 vbgr.sdHerring(jday-herring_hatch+8,100000))*10
  
  if (min(herring.sd)<0){
    herring.sd[herring.sd<0] <-rep(0,length(herring.sd[herring.sd<0]))
  }
  available.herring <- c(floor(rnorm(herring.i[1],herring.s[1],herring.sd[1])/width)*width+width/2,floor(rnorm(herring.i[2],herring.s[2],herring.sd[2])/width)*width+width/2,
                         floor(rnorm(herring.i[3],herring.s[3],herring.sd[3])/width)*width+width/2,floor(rnorm(herring.i[4],herring.s[4],herring.sd[4])/width)*width+width/2,
                         floor(rnorm(herring.i[5],herring.s[5],herring.sd[5])/width)*width+width/2,floor(rnorm(herring.i[6],herring.s[6],herring.sd[6])/width)*width+width/2,
                         floor(rnorm(herring.i[7],herring.s[7],herring.sd[7])/width)*width+width/2,floor(rnorm(herring.i[8],herring.s[8],herring.sd[8])/width)*width+width/2,
                         floor(rnorm(herring.i[9],herring.s[9],herring.sd[9])/width)*width+width/2)
  available.herring[available.herring<0] <- rep(width/2,length(available.herring[available.herring<0]))
  
  herring.av <- as.data.frame(table(available.herring))
  herring.av$N <- as.numeric(herring.av$Freq)
  herring.av$l <- as.numeric(as.character(herring.av$available.herring))/10 # cm
  herring.coef <- LW %>% filter(species=="herring")
  herring.av$B <- herring.av$N*(herring.coef$a*herring.av$l^herring.coef$b) # g
  herring.av$avB <- pmax(0,seal_herring.pref(herring.av$l*10))*herring.av$B
  ####
  
  # flatfish
  ####
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
  fl.i <- Flatfish %>% filter(year==year.i & yday==round(jday*365))
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
  
  ####
  #flounder
  flounder.av <- as.data.frame(table(availableFlounder))
  flounder.av$N <- as.numeric(flounder.av$Freq)
  flounder.av$l <- as.numeric(as.character(flounder.av$availableFlounder))/10 # cm
  flounder.coef <- LW %>% filter(species=="flounder" & year==year.i)
  if(length(flounder.coef$year)==0){
    flounder.coef <- LW %>% filter(species=="flounder" & year==1992)
  }  
  flounder.av$B <- flounder.av$N*(flounder.coef$a*flounder.av$l^flounder.coef$b) # g
  flounder.av$avB <- pmax(0,seal_flatfish.pref(flounder.av$l*10))*flounder.av$B
  
  #plaice
  plaice.av <- as.data.frame(table(availablePlaice))
  plaice.av$N <- as.numeric(plaice.av$Freq)
  plaice.av$l <- as.numeric(as.character(plaice.av$availablePlaice))/10 # cm
  plaice.coef <- LW %>% filter(species=="plaice" & year==year.i)
  if(length(plaice.coef$year)==0){
    plaice.coef <- LW %>% filter(species=="plaice" & year==1994)
  }  
  plaice.av$B <- plaice.av$N*(plaice.coef$a*plaice.av$l^plaice.coef$b) # g
  plaice.av$avB <- pmax(0,seal_flatfish.pref(plaice.av$l*10))*plaice.av$B
  
  #dab
  dab.av <- as.data.frame(table(availableDab))
  dab.av$N <- as.numeric(dab.av$Freq)
  dab.av$l <- as.numeric(as.character(dab.av$availableDab))/10 # cm
  dab.coef <- LW %>% filter(species=="dab" & year==year.i)
  if(length(dab.coef$year)==0){
    dab.coef <- LW %>% filter(species=="dab" & year==1994)
  }  
  dab.av$B <- dab.av$N*(dab.coef$a*dab.av$l^dab.coef$b) # g
  dab.av$avB <- pmax(0,seal_flatfish.pref(dab.av$l*10))*dab.av$B
  sealFood$sampling[((i-1)*n_prey+1):(i*n_prey)] <- samplings[i]
  sealFood$species[((i-1)*n_prey+1):(i*n_prey)] <- prey
  sealFood$biomass[(i-1)*n_prey+1] <- sum(cod.av$avB)
  sealFood$biomass[(i-1)*n_prey+2] <- sum(herring.av$avB)
  sealFood$biomass[(i-1)*n_prey+3] <- sum(flounder.av$avB)
  sealFood$biomass[(i-1)*n_prey+4] <- sum(plaice.av$avB)
  sealFood$biomass[(i-1)*n_prey+5] <- sum(dab.av$avB)
  
  print(i)
}

food.sim <- sealFood[-(166:170),]
write.table(food.sim,"food_sim.csv",row.names = FALSE,sep=';')
#####

# simulate cormorant food index each day 1991-2023
#####
prey <- c("cod","herring","flounder","plaice","dab")
samplings <- paste(rep(1991:2023,each=365),rep(1:365,length(1991:2023)))
n_samplings <- length(samplings)
n_prey <- length(prey)


cormFood <- data.frame(sampling=rep("O",n_samplings*n_prey),
                       biomass=rep(0,n_samplings*n_prey),
                       species=rep("O",n_samplings*n_prey))

cod_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
herring_hatch <- mean(c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365)
Flounder_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Plaice_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Dab_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

for (i in 1:length(samplings)){
  # sampling specific
  jday <- as.numeric(substr(samplings[i],6,nchar(samplings[i])))/365
  year.i <- as.numeric(substr(samplings[i],1,4))
  # cod
  ####
  cod.s <- c(vbgrCod(jday-cod_hatch),vbgrCod(jday-cod_hatch+1),vbgrCod(jday-cod_hatch+2),vbgrCod(jday-cod_hatch+3),
             vbgrCod(jday-cod_hatch+4),vbgrCod(jday-cod_hatch+5),vbgrCod(jday-cod_hatch+6),vbgrCod(jday-cod_hatch+7))*10
  cod.age <- 0:7
  cod.i <- N_ageCod(Cod,cod.age,year.i)[jday*365,1:length(cod.age)]
  cod.sd<- c(vbgr.sdCod(jday-cod_hatch,100000),vbgr.sdCod(jday-cod_hatch+1,100000),vbgr.sdCod(jday-cod_hatch+2,100000),vbgr.sdCod(jday-cod_hatch+3,100000),
             vbgr.sdCod(jday-cod_hatch+4,100000),vbgr.sdCod(jday-cod_hatch+5,100000),vbgr.sdCod(jday-cod_hatch+6,100000),vbgr.sdCod(jday-cod_hatch+7,100000))*10
  
  if (min(cod.sd)<0){
    cod.sd[cod.sd<0] <-rep(0,length(cod.sd[cod.sd<0]))
  }
  available.cod <- c(floor(rnorm(cod.i[1],cod.s[1],cod.sd[1])/width)*width+width/2,floor(rnorm(cod.i[2],cod.s[2],cod.sd[2])/width)*width+width/2,
                     floor(rnorm(cod.i[3],cod.s[3],cod.sd[3])/width)*width+width/2,floor(rnorm(cod.i[4],cod.s[4],cod.sd[4])/width)*width+width/2,
                     floor(rnorm(cod.i[5],cod.s[5],cod.sd[5])/width)*width+width/2,floor(rnorm(cod.i[6],cod.s[6],cod.sd[6])/width)*width+width/2,
                     floor(rnorm(cod.i[7],cod.s[7],cod.sd[7])/width)*width+width/2,floor(rnorm(cod.i[8],cod.s[8],cod.sd[8])/width)*width+width/2)
  available.cod[available.cod<0] <- rep(width/2,length(available.cod[available.cod<0]))
  
  cod.av <- as.data.frame(table(available.cod))
  cod.av$N <- as.numeric(cod.av$Freq)
  cod.av$l <- as.numeric(as.character(cod.av$available.cod))/10 # cm
  cod.coef <- LW %>% filter(species=="cod" & year==year.i)
  if(length(cod.coef$year)==0){
    cod.coef <- LW %>% filter(species=="cod" & year==min(year))
  } 
  cod.av$B <- cod.av$N*(cod.coef$a*cod.av$l^cod.coef$b) # g
  cod.av$avB <- pmax(0,corm_cod.pref(cod.av$l*10))*cod.av$B
  ####
  # herring
  ####
  herring.s <- c(vbgrHerring(jday-herring_hatch),vbgrHerring(jday-herring_hatch+1),
                 vbgrHerring(jday-herring_hatch+2),vbgrHerring(jday-herring_hatch+3),
                 vbgrHerring(jday-herring_hatch+4),vbgrHerring(jday-herring_hatch+5),
                 vbgrHerring(jday-herring_hatch+6),vbgrHerring(jday-herring_hatch+7),
                 vbgrHerring(jday-herring_hatch+8))*10
  herring.age <- 0:8
  herring.i <- N_ageHerring(herring,herring.age,year.i)[jday*365,1:length(herring.age)]
  herring.sd<- c(vbgr.sdHerring(jday-herring_hatch,100000),vbgr.sdHerring(jday-herring_hatch+1,100000),
                 vbgr.sdHerring(jday-herring_hatch+2,100000),vbgr.sdHerring(jday-herring_hatch+3,100000),
                 vbgr.sdHerring(jday-herring_hatch+4,100000),vbgr.sdHerring(jday-herring_hatch+5,100000),
                 vbgr.sdHerring(jday-herring_hatch+6,100000),vbgr.sdHerring(jday-herring_hatch+7,100000),
                 vbgr.sdHerring(jday-herring_hatch+8,100000))*10
  
  if (min(herring.sd)<0){
    herring.sd[herring.sd<0] <-rep(0,length(herring.sd[herring.sd<0]))
  }
  available.herring <- c(floor(rnorm(herring.i[1],herring.s[1],herring.sd[1])/width)*width+width/2,floor(rnorm(herring.i[2],herring.s[2],herring.sd[2])/width)*width+width/2,
                         floor(rnorm(herring.i[3],herring.s[3],herring.sd[3])/width)*width+width/2,floor(rnorm(herring.i[4],herring.s[4],herring.sd[4])/width)*width+width/2,
                         floor(rnorm(herring.i[5],herring.s[5],herring.sd[5])/width)*width+width/2,floor(rnorm(herring.i[6],herring.s[6],herring.sd[6])/width)*width+width/2,
                         floor(rnorm(herring.i[7],herring.s[7],herring.sd[7])/width)*width+width/2,floor(rnorm(herring.i[8],herring.s[8],herring.sd[8])/width)*width+width/2,
                         floor(rnorm(herring.i[9],herring.s[9],herring.sd[9])/width)*width+width/2)
  available.herring[available.herring<0] <- rep(width/2,length(available.herring[available.herring<0]))
  
  herring.av <- as.data.frame(table(available.herring))
  herring.av$N <- as.numeric(herring.av$Freq)
  herring.av$l <- as.numeric(as.character(herring.av$available.herring))/10 # cm
  herring.coef <- LW %>% filter(species=="herring")
  herring.av$B <- herring.av$N*(herring.coef$a*herring.av$l^herring.coef$b) # g
  herring.av$avB <- pmax(0,corm_herring.pref(herring.av$l*10))*herring.av$B
  ####
  
  # flatfish
  ####
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
  fl.i <- Flatfish %>% filter(year==year.i & yday==round(jday*365))
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
  
  ####
  #flounder
  flounder.av <- as.data.frame(table(availableFlounder))
  flounder.av$N <- as.numeric(flounder.av$Freq)
  flounder.av$l <- as.numeric(as.character(flounder.av$availableFlounder))/10 # cm
  flounder.coef <- LW %>% filter(species=="flounder" & year==year.i)
  if(length(flounder.coef$year)==0){
    flounder.coef <- LW %>% filter(species=="flounder" & year==1992)
  }  
  flounder.av$B <- flounder.av$N*(flounder.coef$a*flounder.av$l^flounder.coef$b) # g
  flounder.av$avB <- pmax(0,corm_flatfish.pref(flounder.av$l*10))*flounder.av$B
  
  #plaice
  plaice.av <- as.data.frame(table(availablePlaice))
  plaice.av$N <- as.numeric(plaice.av$Freq)
  plaice.av$l <- as.numeric(as.character(plaice.av$availablePlaice))/10 # cm
  plaice.coef <- LW %>% filter(species=="plaice" & year==year.i)
  if(length(plaice.coef$year)==0){
    plaice.coef <- LW %>% filter(species=="plaice" & year==1994)
  }  
  plaice.av$B <- plaice.av$N*(plaice.coef$a*plaice.av$l^plaice.coef$b) # g
  plaice.av$avB <- pmax(0,corm_flatfish.pref(plaice.av$l*10))*plaice.av$B
  
  #dab
  dab.av <- as.data.frame(table(availableDab))
  dab.av$N <- as.numeric(dab.av$Freq)
  dab.av$l <- as.numeric(as.character(dab.av$availableDab))/10 # cm
  dab.coef <- LW %>% filter(species=="dab" & year==year.i)
  if(length(dab.coef$year)==0){
    dab.coef <- LW %>% filter(species=="dab" & year==1994)
  }  
  dab.av$B <- dab.av$N*(dab.coef$a*dab.av$l^dab.coef$b) # g
  dab.av$avB <- pmax(0,corm_flatfish.pref(dab.av$l*10))*dab.av$B
  cormFood$sampling[((i-1)*n_prey+1):(i*n_prey)] <- samplings[i]
  cormFood$species[((i-1)*n_prey+1):(i*n_prey)] <- prey
  cormFood$biomass[(i-1)*n_prey+1] <- sum(cod.av$avB)
  cormFood$biomass[(i-1)*n_prey+2] <- sum(herring.av$avB)
  cormFood$biomass[(i-1)*n_prey+3] <- sum(flounder.av$avB)
  cormFood$biomass[(i-1)*n_prey+4] <- sum(plaice.av$avB)
  cormFood$biomass[(i-1)*n_prey+5] <- sum(dab.av$avB)
  
  print(i)
}

food.sim <- cormFood[-(166:170),]
write.table(food.sim,"food_sim.csv",row.names = FALSE,sep=';')
#####









