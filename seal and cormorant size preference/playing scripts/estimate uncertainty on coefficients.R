source('sd estimation.R')
rm(list=setdiff(ls(),c('df')))
data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="") # data working directory

###########
# cod 
##########
# Note: I use Model (b) from McQueen et al., 2019, because this has a more realistic growth rate at you ages
# define growth parameters
Linf <- 97.9 # cm. McQueen et al., 2019
L_hatch <- 0.4 # length at hatching, cm (Pepin et al, 1997)
k <- 0.22 # vbgr growth rate, from McQueen
u <- 0.23 # Amplitude of seasonality
w <- 0.91 # Timing of peak growth
#t_hatch <- ((46+135)/2)/365 # time of hatching, Julian day / 365 - spawning from Feb. to May (Fiskeatlas), average set to peak spawning
t_hatch <- ((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to May (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
# simulate larvae distribution by assuming 95 % CI between Jan. 15th and July 15th 
r <- round(rnorm(1000000,mean=t_hatch*365,sd=(197-t_hatch*365)/1.96))
idx <- r<=197 # index below the 97.5% quantile
t <- sort(unique(r[idx]))/365 # larvae hatching days
age <- 197/365-t # larvae age on May 15th
phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
phi_hatch <- phi_t # seasonal variability in growth at hatching
L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf # larvae length on July 15th
df.hatch <- data.frame(age = age,length=L_t,n=as.matrix(table(r[idx]/365))[,1]) # combine in dataframe
mean.weighted <- sum(df.hatch$length*df.hatch$n)/sum(df.hatch$n) # calculate mean larvae length on July 15th
sd.weighted <- sqrt(sum(df.hatch$n*(df.hatch$length-mean.weighted)^2)/sum(df.hatch$n)) # calculate sd for larvae length on July 15th
#sd(rep(df.hatch$length,times=df.hatch$n)) # test sd is calculated correctly
#sd.weighted

# Now estimate growth parameter uncertainties
cod.sd <- df$sd[df$species=="cod"] # cod sd at age
cod.sd <- c(sd.weighted,cod.sd) # add the calculated sd after spawning
cod.age <- df$yd[df$species=="cod"]-t_hatch # get ages
cod.age <- c(197/365-t_hatch,cod.age)  # add age of the calculated sd after spawning
plot(cod.age,cod.sd,ylim=c(0,20),xlab="Age [years]",ylab="SD [cm]]",main="Cod")
points(cod.age[1],cod.sd[1],pch=19,col="red")
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
# Function to calculate sd of L as a function of age. Based on partial derivatives of L (delta method)
calculate_sd <- function(age, k_sd,age_sd,Linf_sd) {
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
   # Partial derivatives
  dL_dLinf <- -exp(-k * (phi_t - phi_hatch +age)) + 1
  dL_dk <- (L_hatch - Linf) * (phi_hatch - phi_t -age ) * exp(-k * (phi_t +age - phi_hatch))
  dL_dage <- -(L_hatch - Linf)*k*exp(-k*(phi_t + age - phi_hatch))
  
  # Variance
  var_L <- (dL_dLinf * Linf_sd)^2 + (dL_dk * k_sd)^2 + (dL_dage * age_sd)^2
  
  # Standard deviation
  sd_L <- sqrt(var_L)
  
  return(sd_L)
} #Original

par <- list(logK_sd = -4,logA_sd=-2,logLinf_sd=3)

obj_fn <- function(par) {
  k_sd <- exp(par[["logK_sd"]])
  age_sd <- exp(par[["logA_sd"]])
  Linf_sd <- exp(par[["logLinf_sd"]])
  pred <- calculate_sd(cod.age, k_sd, age_sd,Linf_sd)
  sum((cod.sd - pred)^2)
}
opt <- nlminb(par, objective = obj_fn)
cod_sd_coef <- exp(opt$par)
exp(opt$par)
par(mfrow=c(1,1))
t.c <- (0:1500)/100
plot(t.c, calculate_sd(t.c,exp(opt$par[1]),exp(opt$par[2]),exp(opt$par[3])),
     type = "l", col = "black", lwd = 2,
     ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,8.5),ylim=c(0,20),main="cod")
points(cod.age,cod.sd)

# simulation test
sd.0 <- matrix(rep(0,length(t.c)*100000),nrow = 100000)
t <- (t.c+t_hatch)
X_t <- t-floor(t) # time of year [0,1]
phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth

for (i in 1:100000){
  k.i <- rnorm(1,mean=k,sd=exp(opt$par[1]))
  Linf.i <- rnorm(1,mean=Linf,sd=exp(opt$par[3]))
  fish.age <- t.c+rnorm(1,0,exp(opt$par[2]))
  L <- (L_hatch-Linf.i)*exp(-k.i*(phi_t+fish.age-phi_hatch))+Linf.i
  #lines(t,L,col='red')
  sd.0[i,] <- L
  
}
lines(t.c,apply(sd.0,2,sd),lwd=2,col="red",lty = 'dashed')

cod_sd_coef
###########
rm(list=setdiff(ls(),c('df','cod_sd_coef')))

###########
# herring 
##########
# define growth parameters
#Parameters 2005-2010 for Western Baltic spring spawning herring (Gröhsler et al., 2013)
Linf <- 30.57 # cm
k <- 0.453 # vbgr growth rate
hatching <- c(92,92,94,96,99,103,106,109,111,118,118,122,126,134,135,143)/365
t_hatch <- mean(hatching) # From figure 4 in Polte et al., 2021
L_hatch <- 0.645 # figure 2 in Bauer et al., 2014
u <- 0.23 # Amplitude of seasonality (same as cod)
w <- 0.91 # Timing of peak growth
hatch_dur <- mean(c(33,37,43,44,46,48,49,52,54,56,57,63,65,66,68,69))/365 # From figure 4B in Polte et al., 2021
r <- round(rnorm(1000000,mean=t_hatch*365,sd=365*((hatch_dur/2)/1.96)))/365
#hist(r)
#quantile(r,c(0.025,0.975))
idx <- r<=t_hatch+hatch_dur/2 # index below the 97.5% quantile
t <- sort(unique(r[idx])) # larvae hatching days

age <- t_hatch+hatch_dur/2-t # larvae age on May 15th
phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf # larvae length on May 15th
df.hatch <- data.frame(age = age,length=L_t,n=as.matrix(table(r[idx]/365))[,1]) # combine in dataframe
mean.weighted <- sum(df.hatch$length*df.hatch$n)/sum(df.hatch$n) # calculate mean larvae length on May 15th
sd.weighted <- sqrt(sum(df.hatch$n*(df.hatch$length-mean.weighted)^2)/sum(df.hatch$n)) # calculate sd for larvae length on May 15th
#sd(rep(df.hatch$length,times=df.hatch$n)) # test sd is calculated correctly
#sd.weighted

# Now estimate growth parameter uncertainties
herring.sd <- df$sd[df$species=="herring"] # herring sd at age
herring.sd <- c(sd.weighted,herring.sd) # add the calculated sd after spawning
herring.age <- df$yd[df$species=="herring"]-t_hatch # get ages
herring.age <- c(hatch_dur/2,herring.age)  # add age of the calculated sd after spawning
plot(herring.age,herring.sd,ylim=c(0,4),xlab="Age [years]",ylab="SD [cm]]",main="Herring")
points(herring.age[1],herring.sd[1],pch=19,col="red")
herring.age <- herring.age[-(16:19)]
herring.sd <- herring.sd[-(16:19)]
# Function to calculate sd of L as a function of age. Based on partial derivatives of L (delta method)
calculate_sd <- function(age, k_sd,age_sd,Linf_sd) {
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  # Partial derivatives
  dL_dLinf <- -exp(-k * (phi_t - phi_hatch +age)) + 1
  dL_dk <- (L_hatch - Linf) * (phi_hatch - phi_t -age ) * exp(-k * (phi_t +age - phi_hatch))
  dL_dage <- -(L_hatch - Linf)*k*exp(-k*(phi_t + age - phi_hatch))
  
  # Variance
  var_L <- (dL_dLinf * Linf_sd)^2 + (dL_dk * k_sd)^2 + (dL_dage * age_sd)^2
  
  # Standard deviation
  sd_L <- sqrt(var_L)
  
  return(sd_L)
} #Original

par <- list(logK_sd = 0,logA_sd=0,logLinf_sd=0)

obj_fn <- function(par) {
  k_sd <- exp(par[["logK_sd"]])
  age_sd <- exp(par[["logA_sd"]])
  Linf_sd <- exp(par[["logLinf_sd"]])
  pred <- calculate_sd(herring.age, k_sd, age_sd,Linf_sd)
  sum((herring.sd - pred)^2)
}
opt <- nlminb(par, objective = obj_fn)
herring_sd_coef <- exp(opt$par)
exp(opt$par)
par(mfrow=c(1,1))
t.c <- (0:1500)/100
plot(t.c, calculate_sd(t.c,exp(opt$par[1]),exp(opt$par[2]),exp(opt$par[3])),
     type = "l", col = "black", lwd = 2,
     ylab = 'SD [cm]',xlab='Age [years]',xlim=c(0,12),ylim=c(0,5),main="Herring")
points(herring.age,herring.sd)
# simulation test
sd.0 <- matrix(rep(0,length(t.c)*100000),nrow = 100000)
phi_t <- u*sin(2*pi*((t.c+t_hatch)-w))/(2*pi) # seasonal variability in growth

for (i in 1:100000){
  k.i <- rnorm(1,mean=k,sd=exp(opt$par[1]))
  Linf.i <- rnorm(1,mean=Linf,sd=exp(opt$par[3]))
  fish.age <- t.c+rnorm(1,0,exp(opt$par[2]))
  L <- (L_hatch-Linf.i)*exp(-k.i*(phi_t+fish.age-phi_hatch))+Linf.i
  #lines(t,L,col='red')
  sd.0[i,] <- L
  
}
lines(t.c,apply(sd.0,2,sd),lwd=2,col="red",lty = 'dashed')
herring_sd_coef
###########
rm(list=setdiff(ls(),c('df','cod_sd_coef','herring_sd_coef')))

###########
# flounder
##########
# define growth parameters
L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
u <- 0.23 # Amplitude of seasonality (same as cod)
w <- 0.91 # Timing of peak growth
Linf <- 34.0031653
k <- 0.4537312

t_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

# simulate larvae distribution by assuming 95 % CI between Mar. 15th and Jul. 15th 
r <- round(rnorm(1000000,mean=t_hatch*365,sd=(197-t_hatch*365)/1.96))/365
#hist(r)
#quantile(r,c(0.025,0.975))*365 # test
idx <- r<= 197/365# index below the 97.5% quantile
t <- sort(unique(r[idx])) # larvae hatching days

age <- 197/365-t # larvae age on May 15th
phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching


L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf # larvae length on July 15th
df.hatch <- data.frame(age = age,length=L_t,n=as.matrix(table(r[idx]/365))[,1]) # combine in dataframe
mean.weighted <- sum(df.hatch$length*df.hatch$n)/sum(df.hatch$n) # calculate mean larvae length on July 15th
sd.weighted <- sqrt(sum(df.hatch$n*(df.hatch$length-mean.weighted)^2)/sum(df.hatch$n)) # calculate sd for larvae length on July 15th

# Now estimate growth parameter uncertainties
flounder.sd <- df$sd[df$species=="flounder"] # flounder sd at age
flounder.sd <- c(sd.weighted,flounder.sd) # add the calculated sd after spawning
flounder.age <- df$yd[df$species=="flounder"]-t_hatch # get ages
flounder.age <- c(197/365-t_hatch,flounder.age)  # add age of the calculated sd after spawning
plot(flounder.age,flounder.sd,ylim=c(0,10),xlab="Age [years]",ylab="SD [cm]]",main="flounder")
points(flounder.age[1],flounder.sd[1],pch=19,col="red")

# Function to calculate sd of L as a function of age. Based on partial derivatives of L (delta method)
calculate_sd <- function(age, k_sd,age_sd,Linf_sd) {
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  # Partial derivatives
  dL_dLinf <- -exp(-k * (phi_t - phi_hatch +age)) + 1
  dL_dk <- (L_hatch - Linf) * (phi_hatch - phi_t -age ) * exp(-k * (phi_t +age - phi_hatch))
  dL_dage <- -(L_hatch - Linf)*k*exp(-k*(phi_t + age - phi_hatch))
  
  # Variance
  var_L <- (dL_dLinf * Linf_sd)^2 + (dL_dk * k_sd)^2 + (dL_dage * age_sd)^2
  
  # Standard deviation
  sd_L <- sqrt(var_L)
  
  return(sd_L)
} #Original

par <- list(logK_sd = 0,logA_sd=0,logLinf_sd=0)

obj_fn <- function(par) {
  k_sd <- exp(par[["logK_sd"]])
  age_sd <- exp(par[["logA_sd"]])
  Linf_sd <- exp(par[["logLinf_sd"]])
  pred <- calculate_sd(flounder.age, k_sd, age_sd,Linf_sd)
  sum((flounder.sd - pred)^2)
}
opt <- nlminb(par, objective = obj_fn)
flounder_sd_coef <- exp(opt$par)

exp(opt$par)
par(mfrow=c(1,1))
t.c <- (0:1500)/100
plot(t.c, calculate_sd(t.c,exp(opt$par[1]),exp(opt$par[2]),exp(opt$par[3])),
     type = "l", col = "black", lwd = 2,
     ylab = 'SD [cm]',xlab='Age [years]',xlim=c(0,12),ylim=c(0,8),main="flounder")
points(flounder.age,flounder.sd)
# simulation test
sd.0 <- matrix(rep(0,length(t.c)*100000),nrow = 100000)
phi_t <- u*sin(2*pi*((t.c+t_hatch)-w))/(2*pi) # seasonal variability in growth

for (i in 1:100000){
  k.i <- rnorm(1,mean=k,sd=exp(opt$par[1]))
  Linf.i <- rnorm(1,mean=Linf,sd=exp(opt$par[3]))
  fish.age <- t.c+rnorm(1,0,exp(opt$par[2]))
  L <- (L_hatch-Linf.i)*exp(-k.i*(phi_t+fish.age-phi_hatch))+Linf.i
  #lines(t,L,col='red')
  sd.0[i,] <- L
  
}
lines(t.c,apply(sd.0,2,sd),lwd=2,col="red",lty = 'dashed')
flounder_sd_coef
###########
rm(list=setdiff(ls(),c('df','cod_sd_coef','herring_sd_coef','flounder_sd_coef')))

###########
# plaice
##########
# define growth parameters
L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
u <- 0.23 # Amplitude of seasonality (same as cod)
w <- 0.91 # Timing of peak growth
Linf <- 41.637260
k <- 0.240703

t_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

# simulate larvae distribution by assuming 95 % CI between Nov. 15th and Mar. 15th 
r <- round(rnorm(1000000,mean=t_hatch*365,sd=(75-t_hatch*365)/1.96))/365
#hist(r)
#quantile(r,c(0.025,0.975))*365 # test
idx <- r<= 75/365# index below the 97.5% quantile
t <- sort(unique(r[idx])) # larvae hatching days

age <- 75/365-t # larvae age on Mar 15th
phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching


L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf # larvae length on May 15th
df.hatch <- data.frame(age = age,length=L_t,n=as.matrix(table(r[idx]/365))[,1]) # combine in dataframe
mean.weighted <- sum(df.hatch$length*df.hatch$n)/sum(df.hatch$n) # calculate mean larvae length on May 15th
sd.weighted <- sqrt(sum(df.hatch$n*(df.hatch$length-mean.weighted)^2)/sum(df.hatch$n)) # calculate sd for larvae length on May 15th

# Now estimate growth parameter uncertainties
plaice.sd <- df$sd[df$species=="plaice"] # plaice sd at age
plaice.sd <- c(sd.weighted,plaice.sd) # add the calculated sd after spawning
plaice.age <- df$yd[df$species=="plaice"]-t_hatch # get ages
plaice.age <- c(75/365-t_hatch,plaice.age)  # add age of the calculated sd after spawning
plot(plaice.age,plaice.sd,ylim=c(0,10),xlab="Age [years]",ylab="SD [cm]]",main="plaice")
points(plaice.age[1],plaice.sd[1],pch=19,col="red")

plaice.sd <- plaice.sd[-(23:27)]
plaice.age <- plaice.age[-(23:27)]
# Function to calculate sd of L as a function of age. Based on partial derivatives of L (delta method)
calculate_sd <- function(age, k_sd,age_sd,Linf_sd) {
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  # Partial derivatives
  dL_dLinf <- -exp(-k * (phi_t - phi_hatch +age)) + 1
  dL_dk <- (L_hatch - Linf) * (phi_hatch - phi_t -age ) * exp(-k * (phi_t +age - phi_hatch))
  dL_dage <- -(L_hatch - Linf)*k*exp(-k*(phi_t + age - phi_hatch))
  
  # Variance
  var_L <- (dL_dLinf * Linf_sd)^2 + (dL_dk * k_sd)^2 + (dL_dage * age_sd)^2
  
  # Standard deviation
  sd_L <- sqrt(var_L)
  
  return(sd_L)
} #Original

par <- list(logK_sd = 0,logA_sd=0,logLinf_sd=0)

obj_fn <- function(par) {
  k_sd <- exp(par[["logK_sd"]])
  age_sd <- exp(par[["logA_sd"]])
  Linf_sd <- exp(par[["logLinf_sd"]])
  pred <- calculate_sd(plaice.age, k_sd, age_sd,Linf_sd)
  sum((plaice.sd - pred)^2)
}
opt <- nlminb(par, objective = obj_fn)
plaice_sd_coef <- exp(opt$par)
exp(opt$par)
par(mfrow=c(1,1))
t.c <- (0:1500)/100
plot(t.c, calculate_sd(t.c,exp(opt$par[1]),exp(opt$par[2]),exp(opt$par[3])),
     type = "l", col = "black", lwd = 2,
     ylab = 'SD [cm]',xlab='Age [years]',xlim=c(0,15),ylim=c(0,10),main="plaice")
points(plaice.age,plaice.sd)
# simulation test
sd.0 <- matrix(rep(0,length(t.c)*100000),nrow = 100000)
phi_t <- u*sin(2*pi*((t.c+t_hatch)-w))/(2*pi) # seasonal variability in growth

for (i in 1:100000){
  k.i <- rnorm(1,mean=k,sd=exp(opt$par[1]))
  Linf.i <- rnorm(1,mean=Linf,sd=exp(opt$par[3]))
  fish.age <- t.c+rnorm(1,0,exp(opt$par[2]))
  L <- (L_hatch-Linf.i)*exp(-k.i*(phi_t+fish.age-phi_hatch))+Linf.i
  #lines(t,L,col='red')
  sd.0[i,] <- L
  
}
lines(t.c,apply(sd.0,2,sd),lwd=2,col="red",lty = 'dashed')
plaice_sd_coef
###########
rm(list=setdiff(ls(),c('df','cod_sd_coef','herring_sd_coef','flounder_sd_coef','plaice_sd_coef')))

###########
# dab
##########
# define growth parameters
L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
u <- 0.23 # Amplitude of seasonality (same as cod)
w <- 0.91 # Timing of peak growth
Linf <- 26.5093927
k <- 0.4846326

t_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

# simulate larvae distribution by assuming 95 % CI between Apr. 15th and Aug. 15th 
r <- round(rnorm(1000000,mean=t_hatch*365,sd=(228-t_hatch*365)/1.96))/365
#hist(r)
#quantile(r,c(0.025,0.975))*365 # test
idx <- r<= 228/365# index below the 97.5% quantile
t <- sort(unique(r[idx])) # larvae hatching days

age <- 228/365-t # larvae age on Mar 15th
phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching


L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf # larvae length on Aug. 15th
df.hatch <- data.frame(age = age,length=L_t,n=as.matrix(table(r[idx]/365))[,1]) # combine in dataframe
mean.weighted <- sum(df.hatch$length*df.hatch$n)/sum(df.hatch$n) # calculate mean larvae length on Aug. 15th
sd.weighted <- sqrt(sum(df.hatch$n*(df.hatch$length-mean.weighted)^2)/sum(df.hatch$n)) # calculate sd for larvae length on May 15th

# Now estimate growth parameter uncertainties
dab.sd <- df$sd[df$species=="dab"] # dab sd at age
dab.sd <- c(sd.weighted,dab.sd) # add the calculated sd after spawning
dab.age <- df$yd[df$species=="dab"]-t_hatch # get ages
dab.age <- c(228/365-t_hatch,dab.age)  # add age of the calculated sd after spawning
plot(dab.age,dab.sd,ylim=c(0,10),xlab="Age [years]",ylab="SD [cm]]",main="dab")
points(dab.age[1],dab.sd[1],pch=19,col="red")

# Function to calculate sd of L as a function of age. Based on partial derivatives of L (delta method)
calculate_sd <- function(age, k_sd,age_sd,Linf_sd) {
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  # Partial derivatives
  dL_dLinf <- -exp(-k * (phi_t - phi_hatch +age)) + 1
  dL_dk <- (L_hatch - Linf) * (phi_hatch - phi_t -age ) * exp(-k * (phi_t +age - phi_hatch))
  dL_dage <- -(L_hatch - Linf)*k*exp(-k*(phi_t + age - phi_hatch))
  
  # Variance
  var_L <- (dL_dLinf * Linf_sd)^2 + (dL_dk * k_sd)^2 + (dL_dage * age_sd)^2
  
  # Standard deviation
  sd_L <- sqrt(var_L)
  
  return(sd_L)
} #Original

par <- list(logK_sd = 0,logA_sd=0,logLinf_sd=0)

obj_fn <- function(par) {
  k_sd <- exp(par[["logK_sd"]])
  age_sd <- exp(par[["logA_sd"]])
  Linf_sd <- exp(par[["logLinf_sd"]])
  pred <- calculate_sd(dab.age, k_sd, age_sd,Linf_sd)
  sum((dab.sd - pred)^2)
}
opt <- nlminb(par, objective = obj_fn)
dab_sd_coef <- exp(opt$par)
exp(opt$par)
par(mfrow=c(1,1))
t.c <- (0:1500)/100
plot(t.c, calculate_sd(t.c,exp(opt$par[1]),exp(opt$par[2]),exp(opt$par[3])),
     type = "l", col = "black", lwd = 2,
     ylab = 'SD [cm]',xlab='Age [years]',xlim=c(0,15),ylim=c(0,10),main="dab")
points(dab.age,dab.sd)
# simulation test
sd.0 <- matrix(rep(0,length(t.c)*100000),nrow = 100000)
phi_t <- u*sin(2*pi*((t.c+t_hatch)-w))/(2*pi) # seasonal variability in growth

for (i in 1:100000){
  k.i <- rnorm(1,mean=k,sd=exp(opt$par[1]))
  Linf.i <- rnorm(1,mean=Linf,sd=exp(opt$par[3]))
  fish.age <- t.c+rnorm(1,0,exp(opt$par[2]))
  L <- (L_hatch-Linf.i)*exp(-k.i*(phi_t+fish.age-phi_hatch))+Linf.i
  #lines(t,L,col='red')
  sd.0[i,] <- L
  
}
lines(t.c,apply(sd.0,2,sd),lwd=2,col="red",lty = 'dashed')
dab_sd_coef
###########
rm(list=setdiff(ls(),c('df','cod_sd_coef','herring_sd_coef','flounder_sd_coef','plaice_sd_coef','dab_sd_coef')))

data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="") # data working directory
coefs <- readRDS(paste(data_wd,"coefficients.R",sep=""))
coefs$cod_sd_coef <- cod_sd_coef
coefs$herring_sd_coef <- herring_sd_coef
coefs$flounder_sd_coef <- flounder_sd_coef
coefs$plaice_sd_coef <- plaice_sd_coef
coefs$dab_sd_coef <- dab_sd_coef

saveRDS(coefs,paste(data_wd,"coefficients.R",sep=""))





