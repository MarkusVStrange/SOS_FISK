library(reshape2)
rm(list=setdiff(ls(),c('hl_N')))

cod <- hl_N %>% filter(species=='cod' & Age>1)
cod <- cod %>% filter(!(LngtClass<20 & Age>5)) # remove a 10 cm 8 year old cod. Must be a mistake
cod <- cod %>% filter(Quarter %in% c(1,4) & Age<10)

cod$sampID <- paste(cod$Year,cod$Quarter)



cod$yd <- cod$Age+cod$jday
cod$cohort <- cod$Year-cod$Age

str(cod)
ages <- sort(unique(cod$yd))
cohorts <- sort(unique(cod$cohort))

#ages <- rep(2:10,each=3)+rep(sort(c(unique(cod$jday))),9)




mat <- matrix(NA,ncol=length(cohorts),nrow = length(ages))

n_mat <- mat
mat_mean <- mat

for (i in 1:length(cohorts)){
  
  for (j in 1:length(ages)){
    dada <- cod %>% filter(yd==ages[j] & cohort==cohorts[i])
    L <- rep(dada$LngtClass,times=round(dada$N_age))
    mat[j,i] <- sd(L)
    mat_mean[j,i] <- mean(L)
    n_mat[j,i] <- sum(dada$N_age)
    
  }
}
mat[mat==0] <- NA
mat <- mat[-(14:16),-c(1:3,41)] #remove the last 3 ages and first hohort with no data
n_mat <-  n_mat[-(14:16),-c(1:3,41)]
mat_mean <- mat_mean[-(14:16),-c(1:3,41)]

dat <- list(obs=mat)

par <- list(logSdObs =0,logMu=rep(0,nrow(dat$obs)), transPhi=c(0.1,0.1), logSigma=0, omega=matrix(0, nrow=nrow(dat$obs), ncol=ncol(dat$obs)))

f <- function(par){
  getAll(par,dat)
  
  # all the good stuff goes here
  phi <- 2*plogis(transPhi)-1
  sigma <- exp(logSigma)
  mu <- exp(matrix(rep(logMu,ncol(omega)),ncol = ncol(par$omega)))
  sdObs <- exp(logSdObs)
  idx <- !(is.na(obs))
  d1 <- function(x)dautoreg(x,phi = phi[1],log=TRUE)
  d2 <- function(x)dautoreg(x,phi = phi[2],log=TRUE)
  nll <- -dseparable(d1,d2)(omega,scale=sigma)
  nll <- nll-sum(dnorm((obs[idx]),mean=(mu+omega)[idx],sd=sdObs,log=TRUE))
  
  est <- mu+omega
  ADREPORT(est)
  nll  
}

obj <- MakeADFun(f, par, random = "omega", silent=TRUE)
fit <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
sdr


pl <- as.list(sdr, "Est",report=TRUE)$est
plsd <- as.list(sdr, "Std",report=TRUE)$est


rownames(pl) <- paste0("age", 1:nrow(pl))
colnames(pl) <- paste0("cohort", 1:ncol(pl))
df <- melt(pl, varnames = c("age", "cohort"), value.name = "sd")
p.sd <- ggplot(df, aes(x = cohort, y = age, fill = sd)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C",limits=c(0,25), na.value = "grey90") +  # better gradient
  labs(
    title = "cod sd",
    x = "Cohort",
    y = "age",
    fill = "sd"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

n_mat[n_mat==0] <- NA
n_mat <- log10(n_mat)
rownames(n_mat) <- paste0("age", 1:nrow(n_mat))
colnames(n_mat) <- paste0("cohort", 1:ncol(n_mat))
df_n <- melt(n_mat, varnames = c("age", "cohort"), value.name = "log_CPUE")
p.n <- ggplot(df_n, aes(x = cohort, y = age, fill = log_CPUE)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C",limits=c(0,5), na.value = "grey90") +  # better gradient
  labs(
    title = "cod CPUE",
    x = "Cohort",
    y = "age",
    fill = "log_CPUE"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


rownames(mat) <- paste0("age", 1:nrow(mat))
colnames(mat) <- paste0("cohort", 1:ncol(mat))
df_obs <- melt(mat, varnames = c("age", "cohort"), value.name = "sd")
p.obs <- ggplot(df_obs, aes(x = cohort, y = age, fill = sd)) +
  geom_tile(color = "white") + 
  scale_fill_viridis_c(option = "C",limits=c(0,25), na.value = "grey90") +  # better gradient
  labs(
    title = "cod sd obs",
    x = "Cohort",
    y = "age",
    fill = "sd"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggarrange(p.sd,p.obs,p.n,nrow=1)


par(mfrow=c(1,3))
sd_means <- rowMeans(pl)
plot(ages[-(14:16)],sd_means,xlim=c(0,9),ylim=c(0,12),xlab="age",ylab="mean sd")

mat.idx <- is.na(mat)
pl.na <- pl
pl.na[mat.idx] <- NA
sd_means2 <- rowMeans(pl.na, na.rm = TRUE)
plot(ages[-(14:16)],sd_means2,xlim=c(0,9),ylim=c(0,12),xlab="age",ylab="mean sd")


#cpue_means <- rowMeans(n_mat, na.rm = TRUE)
#plot(ages[-(14:16)],(cpue_means),xlab="age",ylab="mean sd")


sd_obs <- rowMeans(mat, na.rm = TRUE)
plot(ages[-(14:16)],sd_obs,xlim=c(0,9),ylim=c(0,12),xlab="age",ylab="mean sd")

mean_obs <- rowMeans(mat_mean, na.rm = TRUE)
plot(ages[-(14:16)],mean_obs,xlab="age",ylab="mean lenth")

ag<- aggregate(N_age~LngtClass+yd,data=cod,FUN=sum)


x <- rep(0,length(unique(ag$yd)))
for (i in 1:16){
  d <- ag %>% filter(yd %in% sort(unique(ag$yd))[[i]])
  x[i] <- sum((d$LngtClass*d$N_age))/sum(d$N_age)
}
plot(ages[-(14:16)],x[-(14:16)],xlab="age",ylab="mean length")
