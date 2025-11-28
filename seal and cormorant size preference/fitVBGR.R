########################
# a function to fit a VBGR to cohort mean size-at-age and get the individual SD-at-age within the cohort
# 
#' @param dat = a list containing: cohort, t (Age+time of year), N, Length, Species
#' @param weights = a vector with two elements: 1st element is weight on mean nll and 2nd is on cohort sd nll
#'                  1 is no weight; 2 is sqrt(N)-weight
#' 
#' @return 
#' a dataframe containing consumption in ind. pr. age pr. year of prey for the different predators
#########################

fitVBGR <- function(dat, weights=c(1,1)){
  species <- unique(dat$Species)
  dat$t_hatches <- c(cod=61/365,flounder=((75+197)/2)/365,plaice=0.03972603,dab=0.4561644)
  dat$L_hatches <-c(cod=0.4,flounder=mean(c(0.697,0.662)),plaice=mean(c(0.697,0.662)),dab=mean(c(0.697,0.662)))
  dat$t1s <- c(cod=1,flounder=2,plaice=2,dab=2)
  dat$t3s <- c(cod=5,flounder=10,plaice=10,dab=10)
  
  n_cohorts <- length(unique(dat$cohort))
  par <- list()
  par$logL1.1 <- log(15)
  par$logL2.1 <- log(35)
  par$logL3.1 <- log(50)
  par$logL1.rand <- rep(log(15),n_cohorts-1)
  par$logL2.rand <- rep(log(35),n_cohorts-1)
  par$logL3.rand <- rep(log(50),n_cohorts-1)
  
  par$logSD_m <- 0 # deviation from mean fit
  par$logSD_c <- 0 # year class SD
  par$logSDage <- 0
  par$logRWsd <- 0
  #par$tHatch <- 0.1
  
  f <- function(par){
    getAll(par,dat)
    L1 <- exp(c(logL1.1,logL1.rand))
    L2 <- exp(c(logL2.1,logL2.rand))
    L3 <- exp(c(logL3.1,logL3.rand))
    t_hatch <- as.numeric(t_hatches[unique(Species)])
    L_hatch <- as.numeric(L_hatches[unique(Species)])
    
    r <- (L3-L2)/(L2-L1)
    
    t1 <- as.numeric(t1s[unique(Species)])
    t3 <- as.numeric(t3s[unique(Species)])
    #t2 <- (t1+t3)/2
    SD_m <- exp(logSD_m)
    SD_c <- exp(logSD_c)
    SD_age <- exp(logSDage)
    SD_RW <- exp(logRWsd)
    
    VbgrF <- expression(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(age-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1))
    diffVBGR <- D(VbgrF,"age")
    
    cohorts <- sort(unique(cohort))
    
    df_mean <- matrix(NA,ncol=4)
    ret <- 0
    
    for(i in 1:length(cohorts)){
      cohort.idx <- which(cohort==cohorts[i])
      times.i <- sort(unique(t[cohort.idx]))
      if(i>1){
        ret <- ret-dnorm(log(L1[i]),mean=log(L1[i-1]),sd=SD_RW,log=TRUE)
        ret <- ret-dnorm(log(L2[i]),mean=log(L2[i-1]),sd=SD_RW,log=TRUE)
        ret <- ret-dnorm(log(L3[i]),mean=log(L3[i-1]),sd=SD_RW,log=TRUE)
      }
      
      for(j in 1:length(times.i)){
        age.idx <- which((t[cohort.idx])==times.i[j])
        #mean_idx <- which(t_mean==times.i[j] & cohort_mean==cohorts[i])
        age.j <- times.i[j]-t_hatch
        #      mean_size <- L_mean[mean_idx]
        #weighted mean
        mean_size <- sum(Length[cohort.idx][age.idx]*N[cohort.idx][age.idx])/sum(N[cohort.idx][age.idx])
        #weighted std. dev.
        mean_sd <- sqrt(sum(N[cohort.idx][age.idx]*(Length[cohort.idx][age.idx]-mean_size)^2)/sum(N[cohort.idx][age.idx]))
        
        sd_fac <- ifelse(length(age.idx)>1,1,0)
        
        wgts <- c((sum(N[cohort.idx][age.idx])),sqrt(sum(N[cohort.idx][age.idx])),1)[weights]
        
        weight <- N[cohort.idx][age.idx]/wgts[2]*sd_fac
        
        weight_m <- sum(N[cohort.idx][age.idx])/wgts[1]# or /sqrt(sum(mul[cohort.idx][age.idx]))
        
        dL_dt <- eval(diffVBGR,list(age=age.j,L1=L1[i],L2=L2[i],L3=L3[i],
                                    t1=t1,t3=t3))
        
        L.pred <- eval(VbgrF,list(age=age.j,L1=L1[i],L2=L2[i],L3=L3[i],
                                  t1=t1,t3=t3))
        
        LSD <- sqrt((SD_age*dL_dt)^2+(L.pred*SD_c)^2)
        
        ret <- ret-dnorm(log(mean_size),mean=log(L.pred),sd=SD_m,log=TRUE)*weight_m # fit mean length
        if(length(age.idx)<2) next
        
        ret <- ret-sum(dnorm(Length[cohort.idx][age.idx],
                             mean=mean_size,sd=LSD,
                             log=TRUE)*weight)
        
        df_mean <- rbind(df_mean,c(round(cohorts[i]),times.i[j],mean_size,mean_sd))
      }
      L.Hatch <- eval(VbgrF,list(age=0,L1=L1[i],L2=L2[i],L3=L3[i],
                                 t1=t1,t3=t3))
      ret <- ret-dnorm(L.Hatch,L_hatch,0.01,log=TRUE)#*10^4 play with sd here if the model doesn't converge
      
    }
    df_mean <- df_mean[-1,]
    
    ADREPORT(L1)
    ADREPORT(L2)
    ADREPORT(L3)

    ADREPORT(df_mean)
    
    ret
  }
  obj <- MakeADFun(f,par,random=c("logL1.rand","logL2.rand","logL3.rand"))
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr <- sdreport(obj)
  return(sdr)
}