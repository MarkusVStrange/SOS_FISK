########################
# a function to plot a VBGR fit from a certain species
# 
#' @param dat = a list containing: cohort, t (Age+time of year), N, Length, Species
#' @param sdr = an sdreport element from RTMP
#' @param version = a string: currently "together" or separate, which refers to cohorts
#' 
#' @return 
#' a plot with the specified cohort settings 
#########################

plotVBGR <- function(dat,sdr,version="separate"){
  if(version=="separate"){
    L1 <- as.list(sdr,"Est",report=TRUE)$L1
    L2 <- as.list(sdr,"Est",report=TRUE)$L2
    L3 <- as.list(sdr,"Est",report=TRUE)$L3
    dat_means <- as.list(sdr,"Est",report=TRUE)$df_mean
    
    cohorts <- sort(unique(dat$cohort))
    t <- (0:((max(dat$Age)+1)*100))/100
    
    t_hatches <- c(cod=61/365,flounder=((75+197)/2)/365,plaice=0.03972603,dab=0.4561644)
    
    t1s <- c(cod=1,flounder=2,plaice=2,dab=2)
    t3s <- c(cod=5,flounder=10,plaice=10,dab=10)
    
    t_hatch <- as.numeric(t_hatches[unique(dat$Species)])
    SD_age <- exp(summary(sdr)[rownames(summary(sdr))=="logSDage",1])
    SD_c <- exp(summary(sdr)[rownames(summary(sdr))=="logSD_c",1])
    par(mfrow=c(1,1))
    for(i in 1:length(cohorts)){
      VbgrF <- expression(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(age-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1))
      diffVBGR <- D(VbgrF,"age")
      
      L <- eval(VbgrF,list(age=t-t_hatch,L1=L1[i],L2=L2[i],L3=L3[i],
                           t1=as.numeric(t1s[unique(dat$Species)]),
                           t3=as.numeric(t3s[unique(dat$Species)])))
      dL_dage <- eval(diffVBGR,list(age=t-t_hatch,L1=L1[i],L2=L2[i],L3=L3[i],
                                    t1=as.numeric(t1s[unique(dat$Species)]),
                                    t3=as.numeric(t3s[unique(dat$Species)])))
      
      LSD <- sqrt((SD_age*dL_dage)^2+(L*SD_c)^2)
      m.idx <- which(dat_means[,1]==cohorts[i])
      
      dL_dageData <- eval(diffVBGR,list(age=dat_means[m.idx,2],L1=L1[i],L2=L2[i],L3=L3[i],
                                        t1=as.numeric(t1s[unique(dat$Species)]),
                                        t3=as.numeric(t3s[unique(dat$Species)])))
      LData <- eval(VbgrF,list(age=dat_means[m.idx,2],L1=L1[i],L2=L2[i],L3=L3[i],
                               t1=as.numeric(t1s[unique(dat$Species)]),
                               t3=as.numeric(t3s[unique(dat$Species)])))
      dat_sd <- sqrt((SD_age*dL_dageData)^2+(LData*SD_c)^2)
      
      plot(t,(L) , lwd =3 ,type='l',
           ylim=c(0,max(dat$Length)*1.1),col="black",ylab="Length [cm]",xlab="Age [years]",
           main=paste(unique(dat$Species),"cohort",cohorts[i]))
      lines(t,(L -LSD), lty = "dotted" , lwd =2.5,col="orange")
      lines(t,(L +LSD), lty = "dotted" , lwd =2.5,col="orange")  
      points(dat_means[m.idx,2],dat_means[m.idx,3],pch=19,cex=2)
      points(dat_means[m.idx,2],dat_means[m.idx,3]-dat_sd,pch=19,cex=1,col="orange")
      points(dat_means[m.idx,2],dat_means[m.idx,3]+dat_sd,pch=19,cex=1,col="orange")
      points(dat_means[m.idx,2],dat_means[m.idx,3]-dat_means[m.idx,4],pch=19,cex=1,col="red")
      points(dat_means[m.idx,2],dat_means[m.idx,3]+dat_means[m.idx,4],pch=19,cex=1,col="red")
      lines(c(-10,100),c(0,0),lty="dashed")
      lines(c(t_hatch,t_hatch),c(-10,200),lty="dashed")
      lines(c(0,0),c(-10,200),lty="dotted")
      #lines(c(1,1),c(-10,200),lty="dotted")
      #lines(c(0.25,0.25),c(-10,200),lty="dashed")
      lines(c(-10,100),c(0,0),lwd=1)

      legend(max(dat$Age)*0.7,max(dat$Length)*0.3,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))
      
      #readline() # press enter to continue
    }
  }
  
  if(version=="together"){
    cohorts <- sort(unique(dat$cohort))
    L1 <- as.list(sdr,"Est",report=TRUE)$L1
    L2 <- as.list(sdr,"Est",report=TRUE)$L2
    L3 <- as.list(sdr,"Est",report=TRUE)$L3
    
    colors <- colorRampPalette(c("red", "blue"))(length(cohorts))
    t <- (0:((max(dat$Age)+1)*100))/100
    
    t_hatches <- c(cod=61/365,flounder=((75+197)/2)/365,plaice=0.03972603,dab=0.4561644)
    t1s <- c(cod=1,flounder=2,plaice=2,dab=2)
    t3s <- c(cod=5,flounder=10,plaice=10,dab=10)
    
    t_hatch <- as.numeric(t_hatches[unique(dat$Species)])
    SD_age <- exp(summary(sdr)[rownames(summary(sdr))=="logSDage",1])
    SD_c <- exp(summary(sdr)[rownames(summary(sdr))=="logSD_c",1])
    par(mfrow=c(1,1))
    for(i in 1:length(cohorts)){
      VbgrF <- expression(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(age-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1))
      diffVBGR <- D(VbgrF,"age")
      
      L <- eval(VbgrF,list(age=t-t_hatch,L1=L1[i],L2=L2[i],L3=L3[i],
                           t1=as.numeric(t1s[unique(dat$Species)]),
                           t3=as.numeric(t3s[unique(dat$Species)])))
      #dL_dage <- eval(diffVBGR,list(age=t-t_hatch,L1=L1[i],L2=L2[i],L3=L3[i],
      #                          t1=1,t3=5))
      #
      #LSD <- sqrt((SD_age*dL_dage)^2+(L.pred*SD_c)^2)
      
      if(i==1){
        plot(t,L, lwd =2,type='l',ylim=c(0,max(dat$Length)*1.2),col=colors[i],
             main=paste(unique(dat$Species),"cohorts",min(dat$cohort),"-",max(dat$cohort)))
        lines(c(-10,100),c(0,0),lty="dashed")
        lines(c(t_hatch,t_hatch),c(-10,200),lty="dashed")
      }
      if(i>1){
        lines(t,L, lwd =2,col=colors[i])
      }
      
      #readline() # press enter to continue
    }
  }
  
  
  if(version=="data"){
    cohorts <- sort(unique(dat$cohort))
    dat_means <- as.data.frame(as.list(sdr,"Est",report=TRUE)$df_mean)
    names(dat_means) <- c("year","t","c_mean","sd")
    dat_means$t <- floor(dat_means$t*10/5)*5/10
    
    
    colors <- colorRampPalette(c("red", "blue"))(length(cohorts))
    
    m_m <- aggregate(c_mean~t,data=dat_means,FUN = mean)
    names(m_m) <- c("t","mean")
    
    dat_means <- left_join(dat_means,m_m)
    par(mfrow=c(1,1))
    for(i in 1:length(cohorts)){
      m.idx <- which(dat_means[,1]==cohorts[i])

      if(i==1){
        plot(dat_means$t[m.idx]+3/8,log(dat_means$c_mean[m.idx]/dat_means$mean[m.idx]),pch=19,cex=2,col=colors[i],
               main=paste(unique(dat$Species),"cohorts",min(dat$cohort),"-",max(dat$cohort)),
             ylim=c(min(log(dat_means$c_mean/dat_means$mean))*1.1,max(log(dat_means$c_mean/dat_means$mean))*1.1),
             xlim=c(0,max(dat_means$t)+1),xlab="time [years]",ylab="log(cohort mean : overall mean)")

      }
      if(i>1){
        points(dat_means$t[m.idx]+3/8,log(dat_means$c_mean[m.idx]/dat_means$mean[m.idx]),pch=19,cex=2,col=colors[i])
      }
      
      #readline() # press enter to continue
    }
  }
  
  if(version=="coefficients"){
    cohorts <- sort(unique(dat$cohort))
    L1 <- as.list(sdr,"Est",report=TRUE)$L1
    L2 <- as.list(sdr,"Est",report=TRUE)$L2
    L3 <- as.list(sdr,"Est",report=TRUE)$L3
    
    colors <- colorRampPalette(c("red", "blue"))(length(cohorts))
    par(mfrow=c(1,3))
    plot(cohorts,L1,pch=19,cex=1,col=colors,
         main="L1",
         xlab="cohort [years]",ylab="cm")
    
    plot(cohorts,L2,pch=19,cex=1,col=colors,
         main="L2",
         xlab="cohort [years]",ylab="cm")
    plot(cohorts,L3,pch=19,cex=1,col=colors,
         main="L3",
         xlab="cohort [years]",ylab="cm")
  }
}