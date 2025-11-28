#remotes::install_github("DTUAqua/DATRAS/DATRAS")
#remotes::install_github("casperwberg/surveyIndex/surveyIndex")

library(surveyIndex)
library(sf)
library(reshape2)
library(tidyverse)

data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="")

#dAll <- readRDS(paste(data_wd,"DatrasExchange.R",sep=""))
dAll <- readRDS("datras.R")
WB <- st_read(paste(data_wd,"shapefiles/ICES_areas.shp",sep=""))
plot(WB %>% filter(ICES_SUB %in% 22:24))
WB <- as(WB,"Spatial")

dAll <- addSpatialData(dAll,WB)

# Cod ALK Q1
#####
mc.cores<-1; library(parallel)
years <- 1991:2025

codQ1 <- data.frame(LngtCm=NA,Age=NA,N=NA,Year=NA,Quarter=NA,qday=NA,Species="O",FishingDur=NA)

for(i in years){
  d<-subset(dAll, Species=="Gadus morhua",Quarter==1,ICES_SUB %in% 22:24,Year==i)
  d<-addSpectrum(d,by=1)
  ## get idea about number of age groups to include
  agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])

  ages <- sort(unique(d[["CA"]]$Age))
  #if(0 %in% ages) ages <- ages[-1]
  
  if(sum(diff(ages))!=(length(ages)-1)){
    ages <- ages[1:which(diff(ages)>1)]
  }
  d <- subset(d,Age %in% ages)
  d<-addSpectrum(d,by=1)

  d.ysplit <- split(d, d$Year)
  ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
                varCof=FALSE,maxK=50,mc.cores=mc.cores)
  
  hh <- d[["HH"]]
  hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
  hl <- d[["HL"]]
  hl <- left_join(hl,hh %>% dplyr::select(haul.id,Jday))
  
  
  dat.i = attr(ALK[[1]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  fishing.dur <- 0
  N.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  j.day <- matrix(0,nrow=n_hauls,ncol=2)
  for(j in 1:n_hauls){
    dada <- hl %>% filter(haul.id == as.character(unique(dat.i$haul.id)[j]) & !is.na(LngtCm | HLNoAtLngt)) %>% 
      dplyr::select(LngtCm,HLNoAtLngt,Jday,HaulDur)
    if(length(dada$LngtCm)<1) next
    
    alk.j <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    if(j==1) plotALKfit(ALK[[1]],row=j,main=i)
    
    N.i[match(floor(dada$LngtCm),cm.breaks),,j] <- (alk.j[match(floor(dada$LngtCm),cm.breaks),])*dada$HLNoAtLngt
    j.day[j,] <- c(unique(dada$Jday)/365,sum(dada$HLNoAtLngt))
    fishing.dur <- fishing.dur+dada$HaulDur[1]
  }
  a <- melt(apply(N.i,c(1,2),FUN=sum))
  names(a) <- c("LngtCm","Age","N")
  a$LngtCm <- cm.breaks[a$LngtCm]
  a$Age <- ages[a$Age]
  a$Year <- as.numeric(i)
  a$Quarter <- 1
  a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
  a$Species <- "cod"
  a$FishingDur <- fishing.dur
  
  codQ1 <- rbind(codQ1,a)
  print(i)
}
codQ1 <- codQ1[-1,]

saveRDS(codQ1, file="N_cod1.RData")
readRDS("N_cod1.RData")

#####

# Cod ALK Q4
#####
mc.cores<-1; library(parallel)
years <- 1991:2024

codQ4 <- data.frame(LngtCm=NA,Age=NA,N=NA,Year=NA,Quarter=NA,qday=NA,Species="O",FishingDur=NA)

for(i in years){
  d<- subset(dAll, Species=="Gadus morhua",Quarter==4,ICES_SUB %in% 22:24,Year==i)
  d<- addSpectrum(d,by=1)
  if(i %in% c(2015,2018)){
    d <- subset(d,Country!="DE")
    d<-addSpectrum(d,by=1)
  }

  ## get idea about number of age groups to include
  agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
  
  ages <- sort(unique(d[["CA"]]$Age))

  if(sum(diff(ages))!=(length(ages)-1)){
    ages <- ages[1:which(diff(ages)>1)]
  }
  
  d <- subset(d,Age %in% ages)
  d<-addSpectrum(d,by=1)
  
  d.ysplit <- split(d, d$Year)
  ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
                varCof=FALSE,maxK=50,mc.cores=mc.cores)
  
  hh <- d[["HH"]]
  hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
  hl <- d[["HL"]]
  hl <- left_join(hl,hh %>% dplyr::select(haul.id,Jday))

  dat.i = attr(ALK[[1]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  fishing.dur <- 0
  N.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  j.day <- matrix(0,nrow=n_hauls,ncol=2)
  for(j in 1:n_hauls){
    dada <- hl %>% filter(haul.id == as.character(unique(dat.i$haul.id)[j]) & !is.na(LngtCm | HLNoAtLngt)) %>% 
      dplyr::select(LngtCm,HLNoAtLngt,Jday,HaulDur)

    if(length(dada$LngtCm)<1) next
    
    alk.j <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    if(j==1) plotALKfit(ALK[[1]],row=j,main=i)
    
    N.i[match(floor(dada$LngtCm),cm.breaks),,j] <- (alk.j[match(floor(dada$LngtCm),cm.breaks),])*dada$HLNoAtLngt
    j.day[j,] <- c(unique(dada$Jday)/365,sum(dada$HLNoAtLngt))
    fishing.dur <- fishing.dur+dada$HaulDur[1]
  }
  a <- melt(apply(N.i,c(1,2),FUN=sum))
  names(a) <- c("LngtCm","Age","N")
  a$LngtCm <- cm.breaks[a$LngtCm]
  a$Age <- ages[a$Age]
  a$Year <- as.numeric(i)
  a$Quarter <- 4
  a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
  a$Species <- "cod"
  a$FishingDur <- fishing.dur
  
  codQ4 <- rbind(codQ4,a)
  print(i)
}
codQ4 <- codQ4[-1,]

saveRDS(codQ4, file="N_cod4.RData")
readRDS("N_cod4.RData")

#####


#rm(list=setdiff(ls(),c('dAll')))
# Flounder ALK Q1
#####
mc.cores<-1; library(parallel)
years <- c(1993,1991,1992,1994:2025)

flounderQ1 <- data.frame(LngtCm=NA,Age=NA,N=NA,Year=NA,Quarter=NA,qday=NA,Species="O")

for(i in years){
  d<-subset(dAll, Species=="Platichthys flesus",Quarter==1,ICES_SUB %in% 22:24,Year==i,LngtCm<100)
  d<-addSpectrum(d,by=1)
  ## get idea about number of age groups to include
  agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
  
  if(length(agetab)<1){
    hh <- d[["HH"]]
    hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
    hl <- d[["HL"]]
    hl <- left_join(hl,hh %>% select(haul.id,Jday))
    hl <- hl %>% filter(LngtCm<70)
    
    dat.i = attr(ALK[[1]], "data")
    cm.breaks = attr(dat.i, "cm.breaks")
    n_hauls <- length(unique(dat.i$haul.id))
    alk.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
    for(j in 1:n_hauls){
      alk.i[,,j] <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    }
    alk_mean <- apply(alk.i,c(1,2),FUN=mean)
    matplot(cm.breaks[-length(cm.breaks)], alk_mean, lwd = 2, type = 'l',main=i)    
    hl.ag <- aggregate(HLNoAtLngt~LngtCm,data=hl,FUN=sum)
    
    N.i <- alk_mean[match(sort(unique(hl$LngtCm)),cm.breaks),]*hl.ag$HLNoAtLngt
    
    a <- melt(apply(N.i,c(1,2),FUN=sum))
    names(a) <- c("LngtCm","Age","N")
    a$LngtCm <- hl.ag$LngtCm[a$LngtCm]
    a$Year <- as.numeric(i)
    a$Quarter <- 1
    a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
    a$Species <- "flounder"
    
    flounderQ1 <- rbind(flounderQ1,a)
    print(i)
    next
  }
  
  ages <- sort(unique(d[["CA"]]$Age))
  #if(0 %in% ages) ages <- ages[-1]
  
  if(sum(diff(ages))!=(length(ages)-1)){
    ages <- ages[1:which(diff(ages)>1)]
  }
  d <- subset(d,Age %in% ages)
  d<-addSpectrum(d,by=1)

  d.ysplit <- split(d, d$Year)
  ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
                varCof=FALSE,maxK=50,mc.cores=mc.cores)
  
  hh <- d[["HH"]]
  hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
  hl <- d[["HL"]]
  hl <- left_join(hl,hh %>% select(haul.id,Jday))
  hl <- hl %>% filter(LngtCm<100)
  
  dat.i = attr(ALK[[1]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  N.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  j.day <- matrix(0,nrow=n_hauls,ncol=2)
  for(j in 1:n_hauls){
    dada <- hl %>% filter(haul.id == as.character(unique(dat.i$haul.id)[j]) & !is.na(LngtCm | HLNoAtLngt)) %>% 
      select(LngtCm,HLNoAtLngt,Jday)
    if(length(dada$LngtCm)<1) next

    alk.j <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    if(j==1) plotALKfit(ALK[[1]],row=j,main=i)
    
    N.i[match(floor(dada$LngtCm),cm.breaks),,j] <- (alk.j[match(floor(dada$LngtCm),cm.breaks),])*dada$HLNoAtLngt
    j.day[j,] <- c(unique(dada$Jday)/365,sum(dada$HLNoAtLngt))
  }
  a <- melt(apply(N.i,c(1,2),FUN=sum))
  names(a) <- c("LngtCm","Age","N")
  a$LngtCm <- cm.breaks[a$LngtCm]
  a$Age <- ages[a$Age]
  a$Year <- as.numeric(i)
  a$Quarter <- 1
  a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
  a$Species <- "flounder"
  
  flounderQ1 <- rbind(flounderQ1,a)
  print(i)
}
flounderQ1 <- flounderQ1[-1,]

saveRDS(flounderQ1, file="N_flounder1.RData")
readRDS("N_flounder1.RData")

#####



# Flounder ALK Q4
#####
mc.cores<-1; library(parallel)
years <- 1991:2024

flounderQ4 <- data.frame(LngtCm=NA,Age=NA,N=NA,Year=NA,Quarter=NA,qday=NA,Species="O")

for(i in years){
  d<-subset(dAll, Species=="Platichthys flesus",Quarter==4,ICES_SUB %in% 22:24, Year==i,LngtCm<100)
  d<-addSpectrum(d,by=1)
  ## get idea about number of age groups to include
  agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])

  ages <- sort(unique(d[["CA"]]$Age))

  if(sum(diff(ages))!=(length(ages)-1)){
    ages <- ages[1:which(diff(ages)>1)]
    d <- subset(d,Age %in% ages)
    d<-addSpectrum(d,by=1)
  }
  
  d.ysplit <- split(d, d$Year)
  ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
                varCof=FALSE,maxK=50,mc.cores=mc.cores)
  
  hh <- d[["HH"]]
  hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
  hl <- d[["HL"]]
  hl <- left_join(hl,hh %>% select(haul.id,Jday))

  dat.i = attr(ALK[[1]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  N.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  j.day <- matrix(0,nrow=n_hauls,ncol=2)
  for(j in 1:n_hauls){
    dada <- hl %>% filter(haul.id == as.character(unique(dat.i$haul.id)[j]) & !is.na(LngtCm | HLNoAtLngt)) %>% 
      select(LngtCm,HLNoAtLngt,Jday)
    if(length(dada$LngtCm)<1) next
    
    alk.j <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    if(j==1) plotALKfit(ALK[[1]],row=j,main=i)
    
    N.i[match(floor(dada$LngtCm),cm.breaks),,j] <- (alk.j[match(floor(dada$LngtCm),cm.breaks),])*dada$HLNoAtLngt
    j.day[j,] <- c(unique(dada$Jday)/365,sum(dada$HLNoAtLngt))
  }
  a <- melt(apply(N.i,c(1,2),FUN=sum))
  names(a) <- c("LngtCm","Age","N")
  a$LngtCm <- cm.breaks[a$LngtCm]
  a$Age <- ages[a$Age]
  a$Year <- as.numeric(i)
  a$Quarter <- 4
  a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
  a$Species <- "flounder"
  
  flounderQ4 <- rbind(flounderQ4,a)
  print(i)
}
flounderQ4 <- flounderQ4[-1,]

saveRDS(flounderQ4, file="N_flounder4.RData")
t <- readRDS("N_flounder4.RData")

#####

# Plaice ALK Q1
#####
mc.cores<-1; library(parallel)
years <- c(1994,1991:1993,1995:2025)

plaiceQ1 <- data.frame(LngtCm=NA,Age=NA,N=NA,Year=NA,Quarter=NA,qday=NA,Species="O")

for(i in years){
  d<-subset(dAll, Species=="Pleuronectes platessa",Quarter==1,ICES_SUB %in% 22:24,Year==i,LngtCm<100)
  d<-addSpectrum(d,by=1)
  ## get idea about number of age groups to include
  agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
  
  if(length(agetab)<1){
    d <- subset(d,LngtCm %in% (cm.breaks-1))
    d<-addSpectrum(d,by=1)
    hh <- d[["HH"]]
    hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
    hl <- d[["HL"]]
    hl <- left_join(hl,hh %>% select(haul.id,Jday))
    hl <- hl %>% filter(LngtCm<100)
    
    dat.i = attr(ALK[[1]], "data")
    cm.breaks = attr(dat.i, "cm.breaks")
    n_hauls <- length(unique(dat.i$haul.id))
    alk.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
    for(j in 1:n_hauls){
      alk.i[,,j] <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    }
    alk_mean <- apply(alk.i,c(1,2),FUN=mean)
    matplot(cm.breaks[-length(cm.breaks)], alk_mean, lwd = 2, type = 'l',main=i)    
    hl.ag <- aggregate(HLNoAtLngt~LngtCm,data=hl,FUN=sum)
    
    N.i <- alk_mean[match(sort(unique(hl$LngtCm)),cm.breaks),]*hl.ag$HLNoAtLngt
    
    a <- melt(apply(N.i,c(1,2),FUN=sum))
    names(a) <- c("LngtCm","Age","N")
    a$LngtCm <- hl.ag$LngtCm[a$LngtCm]
    a$Year <- as.numeric(i)
    a$Quarter <- 1
    a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
    a$Species <- "plaice"
    
    plaiceQ1 <- rbind(plaiceQ1,a)
    print(i)
    next
  }
  
  ages <- sort(unique(d[["CA"]]$Age))
  #if(0 %in% ages) ages <- ages[-1]
  
  if(sum(diff(ages))!=(length(ages)-1)){
    ages <- ages[1:which(diff(ages)>1)]
  }
  d <- subset(d,Age %in% ages)
  d<-addSpectrum(d,by=1)
  
  d.ysplit <- split(d, d$Year)
  ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
                varCof=FALSE,maxK=50,mc.cores=mc.cores)
  
  hh <- d[["HH"]]
  hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
  hl <- d[["HL"]]
  hl <- left_join(hl,hh %>% select(haul.id,Jday))
  hl <- hl %>% filter(LngtCm<100)
  
  dat.i = attr(ALK[[1]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  N.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  j.day <- matrix(0,nrow=n_hauls,ncol=2)
  for(j in 1:n_hauls){
    dada <- hl %>% filter(haul.id == as.character(unique(dat.i$haul.id)[j]) & !is.na(LngtCm | HLNoAtLngt)) %>% 
      select(LngtCm,HLNoAtLngt,Jday)
    if(length(dada$LngtCm)<1) next
    
    alk.j <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    if(j==1) plotALKfit(ALK[[1]],row=j,main=i)
    
    N.i[match(floor(dada$LngtCm),cm.breaks),,j] <- (alk.j[match(floor(dada$LngtCm),cm.breaks),])*dada$HLNoAtLngt
    j.day[j,] <- c(unique(dada$Jday)/365,sum(dada$HLNoAtLngt))
  }
  a <- melt(apply(N.i,c(1,2),FUN=sum))
  names(a) <- c("LngtCm","Age","N")
  a$LngtCm <- cm.breaks[a$LngtCm]
  a$Age <- ages[a$Age]
  a$Year <- as.numeric(i)
  a$Quarter <- 1
  a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
  a$Species <- "plaice"
  
  plaiceQ1 <- rbind(plaiceQ1,a)
  print(i)
}
plaiceQ1 <- plaiceQ1[-1,]

saveRDS(plaiceQ1, file="N_plaice1.RData")
readRDS("N_plaice1.RData")

#####

# Plaice ALK Q4
#####
mc.cores<-1; library(parallel)
years <- c(1994,1991:1993,1995:2024)

plaiceQ4 <- data.frame(LngtCm=NA,Age=NA,N=NA,Year=NA,Quarter=NA,qday=NA,Species="O")

for(i in years){
  d<-subset(dAll, Species=="Pleuronectes platessa",Quarter==4,ICES_SUB %in% 22:24,Year==i,LngtCm<100)
  d<-addSpectrum(d,by=1)
  ## get idea about number of age groups to include
  agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
  
  if(length(agetab)<1){
    d <- subset(d,LngtCm %in% (cm.breaks-1))
    d<-addSpectrum(d,by=1)
    hh <- d[["HH"]]
    hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
    hl <- d[["HL"]]
    hl <- left_join(hl,hh %>% select(haul.id,Jday))
    hl <- hl %>% filter(LngtCm<100)
    
    dat.i = attr(ALK[[1]], "data")
    cm.breaks = attr(dat.i, "cm.breaks")
    n_hauls <- length(unique(dat.i$haul.id))
    alk.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
    for(j in 1:n_hauls){
      alk.i[,,j] <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    }
    alk_mean <- apply(alk.i,c(1,2),FUN=mean)
    matplot(cm.breaks[-length(cm.breaks)], alk_mean, lwd = 2, type = 'l',main=i)    
    hl.ag <- aggregate(HLNoAtLngt~LngtCm,data=hl,FUN=sum)
    
    N.i <- alk_mean[match(sort(unique(hl$LngtCm)),cm.breaks),]*hl.ag$HLNoAtLngt
    
    a <- melt(apply(N.i,c(1,2),FUN=sum))
    names(a) <- c("LngtCm","Age","N")
    a$LngtCm <- hl.ag$LngtCm[a$LngtCm]
    a$Year <- as.numeric(i)
    a$Quarter <- 4
    a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
    a$Species <- "plaice"
    
    plaiceQ4 <- rbind(plaiceQ4,a)
    print(i)
    next
  }
  
  ages <- sort(unique(d[["CA"]]$Age))

  if(sum(diff(ages))!=(length(ages)-1)){
    ages <- ages[1:which(diff(ages)>1)]
  }
  d <- subset(d,Age %in% ages)
  d<-addSpectrum(d,by=1)
  
  d.ysplit <- split(d, d$Year)
  ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
                varCof=FALSE,maxK=50,mc.cores=mc.cores)
  
  hh <- d[["HH"]]
  hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
  hl <- d[["HL"]]
  hl <- left_join(hl,hh %>% select(haul.id,Jday))
  hl <- hl %>% filter(LngtCm<70)
  
  dat.i = attr(ALK[[1]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  N.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  j.day <- matrix(0,nrow=n_hauls,ncol=2)
  for(j in 1:n_hauls){
    dada <- hl %>% filter(haul.id == as.character(unique(dat.i$haul.id)[j]) & !is.na(LngtCm | HLNoAtLngt)) %>% 
      select(LngtCm,HLNoAtLngt,Jday)
    if(length(dada$LngtCm)<1) next
    
    alk.j <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    if(j==1) plotALKfit(ALK[[1]],row=j,main=i)
    
    N.i[match(floor(dada$LngtCm),cm.breaks),,j] <- (alk.j[match(floor(dada$LngtCm),cm.breaks),])*dada$HLNoAtLngt
    j.day[j,] <- c(unique(dada$Jday)/365,sum(dada$HLNoAtLngt))
  }
  a <- melt(apply(N.i,c(1,2),FUN=sum))
  names(a) <- c("LngtCm","Age","N")
  a$LngtCm <- cm.breaks[a$LngtCm]
  a$Age <- ages[a$Age]
  a$Year <- as.numeric(i)
  a$Quarter <- 4
  a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
  a$Species <- "plaice"
  
  plaiceQ4 <- rbind(plaiceQ4,a)
  print(i)
}
plaiceQ4 <- plaiceQ4[-1,]

saveRDS(plaiceQ4, file="N_plaice4.RData")
readRDS("N_plaice4.RData")

#####

# dab ALK Q1
#####
mc.cores<-1; library(parallel)
years <- c(1994,1991:1993,1995:2025)

dabQ1 <- data.frame(LngtCm=NA,Age=NA,N=NA,Year=NA,Quarter=NA,qday=NA,Species="O")

for(i in years){
  d<-subset(dAll, Species=="Limanda limanda",Quarter==1,ICES_SUB %in% 22:24,Year==i,LngtCm<100)
  d<-addSpectrum(d,by=1)
  ## get idea about number of age groups to include
  agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
  
  if(length(agetab)<1){
    d <- subset(d,LngtCm %in% (cm.breaks-1))
    d<-addSpectrum(d,by=1)
    hh <- d[["HH"]]
    hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
    hl <- d[["HL"]]
    hl <- left_join(hl,hh %>% select(haul.id,Jday))
    hl <- hl %>% filter(LngtCm<100)
    
    dat.i = attr(ALK[[1]], "data")
    cm.breaks = attr(dat.i, "cm.breaks")
    n_hauls <- length(unique(dat.i$haul.id))
    alk.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
    for(j in 1:n_hauls){
      alk.i[,,j] <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    }
    alk_mean <- apply(alk.i,c(1,2),FUN=mean)
    matplot(cm.breaks[-length(cm.breaks)], alk_mean, lwd = 2, type = 'l',main=i)    
    hl.ag <- aggregate(HLNoAtLngt~LngtCm,data=hl,FUN=sum)
    
    N.i <- alk_mean[match(sort(unique(hl$LngtCm)),cm.breaks),]*hl.ag$HLNoAtLngt
    
    a <- melt(apply(N.i,c(1,2),FUN=sum))
    names(a) <- c("LngtCm","Age","N")
    a$LngtCm <- hl.ag$LngtCm[a$LngtCm]
    a$Year <- as.numeric(i)
    a$Quarter <- 1
    a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
    a$Species <- "dab"
    
    dabQ1 <- rbind(dabQ1,a)
    print(i)
    next
  }
  
  ages <- sort(unique(d[["CA"]]$Age))
  #if(0 %in% ages) ages <- ages[-1]
  
  if(sum(diff(ages))!=(length(ages)-1)){
    ages <- ages[1:which(diff(ages)>1)]
    
  }
  d <- subset(d,Age %in% ages)
  d<-addSpectrum(d,by=1)
  
  d.ysplit <- split(d, d$Year)
  ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
                varCof=FALSE,maxK=50,mc.cores=mc.cores)
  
  hh <- d[["HH"]]
  hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
  hl <- d[["HL"]]
  hl <- left_join(hl,hh %>% select(haul.id,Jday))
  hl <- hl %>% filter(LngtCm<100)
  
  dat.i = attr(ALK[[1]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  N.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  j.day <- matrix(0,nrow=n_hauls,ncol=2)
  for(j in 1:n_hauls){
    dada <- hl %>% filter(haul.id == as.character(unique(dat.i$haul.id)[j]) & !is.na(LngtCm | HLNoAtLngt)) %>% 
      select(LngtCm,HLNoAtLngt,Jday)
    if(length(dada$LngtCm)<1) next
    
    alk.j <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    if(j==1) plotALKfit(ALK[[1]],row=j,main=i)
    
    N.i[match(floor(dada$LngtCm),cm.breaks),,j] <- (alk.j[match(floor(dada$LngtCm),cm.breaks),])*dada$HLNoAtLngt
    j.day[j,] <- c(unique(dada$Jday)/365,sum(dada$HLNoAtLngt))
  }
  a <- melt(apply(N.i,c(1,2),FUN=sum))
  names(a) <- c("LngtCm","Age","N")
  a$LngtCm <- cm.breaks[a$LngtCm]
  a$Age <- ages[a$Age]
  a$Year <- as.numeric(i)
  a$Quarter <- 1
  a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
  a$Species <- "dab"
  
  dabQ1 <- rbind(dabQ1,a)
  print(i)
}
dabQ1 <- dabQ1[-1,]

saveRDS(dabQ1, file="N_dab1.RData")
readRDS("N_dab1.RData")

#####

# dab ALK Q4
#####
mc.cores<-1; library(parallel)
years <- c(1994,1991:1993,1995:2024)

dabQ4 <- data.frame(LngtCm=NA,Age=NA,N=NA,Year=NA,Quarter=NA,qday=NA,Species="O")

for(i in years){
  d<-subset(dAll, Species=="Limanda limanda",Quarter==4,ICES_SUB %in% 22:24)
  d<-addSpectrum(d,by=1)
  ## get idea about number of age groups to include
  agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
  
  if(length(agetab)<1){
    d <- subset(d,LngtCm %in% (cm.breaks-1))
    d<-addSpectrum(d,by=1)
    hh <- d[["HH"]]
    hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
    hl <- d[["HL"]]
    hl <- left_join(hl,hh %>% select(haul.id,Jday))
    hl <- hl %>% filter(LngtCm<100)
    
    dat.i = attr(ALK[[1]], "data")
    cm.breaks = attr(dat.i, "cm.breaks")
    n_hauls <- length(unique(dat.i$haul.id))
    alk.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
    for(j in 1:n_hauls){
      alk.i[,,j] <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    }
    alk_mean <- apply(alk.i,c(1,2),FUN=mean)
    matplot(cm.breaks[-length(cm.breaks)], alk_mean, lwd = 2, type = 'l',main=i)    
    hl.ag <- aggregate(HLNoAtLngt~LngtCm,data=hl,FUN=sum)
    
    N.i <- alk_mean[match(sort(unique(hl$LngtCm)),cm.breaks),]*hl.ag$HLNoAtLngt
    
    a <- melt(apply(N.i,c(1,2),FUN=sum))
    names(a) <- c("LngtCm","Age","N")
    a$LngtCm <- hl.ag$LngtCm[a$LngtCm]
    a$Year <- as.numeric(i)
    a$Quarter <- 4
    a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
    a$Species <- "dab"
    
    dabQ4 <- rbind(dabQ4,a)
    print(i)
    next
  }
  
  ages <- sort(unique(d[["CA"]]$Age))

  if(sum(diff(ages))!=(length(ages)-1)){
    ages <- ages[1:which(diff(ages)>1)]
    
  }
  d <- subset(d,Age %in% ages)
  d<-addSpectrum(d,by=1)
  
  d.ysplit <- split(d, d$Year)
  ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
                varCof=FALSE,maxK=50,mc.cores=mc.cores)
  
  hh <- d[["HH"]]
  hh$Jday <- yday(as.POSIXct(paste(hh$Year,hh$Month,hh$Day,sep="-")))
  hl <- d[["HL"]]
  hl <- left_join(hl,hh %>% select(haul.id,Jday))
  hl <- hl %>% filter(LngtCm<100)
  
  dat.i = attr(ALK[[1]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  N.i <- array(0,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  j.day <- matrix(0,nrow=n_hauls,ncol=2)
  p.count <- 1
  for(j in 1:n_hauls){
    dada <- hl %>% filter(haul.id == as.character(unique(dat.i$haul.id)[j]) & !is.na(LngtCm | HLNoAtLngt)) %>% 
      select(LngtCm,HLNoAtLngt,Jday)
    if(length(dada$LngtCm)<1) next
    
    alk.j <- NageByHaul(j, ALK[[1]], returnALK = TRUE)
    if(p.count==1) plotALKfit(ALK[[1]],row=j,main=i+1)
    p.count <- p.count+1
    N.i[match(floor(dada$LngtCm),cm.breaks),,j] <- (alk.j[match(floor(dada$LngtCm),cm.breaks),])*dada$HLNoAtLngt
    j.day[j,] <- c(unique(dada$Jday)/365,sum(dada$HLNoAtLngt))
  }
  a <- melt(apply(N.i,c(1,2),FUN=sum))
  names(a) <- c("LngtCm","Age","N")
  a$LngtCm <- cm.breaks[a$LngtCm]
  a$Age <- ages[a$Age]
  a$Year <- as.numeric(i)
  a$Quarter <- 4
  a$qday <- sum(j.day[,1]*j.day[,2])/sum(j.day[,2])
  a$Species <- "dab"
  
  dabQ4 <- rbind(dabQ4,a)
  print(i)
}
dabQ4 <- dabQ4[-1,]

saveRDS(dabQ4, file="N_dab4.RData")
readRDS("N_dab4.RData")

#####


# Herring ALK Q1 - Not updated
#####
mc.cores<-1; library(parallel)
d<-subset(dAll, Species=="Clupea harengus",Quarter==1,ICES_SUB %in% 22:24)
#dAll<-NULL; gc(); ## lose dAll because it takes up a lot of memory
d<-addSpectrum(d,by=1)
## get idea about number of age groups to include
agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
agetab.df<-as.data.frame(agetab)
ages<-1:10
## require at least 1 aged individual in each year
for(a in ages){
  if(any(agetab.df$Freq[agetab.df$Age==a]<1))
    d<-fixAgeGroup(d,age=a,fun=ifelse(a==min(ages),"min","mean"))
}
d<-subset(d,Age>=min(ages))

d.ysplit <- split(d, d$Year)
ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
              varCof=FALSE,maxK=50,mc.cores=mc.cores)
years <- levels(d[[1]]$Year)
n_years <- length(years)
ca <- d[[1]]
hh <- d[[2]]
hl <- d[[3]]
hl.ag <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = hl,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haul.id,data = hl,FUN=mean)
hl.ag <- hl.ag %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.ag$cpue <- hl.ag$HLNoAtLngt/hl.ag$HaulDur

herQ1 <- data.frame(LngtCm=NA,Year=NA,Quarter=NA,Species="O",
                    HLNoAtLngt=NA,HaulDur=NA,cpue=NA,Age=NA,alk=NA)
ALK_Her1 <- list()
for(i in 1:n_years){
  dat.i = attr(ALK[[i]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  alk.i <- array(NA,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  for(j in 1:n_hauls){
    row=j
    alk.i[,,j] = NageByHaul(row, ALK[[i]], returnALK = TRUE)
  }
  ma <- apply(alk.i,c(1,2),FUN=mean)
  ALK_Her1[[years[i]]] <- apply(alk.i,c(1,2),FUN=mean)
  
  a <- melt(ALK_Her1[[years[i]]])
  names(a) <- c("LngtCm","Age","alk")
  a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])
  her.i <- hl.ag %>% filter(Year==years[i])
  
  h.i <- left_join(her.i,a)
  herQ1 <- rbind(herQ1,h.i)
  print(years[i])
}
herQ1 <- herQ1[-1,]
#saveRDS(ALK_Cod1, file="ALK_Cod1.RData")
rm(list=setdiff(ls(),c('dAll','codQ1','codQ4','herQ1','herQ4','flounderQ1','flounderQ4'
                       ,'plaiceQ1','plaiceQ4','dabQ1','dabQ4','WB_rects','hl_markus')))
#####

# Herring ALK Q4 - Not updated
#####
mc.cores<-1; library(parallel)
d<-subset(dAll, Species=="Clupea harengus",Quarter==4,
          StatRec %in% WB_rects)
years <- levels(d[[1]]$Year)  # Weird but it works
d<-subset(dAll, Species=="Clupea harengus",Quarter==4,
          StatRec %in% WB_rects,Year %in% years)
#dAll<-NULL; gc(); ## lose dAll because it takes up a lot of memory
d<-addSpectrum(d,by=1)
## get idea about number of age groups to include
agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
agetab.df<-as.data.frame(agetab)
ages<-1:10
## require at least 1 aged individual in each year
for(a in ages){
  if(any(agetab.df$Freq[agetab.df$Age==a]<1))
    d<-fixAgeGroup(d,age=a,fun=ifelse(a==min(ages),"min","mean"))
}
d<-subset(d,Age>=min(ages))

d.ysplit <- split(d, d$Year)
ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
              varCof=FALSE,maxK=50,mc.cores=mc.cores)
years <- levels(d[[1]]$Year)
n_years <- length(years)
ca <- d[[1]]
hh <- d[[2]]
hl <- d[[3]]
hl.ag <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = hl,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haul.id,data = hl,FUN=mean)
hl.ag <- hl.ag %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.ag$cpue <- hl.ag$HLNoAtLngt/hl.ag$HaulDur

herQ4 <- data.frame(LngtCm=NA,Year=NA,Quarter=NA,Species="O",
                    HLNoAtLngt=NA,HaulDur=NA,cpue=NA,Age=NA,alk=NA)
ALK_Her4 <- list()
for(i in 1:n_years){
  dat.i = attr(ALK[[i]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  alk.i <- array(NA,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  for(j in 1:n_hauls){
    row=j
    alk.i[,,j] = NageByHaul(row, ALK[[i]], returnALK = TRUE)
  }
  ma <- apply(alk.i,c(1,2),FUN=mean)
  ALK_Her4[[years[i]]] <- apply(alk.i,c(1,2),FUN=mean)
  
  a <- melt(ALK_Her4[[years[i]]])
  names(a) <- c("LngtCm","Age","alk")
  a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])
  her.i <- hl.ag %>% filter(Year==years[i])
  
  h.i <- left_join(her.i,a)
  herQ4 <- rbind(herQ4,h.i)
  print(years[i])
}
prior_years <-2021
alk1 <- ALK_Her4[["2020"]]
alk2 <- ALK_Her4[["2022"]]
alk <- (alk1+alk2)/2
matplot(cm.breaks[-length(cm.breaks)], alk, lwd = 2, type = 'l')
a <- melt(alk)
names(a) <- c("LngtCm","Age","alk")
a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])

dada <- hl_markus %>% filter(Year %in% prior_years & species=="herring" & Quarter==4)
names(dada)[c(3,5)] <- c("LngtCm","Species")
dada$Year <- factor(dada$Year)
hl.prior <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = dada,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haulID,data =  dada,FUN=mean)
hl.prior <- hl.prior %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.prior$cpue <- hl.prior$HLNoAtLngt/hl.prior$HaulDur
for(i in 1:length(prior_years)){
  her.i <- hl.prior %>% filter(Year==prior_years[i])
  h.i <- left_join(her.i,a)
  herQ4 <- rbind(herQ4,h.i)
  print(prior_years[i])
}
herQ4 <- herQ4[-1,]
#saveRDS(ALK_Cod1, file="ALK_Cod1.RData")
rm(list=setdiff(ls(),c('dAll','codQ1','codQ4','herQ1','herQ4','flounderQ1','flounderQ4'
                       ,'plaiceQ1','plaiceQ4','dabQ1','dabQ4','WB_rects','hl_markus')))
#####

# Flounder ALK QI
#####
mc.cores<-1; library(parallel)
d<-subset(dAll, Species=="Platichthys flesus",
          Year %in% factor(1992:2024),Quarter==1)
#dAll<-NULL; gc(); ## lose dAll because it takes up a lot of memory
d<-addSpectrum(d,by=1)
## get idea about number of age groups to include
agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
agetab.df<-as.data.frame(agetab)
ages<-1:15
## require at least 1 aged individual in each year
for(a in ages){
  if(any(agetab.df$Freq[agetab.df$Age==a]<1))
    d<-fixAgeGroup(d,age=a,fun=ifelse(a==min(ages),"min","mean"))
}
d<-subset(d,Age>=min(ages))

d.ysplit <- split(d, d$Year)
ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
              varCof=FALSE,maxK=50,mc.cores=mc.cores)
years <- levels(d[[1]]$Year)
n_years <- length(years)
ca <- d[[1]]
hh <- d[[2]]
hl <- d[[3]]
hl.ag <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = hl,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haul.id,data = hl,FUN=mean)
hl.ag <- hl.ag %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.ag$cpue <- hl.ag$HLNoAtLngt/hl.ag$HaulDur

flounderQ1 <- data.frame(LngtCm=NA,Year=NA,Quarter=NA,Species="O",
                    HLNoAtLngt=NA,HaulDur=NA,cpue=NA,Age=NA,alk=NA)
ALK_Flounder1 <- list()
for(i in 1:n_years){
  dat.i = attr(ALK[[i]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  alk.i <- array(NA,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  for(j in 1:n_hauls){
    row=j
    alk.i[,,j] = NageByHaul(row, ALK[[i]], returnALK = TRUE)
  }
  ma <- apply(alk.i,c(1,2),FUN=mean)
  ALK_Flounder1[[years[i]]] <- apply(alk.i,c(1,2),FUN=mean)
  
  a <- melt(ALK_Flounder1[[years[i]]])
  names(a) <- c("LngtCm","Age","alk")
  a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])
  flounder.i <- hl.ag %>% filter(Year==years[i])
  
  f.i <- left_join(flounder.i,a)
  flounderQ1 <- rbind(flounderQ1,f.i)
  print(years[i])
}
#prior_years <-1991
#alk1 <- ALK_Flounder1[["1992"]]
#alk2 <- ALK_Flounder1[["1993"]]
##alk3 <- ALK_Flounder1[["1994"]]
#alk <- (alk1+alk2+alk3)/3
#matplot(cm.breaks[-length(cm.breaks)], alk, lwd = 2, type = 'l')
#a <- melt(alk)
#names(a) <- c("LngtCm","Age","alk")
#a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])

#dada <- hl_markus %>% filter(Year %in% prior_years & species=="flounder" & Quarter==1)
#names(dada)[c(3,5)] <- c("LngtCm","Species")
#dada$Year <- factor(dada$Year)
#hl.prior <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = dada,FUN=sum)
#totHdur <- aggregate(HaulDur~Year+haulID,data =  dada,FUN=mean)
#hl.prior <- hl.prior %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
#hl.prior$cpue <- hl.prior$HLNoAtLngt/hl.prior$HaulDur
#for(i in 1:length(prior_years)){
#  flounder.i <- hl.prior %>% filter(Year==prior_years[i])
#  
#  f.i <- left_join(flounder.i,a)
#  flounderQ1 <- rbind(flounderQ1,f.i)
#  print(prior_years[i])
#}
flounderQ1 <- flounderQ1[-1,]
#saveRDS(ALK_Cod1, file="ALK_Cod1.RData")
rm(list=setdiff(ls(),c('dAll','codQ1','codQ4','herQ1','herQ4','flounderQ1','flounderQ4'
                       ,'plaiceQ1','plaiceQ4','dabQ1','dabQ4','WB_rects','hl_markus')))

#####

# Flounder ALK Q4
#####
mc.cores<-1; library(parallel)
d<-subset(dAll, Species=="Platichthys flesus",Quarter==4,StatRec %in% WB_rects)
#dAll<-NULL; gc(); ## lose dAll because it takes up a lot of memory
d<-addSpectrum(d,by=1)
## get idea about number of age groups to include
agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
agetab.df<-as.data.frame(agetab)
ages<-1:15
## require at least 1 aged individual in each year
for(a in ages){
  if(any(agetab.df$Freq[agetab.df$Age==a]<1))
    d<-fixAgeGroup(d,age=a,fun=ifelse(a==min(ages),"min","mean"))
}
d<-subset(d,Age>=min(ages))

d.ysplit <- split(d, d$Year)
ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
              varCof=FALSE,maxK=50,mc.cores=mc.cores)
years <- levels(d[[1]]$Year)
n_years <- length(years)
ca <- d[[1]]
hh <- d[[2]]
hl <- d[[3]]
hl.ag <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = hl,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haul.id,data = hl,FUN=mean)
hl.ag <- hl.ag %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.ag$cpue <- hl.ag$HLNoAtLngt/hl.ag$HaulDur

flounderQ4 <- data.frame(LngtCm=NA,Year=NA,Quarter=NA,Species="O",
                         HLNoAtLngt=NA,HaulDur=NA,cpue=NA,Age=NA,alk=NA)
ALK_Flounder4 <- list()
for(i in 1:n_years){
  dat.i = attr(ALK[[i]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  alk.i <- array(NA,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  for(j in 1:n_hauls){
    row=j
    alk.i[,,j] = NageByHaul(row, ALK[[i]], returnALK = TRUE)
  }
  ma <- apply(alk.i,c(1,2),FUN=mean)
  ALK_Flounder4[[years[i]]] <- apply(alk.i,c(1,2),FUN=mean)
  
  a <- melt(ALK_Flounder4[[years[i]]])
  names(a) <- c("LngtCm","Age","alk")
  a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])
  flounder.i <- hl.ag %>% filter(Year==years[i])
  
  f.i <- left_join(flounder.i,a)
  flounderQ4 <- rbind(flounderQ4,f.i)
  print(years[i])
}
flounderQ4 <- flounderQ4[-1,]
#saveRDS(ALK_Cod1, file="ALK_Cod1.RData")
rm(list=setdiff(ls(),c('dAll','codQ1','codQ4','herQ1','herQ4','flounderQ1','flounderQ4'
                       ,'plaiceQ1','plaiceQ4','dabQ1','dabQ4','WB_rects','hl_markus')))
#####

# Plaice ALK QI
#####
mc.cores<-1; library(parallel)
d<-subset(dAll, Species=="Pleuronectes platessa",Quarter==1,StatRec %in% WB_rects)
years <- levels(d[[1]]$Year)
d<-subset(dAll, Species=="Pleuronectes platessa",Quarter==1,StatRec %in% WB_rects,
          Year %in% years)
d<-addSpectrum(d,by=1)
## get idea about number of age groups to include
agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
agetab.df<-as.data.frame(agetab)
ages<-1:15
## require at least 1 aged individual in each year
for(a in ages){
  if(any(agetab.df$Freq[agetab.df$Age==a]<1))
    d<-fixAgeGroup(d,age=a,fun=ifelse(a==min(ages),"min","mean"))
}
d<-subset(d,Age>=min(ages))

d.ysplit <- split(d, d$Year)
ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
              varCof=FALSE,maxK=50,mc.cores=mc.cores)
years <- levels(d[[1]]$Year)
n_years <- length(years)
ca <- d[[1]]
hh <- d[[2]]
hl <- d[[3]]
hl.ag <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = hl,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haul.id,data = hl,FUN=mean)
hl.ag <- hl.ag %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.ag$cpue <- hl.ag$HLNoAtLngt/hl.ag$HaulDur

plaiceQ1 <- data.frame(LngtCm=NA,Year=NA,Quarter=NA,Species="O",
                         HLNoAtLngt=NA,HaulDur=NA,cpue=NA,Age=NA,alk=NA)
ALK_Plaice1 <- list()
for(i in 1:n_years){
  dat.i = attr(ALK[[i]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  alk.i <- array(NA,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  for(j in 1:n_hauls){
    row=j
    alk.i[,,j] = NageByHaul(row, ALK[[i]], returnALK = TRUE)
  }
  ma <- apply(alk.i,c(1,2),FUN=mean)
  ALK_Plaice1[[years[i]]] <- apply(alk.i,c(1,2),FUN=mean)
  
  a <- melt(ALK_Plaice1[[years[i]]])
  names(a) <- c("LngtCm","Age","alk")
  a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])
  plaice.i <- hl.ag %>% filter(Year==years[i])
  
  p.i <- left_join(plaice.i,a)
  plaiceQ1 <- rbind(plaiceQ1,p.i)
  print(years[i])
}

prior_years <-1991:1993
alk1 <- ALK_Plaice1[["1994"]]
alk2 <- ALK_Plaice1[["1995"]]
alk3 <- ALK_Plaice1[["1996"]]
alk <- (alk1+alk2+alk3)/3
matplot(cm.breaks[-length(cm.breaks)], alk, lwd = 2, type = 'l')
a <- melt(alk)
names(a) <- c("LngtCm","Age","alk")
a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])

dada <- hl_markus %>% filter(Year %in% prior_years & species=="plaice" & Quarter==1)
names(dada)[c(3,5)] <- c("LngtCm","Species")
dada$Year <- factor(dada$Year)
hl.prior <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = dada,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haulID,data =  dada,FUN=mean)
hl.prior <- hl.prior %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.prior$cpue <- hl.prior$HLNoAtLngt/hl.prior$HaulDur
for(i in 1:length(prior_years)){
  plaice.i <- hl.prior %>% filter(Year==prior_years[i])
  
  p.i <- left_join(plaice.i,a)
  plaiceQ1 <- rbind(plaiceQ1,p.i)
  print(prior_years[i])
}
plaiceQ1 <- plaiceQ1[-1,]
#saveRDS(ALK_Cod1, file="ALK_Cod1.RData")
rm(list=setdiff(ls(),c('dAll','codQ1','codQ4','herQ1','herQ4','flounderQ1','flounderQ4'
                       ,'plaiceQ1','plaiceQ4','dabQ1','dabQ1','WB_rects','hl_markus')))

#####

# Plaice ALK Q4
#####
mc.cores<-1; library(parallel)
d<-subset(dAll, Species=="Pleuronectes platessa",Quarter==4,StatRec %in% WB_rects)
years <- levels(d[[1]]$Year)
d<-subset(dAll, Species=="Pleuronectes platessa",Quarter==4,StatRec %in% WB_rects,
          Year %in% years)
d<-addSpectrum(d,by=1)
## get idea about number of age groups to include
agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
agetab.df<-as.data.frame(agetab)
ages<-1:15
## require at least 1 aged individual in each year
for(a in ages){
  if(any(agetab.df$Freq[agetab.df$Age==a]<1))
    d<-fixAgeGroup(d,age=a,fun=ifelse(a==min(ages),"min","mean"))
}
d<-subset(d,Age>=min(ages))

d.ysplit <- split(d, d$Year)
ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
              varCof=FALSE,maxK=50,mc.cores=mc.cores)
years <- levels(d[[1]]$Year)
n_years <- length(years)
ca <- d[[1]]
hh <- d[[2]]
hl <- d[[3]]
hl.ag <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = hl,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haul.id,data = hl,FUN=mean)
hl.ag <- hl.ag %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.ag$cpue <- hl.ag$HLNoAtLngt/hl.ag$HaulDur

plaiceQ4 <- data.frame(LngtCm=NA,Year=NA,Quarter=NA,Species="O",
                       HLNoAtLngt=NA,HaulDur=NA,cpue=NA,Age=NA,alk=NA)
ALK_Plaice4 <- list()
for(i in 1:n_years){
  dat.i = attr(ALK[[i]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  alk.i <- array(NA,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  for(j in 1:n_hauls){
    row=j
    alk.i[,,j] = NageByHaul(row, ALK[[i]], returnALK = TRUE)
  }
  ma <- apply(alk.i,c(1,2),FUN=mean)
  ALK_Plaice4[[years[i]]] <- apply(alk.i,c(1,2),FUN=mean)
  
  a <- melt(ALK_Plaice4[[years[i]]])
  names(a) <- c("LngtCm","Age","alk")
  a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])
  plaice.i <- hl.ag %>% filter(Year==years[i])
  
  p.i <- left_join(plaice.i,a)
  plaiceQ4 <- rbind(plaiceQ4,p.i)
  print(years[i])
}
prior_years <-1991:1993
alk1 <- ALK_Plaice4[["1994"]]
alk2 <- ALK_Plaice4[["1995"]]
alk3 <- ALK_Plaice4[["1996"]]
alk <- (alk1+alk2+alk3)/3
matplot(cm.breaks[-length(cm.breaks)], alk, lwd = 2, type = 'l')
a <- melt(alk)
names(a) <- c("LngtCm","Age","alk")
a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])

dada <- hl_markus %>% filter(Year %in% prior_years & species=="plaice" & Quarter==4)
names(dada)[c(3,5)] <- c("LngtCm","Species")
dada$Year <- factor(dada$Year)
hl.prior <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = dada,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haulID,data =  dada,FUN=mean)
hl.prior <- hl.prior %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.prior$cpue <- hl.prior$HLNoAtLngt/hl.prior$HaulDur
for(i in 1:length(prior_years)){
  plaice.i <- hl.prior %>% filter(Year==prior_years[i])
  
  p.i <- left_join(plaice.i,a)
  plaiceQ4 <- rbind(plaiceQ4,p.i)
  print(prior_years[i])
}
plaiceQ4 <- plaiceQ4[-1,]
#saveRDS(ALK_Cod1, file="ALK_Cod1.RData")
rm(list=setdiff(ls(),c('dAll','codQ1','codQ4','herQ1','herQ4','flounderQ1','flounderQ4'
                       ,'plaiceQ1','plaiceQ4','dabQ1','dabQ4','WB_rects','hl_markus')))
#####

# Dab ALK QI
#####
mc.cores<-1; library(parallel)
d<-subset(dAll, Species=="Limanda limanda",Quarter==1,StatRec %in% WB_rects)
years <- levels(d[[1]]$Year)
d<-subset(dAll, Species=="Limanda limanda",Quarter==1,StatRec %in% WB_rects,
          Year %in% years)
d<-addSpectrum(d,by=1)
## get idea about number of age groups to include
agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
agetab.df<-as.data.frame(agetab)
ages<-1:12
## require at least 1 aged individual in each year
for(a in ages){
  if(any(agetab.df$Freq[agetab.df$Age==a]<1))
    d<-fixAgeGroup(d,age=a,fun=ifelse(a==min(ages),"min","mean"))
}
d<-subset(d,Age>=min(ages))

d.ysplit <- split(d, d$Year)
ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
              varCof=FALSE,maxK=50,mc.cores=mc.cores)
years <- levels(d[[1]]$Year)
n_years <- length(years)
ca <- d[[1]]
hh <- d[[2]]
hl <- d[[3]]
hl.ag <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = hl,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haul.id,data = hl,FUN=mean)
hl.ag <- hl.ag %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.ag$cpue <- hl.ag$HLNoAtLngt/hl.ag$HaulDur

dabQ1 <- data.frame(LngtCm=NA,Year=NA,Quarter=NA,Species="O",
                       HLNoAtLngt=NA,HaulDur=NA,cpue=NA,Age=NA,alk=NA)
ALK_Dab1 <- list()
for(i in 1:n_years){
  dat.i = attr(ALK[[i]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  alk.i <- array(NA,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  for(j in 1:n_hauls){
    row=j
    alk.i[,,j] = NageByHaul(row, ALK[[i]], returnALK = TRUE)
  }
  ma <- apply(alk.i,c(1,2),FUN=mean)
  ALK_Dab1[[years[i]]] <- apply(alk.i,c(1,2),FUN=mean)
  
  a <- melt(ALK_Dab1[[years[i]]])
  names(a) <- c("LngtCm","Age","alk")
  a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])
  dab.i <- hl.ag %>% filter(Year==years[i])
  
  d.i <- left_join(dab.i,a)
  dabQ1 <- rbind(dabQ1,d.i)
  print(years[i])
}
prior_years <-1991:2007
alk1 <- ALK_Dab1[["2008"]]
alk2 <- ALK_Dab1[["2009"]]
alk3 <- ALK_Dab1[["2010"]]
alk <- (alk1+alk2+alk3)/3
matplot(cm.breaks[-length(cm.breaks)], alk, lwd = 2, type = 'l')
a <- melt(alk)
names(a) <- c("LngtCm","Age","alk")
a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])

dada <- hl_markus %>% filter(Year %in% prior_years & species=="dab" & Quarter==1)
names(dada)[c(3,5)] <- c("LngtCm","Species")
dada$Year <- factor(dada$Year)
hl.prior <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = dada,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haulID,data =  dada,FUN=mean)
hl.prior <- hl.prior %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.prior$cpue <- hl.prior$HLNoAtLngt/hl.prior$HaulDur
for(i in 1:length(prior_years)){
  dab.i <- hl.prior %>% filter(Year==prior_years[i])
  
  d.i <- left_join(dab.i,a)
  dabQ1 <- rbind(dabQ1,d.i)
  print(prior_years[i])
}
dabQ1 <- dabQ1[-1,]
#saveRDS(ALK_Cod1, file="ALK_Cod1.RData")
rm(list=setdiff(ls(),c('dAll','codQ1','codQ4','herQ1','herQ4','flounderQ1','flounderQ4'
                       ,'plaiceQ1','plaiceQ4','dabQ1','dabQ4','WB_rects','hl_markus')))

#####

# Dab ALK Q4
#####
mc.cores<-1; library(parallel)
d<-subset(dAll, Species=="Limanda limanda",Quarter==4,StatRec %in% WB_rects)
years <- levels(d[[1]]$Year)
d<-subset(dAll, Species=="Limanda limanda",Quarter==4,StatRec %in% WB_rects,
          Year %in% years)
d<-addSpectrum(d,by=1)
## get idea about number of age groups to include
agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
agetab.df<-as.data.frame(agetab)
ages<-1:12
## require at least 1 aged individual in each year
for(a in ages){
  if(any(agetab.df$Freq[agetab.df$Age==a]<1))
    d<-fixAgeGroup(d,age=a,fun=ifelse(a==min(ages),"min","mean"))
}
d<-subset(d,Age>=min(ages))

d.ysplit <- split(d, d$Year)
ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
              varCof=FALSE,maxK=50,mc.cores=mc.cores)
years <- levels(d[[1]]$Year)
n_years <- length(years)
ca <- d[[1]]
hh <- d[[2]]
hl <- d[[3]]
hl.ag <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = hl,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haul.id,data = hl,FUN=mean)
hl.ag <- hl.ag %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.ag$cpue <- hl.ag$HLNoAtLngt/hl.ag$HaulDur

dabQ4 <- data.frame(LngtCm=NA,Year=NA,Quarter=NA,Species="O",
                    HLNoAtLngt=NA,HaulDur=NA,cpue=NA,Age=NA,alk=NA)
ALK_Dab4 <- list()
for(i in 1:n_years){
  dat.i = attr(ALK[[i]], "data")
  cm.breaks = attr(dat.i, "cm.breaks")
  n_hauls <- length(unique(dat.i$haul.id))
  alk.i <- array(NA,dim=c(length(cm.breaks)-1,length(ages),n_hauls))
  for(j in 1:n_hauls){
    row=j
    alk.i[,,j] = NageByHaul(row, ALK[[i]], returnALK = TRUE)
  }
  ma <- apply(alk.i,c(1,2),FUN=mean)
  ALK_Dab4[[years[i]]] <- apply(alk.i,c(1,2),FUN=mean)
  
  a <- melt(ALK_Dab4[[years[i]]])
  names(a) <- c("LngtCm","Age","alk")
  a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])
  dab.i <- hl.ag %>% filter(Year==years[i])
  
  d.i <- left_join(dab.i,a)
  dabQ4 <- rbind(dabQ4,d.i)
  print(years[i])
}
prior_years <-1991:2007
alk1 <- ALK_Dab4[["2008"]]
alk2 <- ALK_Dab4[["2009"]]
alk3 <- ALK_Dab4[["2010"]]
alk <- (alk1+alk2+alk3)/3
matplot(cm.breaks[-length(cm.breaks)], alk, lwd = 2, type = 'l')
a <- melt(alk)
names(a) <- c("LngtCm","Age","alk")
a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])

dada <- hl_markus %>% filter(Year %in% prior_years & species=="dab" & Quarter==4)
names(dada)[c(3,5)] <- c("LngtCm","Species")
dada$Year <- factor(dada$Year)
hl.prior <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = dada,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haulID,data =  dada,FUN=mean)
hl.prior <- hl.prior %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.prior$cpue <- hl.prior$HLNoAtLngt/hl.prior$HaulDur
for(i in 1:length(prior_years)){
  dab.i <- hl.prior %>% filter(Year==prior_years[i])
  d.i <- left_join(dab.i,a)
  dabQ4 <- rbind(dabQ4,d.i)
  print(prior_years[i])
}
dabQ4 <- dabQ4[-1,]
#saveRDS(ALK_Cod1, file="ALK_Cod1.RData")
rm(list=setdiff(ls(),c('dAll','codQ1','codQ4','herQ1','herQ4','flounderQ1','flounderQ4'
                       ,'plaiceQ1','plaiceQ4','dabQ1','dabQ4','WB_rects','hl_markus')))

#####

hl_N.CB <- rbind(codQ1,codQ4,herQ1,herQ4,flounderQ1,flounderQ4,
                 plaiceQ1,plaiceQ4,dabQ1,dabQ4)

hl_N.CB$N <- hl_N.CB$cpue*hl_N.CB$alk*60 # CPUE - ind. pr. hour
write.table(hl_N.CB,"hl_N.CB.csv",row.names = FALSE,sep=';')
plotALKfit(ALK[[1]],row=5)
getAnywhere(plotALKfit)

