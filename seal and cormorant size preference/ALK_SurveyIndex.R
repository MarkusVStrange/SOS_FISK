remotes::install_github("DTUAqua/DATRAS/DATRAS")
remotes::install_github("casperwberg/surveyIndex/surveyIndex")

library(surveyIndex)
library(DATRAS)

downloadExchange("BITS",1985)
hl_markus <- read.table("hl_markus.csv",header=TRUE,sep=';')
dAll<-readExchangeDir(".",strict=TRUE)
rm(list=setdiff(ls(),c('dAll','hl_markus')))

WB_rects <- c("37G2", "37G3", "37G4", "37G5", "38G2", "38G3", "38G4", 
                 "38G5", "39G4", "39G5", "39G6", "40G5", "40G6")

# Cod ALK QI
#####
mc.cores<-1; library(parallel)
d<-subset(dAll, Species=="Gadus morhua",Quarter==1,StatRec %in% WB_rects)
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

codQ1 <- data.frame(LngtCm=NA,Year=NA,Quarter=NA,Species="O",
                    HLNoAtLngt=NA,HaulDur=NA,cpue=NA,Age=NA,alk=NA)
ALK_Cod1 <- list()
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
  ALK_Cod1[[years[i]]] <- apply(alk.i,c(1,2),FUN=mean)
  
  a <- melt(ALK_Cod1[[years[i]]])
  names(a) <- c("LngtCm","Age","alk")
  a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])
  cod.i <- hl.ag %>% filter(Year==years[i])
  
  c.i <- left_join(cod.i,a)
  codQ1 <- rbind(codQ1,c.i)
  print(years[i])
}
codQ1 <- codQ1[-1,]
#saveRDS(ALK_Cod1, file="ALK_Cod1.RData")
rm(list=setdiff(ls(),c('dAll','codQ1','codQ4','herQ1','herQ4','flounderQ1','flounderQ4'
                       ,'plaiceQ1','plaiceQ4','dabQ1','dabQ4','WB_rects','hl_markus')))
#####

# Cod ALK Q4
#####
mc.cores<-1; library(parallel)
d<-subset(dAll, Species=="Gadus morhua",Quarter==4,StatRec %in% WB_rects)
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

codQ4 <- data.frame(LngtCm=NA,Year=NA,Quarter=NA,Species="O",
                    HLNoAtLngt=NA,HaulDur=NA,cpue=NA,Age=NA,alk=NA)
ALK_Cod4 <- list()
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
  ALK_Cod4[[years[i]]] <- apply(alk.i,c(1,2),FUN=mean)
  
  a <- melt(ALK_Cod4[[years[i]]])
  names(a) <- c("LngtCm","Age","alk")
  a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])
  cod.i <- hl.ag %>% filter(Year==years[i])
  
  c.i <- left_join(cod.i,a)
  codQ4 <- rbind(codQ4,c.i)
  print(years[i])
}
codQ4 <- codQ4[-1,]
#saveRDS(ALK_Cod1, file="ALK_Cod1.RData")
rm(list=setdiff(ls(),c('dAll','codQ1','codQ4','herQ1','herQ4','flounderQ1','flounderQ4'
                       ,'plaiceQ1','plaiceQ4','dabQ1','dabQ4','WB_rects','hl_markus')))
#####

# Herring ALK Q1
#####
mc.cores<-1; library(parallel)
d<-subset(dAll, Species=="Clupea harengus",Quarter==1,StatRec %in% WB_rects)
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

# Herring ALK Q4
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
d<-subset(dAll, Species=="Platichthys flesus",StatRec %in% WB_rects,
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
prior_years <-1991
alk1 <- ALK_Flounder1[["1992"]]
alk2 <- ALK_Flounder1[["1993"]]
alk3 <- ALK_Flounder1[["1994"]]
alk <- (alk1+alk2+alk3)/3
matplot(cm.breaks[-length(cm.breaks)], alk, lwd = 2, type = 'l')
a <- melt(alk)
names(a) <- c("LngtCm","Age","alk")
a$LngtCm <- round(cm.breaks[-1][round(a$LngtCm)])

dada <- hl_markus %>% filter(Year %in% prior_years & species=="flounder" & Quarter==1)
names(dada)[c(3,5)] <- c("LngtCm","Species")
dada$Year <- factor(dada$Year)
hl.prior <- aggregate(HLNoAtLngt~LngtCm+Year+Quarter+Species,data = dada,FUN=sum)
totHdur <- aggregate(HaulDur~Year+haulID,data =  dada,FUN=mean)
hl.prior <- hl.prior %>% left_join(aggregate(HaulDur~Year,data = totHdur,FUN=sum))
hl.prior$cpue <- hl.prior$HLNoAtLngt/hl.prior$HaulDur
for(i in 1:length(prior_years)){
  flounder.i <- hl.prior %>% filter(Year==prior_years[i])
  
  f.i <- left_join(flounder.i,a)
  flounderQ1 <- rbind(flounderQ1,f.i)
  print(prior_years[i])
}
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

