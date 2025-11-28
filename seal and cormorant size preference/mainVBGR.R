library(RTMB)
library(tidyverse)
########################
# Fit growth functions
########################

# cod growth functions ALK
#####
d <- rbind(readRDS("N_Cod1.RData"),readRDS("N_Cod4.RData"))
ca_cod <- d %>% filter(N!=0)
ca_cod$cohort <- ca_cod$Year-ca_cod$Age
ca_cod$t <- ca_cod$Age+ca_cod$qday



names(ca_cod)[names(ca_cod)=="LngtCm"] <- "Length"

df <- ca_cod %>% select(Length,t,cohort,N,Year,Species,Age)

tplot <- ca_cod %>%
  group_by(t,Age,Quarter,Year) %>%
  summarise(
    w_mean = sum(Length * N) / sum(N)
  )

dag <- aggregate(w_mean~Age+Quarter,data=tplot,FUN=mean)
plot(tplot$t,tplot$w_mean,xlim=c(0,20),ylim=c(0,120))
axis(1,
     at = seq(0, 20, by = 1),
     tck = 1, lty = 2, col = "gray")

plot(dag$Age+dag$Quarter/4-1/8,dag$w_mean,xlim=c(0,20),ylim=c(0,120))
axis(1,
     at = seq(0, 20, by = 1),
     tck = 1, lty = 2, col = "gray")

table(ca_cod$cohort)
cohortX <- 1988:2022

dat <- as.list(df %>% filter(cohort %in% cohortX & t>0.5))

rm(list=setdiff(ls(),c('dat','ca_cod')))
source('fitVBGR.R')

sdr <- fitVBGR(dat,weights=c(1,2)) # weights are for mean and sd, respectively
sdr

#####

#cod growth functions CA
#####
# working directory to coefficient-file
library(sf)
library(DATRAS)
wd_coef <- "C:/Users/mavast/Documents/GitHub/SOS data/"

#dAll <- readRDS(paste(wd_coef,"DatrasExchange.R",sep=""))
dAll <- readRDS("datras.R")
WB <- st_read(paste(wd_coef,"shapefiles/ICES_areas.shp",sep=""))
#plot(WB %>% filter(ICES_SUB %in% 22:24))
WB <- as(WB,"Spatial")

dAll <- addSpatialData(dAll,WB)

# cod
d<-subset(dAll, Species=="Gadus morhua",ICES_SUB %in% 22:24)
d<-addSpectrum(d,by=1)
options(warn = -1)
ca <- d[['CA']]
hh <- d[['HH']]
hl <- d[['HL']]

hl <- aggregate(HLNoAtLngt~haul.id+LngtCm,data=hl,FUN=sum)

ca$lengthID <- paste(ca$haul.id,ca$LngtCm)
hl$lengthID <- paste(hl$haul.id,hl$LngtCm)
ca_cod <- left_join(ca,hl %>% select(lengthID,HLNoAtLngt))
ca_cod$cohort <- as.numeric(as.character(ca_cod$Year))-ca_cod$Age
#ca_cod <- ca_cod %>% filter(Age>0 & !is.na(HLNoAtLngt))
ca_cod <- ca_cod %>% filter(!is.na(HLNoAtLngt))
hh$Jday <- yday(make_date(year=hh$Year,month=hh$Month,day=hh$Day))/yday(make_date(year=hh$Year,month=12,day=31))
ca_cod <- left_join(ca_cod,hh %>% select(haul.id,Jday))
ca_cod <- ca_cod %>% filter(!(is.na(Jday) | is.na(Age)))
ca_cod <- ca_cod %>% filter(!(Age<1 & Quarter==1))
#hist(ca_cod$LngtCm[ca_cod$Age==0])
ca_cod <- ca_cod %>% filter(!(Age<1 & Country=="DE" & Year %in% c(2015,2018)))
#hist(ca_cod$LngtCm[ca_cod$Age==0])

cohorts <- 1988:2022#sort(unique(ca_cod$cohort))[-c(1:6)] # cohort with sufficient observations to fit
ca_cod <- ca_cod %>% filter(cohort %in% cohorts)
ca_cod$logLength <- log(ca_cod$LngtCm)
ca_cod$aq <- paste(ca_cod$Age,ca_cod$Quarter)
ca_cod$yq <- paste(ca_cod$Year,ca_cod$Quarter)
t <- aggregate(Jday~yq,data=ca_cod,FUN=mean)
ca_cod$qday <- t$Jday[match(ca_cod$yq,t$yq)]
ca_cod$Length <- ca_cod$LngtCm
ca_cod$t <- ca_cod$Age+ca_cod$qday

ca_cod$Species <- "cod"

df <- ca_cod %>% select(Length,t,cohort,HLNoAtLngt,haul.id,Year,Species,Age,Quarter)

ca_cod$n <- 1
df.n <- aggregate(n~haul.id+Length,data=ca_cod,FUN=sum)

df <- left_join(df,df.n,)
df$mul <- df$HLNoAtLngt/df$n

dat <- as.list(df)
rm(list=setdiff(ls(),c('dat')))

source('fitVBGR_ca.R')

sdr <- fitVBGR_ca(dat,weights=c(1,2)) # weights are for mean and sd, respectively
sdr
#####

# flounder growth functions
#####
d <- rbind(readRDS("N_flounder1.RData"),readRDS("N_flounder4.RData"))
ca_flounder <- d %>% filter(N!=0)
ca_flounder$cohort <- ca_flounder$Year-ca_flounder$Age
ca_flounder$t <- ca_flounder$Age+ca_flounder$qday

names(ca_flounder)[names(ca_flounder)=="LngtCm"] <- "Length"

df <- ca_flounder %>% select(Length,t,cohort,N,Year,Species,Quarter,Age)
table(df$cohort)
cohortX <- 1984:2022

tplot <- df %>% filter(cohort %in% cohortX) %>%
  group_by(t,Age,Quarter,Year) %>%
  summarise(
    w_mean = sum(Length * N) / sum(N)
  )

dag <- aggregate(w_mean~Age+Quarter,data=tplot,FUN=mean,ylim=c(0,60))
plot(tplot$t,tplot$w_mean,xlim=c(0,20))
axis(1,
     at = seq(0, 20, by = 1),
     tck = 1, lty = 2, col = "gray")

lines(c(10,10),c(-10,100),lwd=2,col="red")
lines(c(2,2),c(-10,100),lwd=2,col="red")
lines(c(12/2,12/2),c(-10,100),lwd=2,col="red")

plot(dag$Age+dag$Quarter/4-1/8,dag$w_mean,xlim=c(0,20),ylim=c(0,60))
axis(1,
     at = seq(0, 20, by = 1),
     tck = 1, lty = 2, col = "gray")
lines(c(10,10),c(-10,100),lwd=2,col="red")
lines(c(2,2),c(-10,100),lwd=2,col="red")
lines(c(12/2,12/2),c(-10,100),lwd=2,col="red")
dat <- as.list(df %>% filter(cohort %in% cohortX & t>1.5))
dat$cohortX <- 2015

rm(list=setdiff(ls(),c('dat','ca_flounder')))
source('fitVBGR.R')

sdr <- fitVBGR(dat,weights=c(2,2))
sdr
#####

#flounder growth functions CA
#####
# working directory to coefficient-file
library(sf)
library(DATRAS)
wd_coef <- "C:/Users/mavast/Documents/GitHub/SOS data/"

#dAll <- readRDS(paste(wd_coef,"DatrasExchange.R",sep=""))
dAll <- readRDS("datras.R")
WB <- st_read(paste(wd_coef,"shapefiles/ICES_areas.shp",sep=""))
#plot(WB %>% filter(ICES_SUB %in% 22:24))
WB <- as(WB,"Spatial")

dAll <- addSpatialData(dAll,WB)

# flounder
d<-subset(dAll, Species=="Platichthys flesus",ICES_SUB %in% 22:24)
d<-addSpectrum(d,by=1)
options(warn = -1)
ca <- d[['CA']]
hh <- d[['HH']]
hl <- d[['HL']]

hl <- aggregate(HLNoAtLngt~haul.id+LngtCm,data=hl,FUN=sum)

ca$lengthID <- paste(ca$haul.id,ca$LngtCm)
hl$lengthID <- paste(hl$haul.id,hl$LngtCm)
ca_flounder <- left_join(ca,hl %>% select(lengthID,HLNoAtLngt))
ca_flounder$cohort <- as.numeric(as.character(ca_flounder$Year))-ca_flounder$Age
#ca_flounder <- ca_flounder %>% filter(Age>0 & !is.na(HLNoAtLngt))
ca_flounder <- ca_flounder %>% filter(!is.na(HLNoAtLngt))
hh$Jday <- yday(make_date(year=hh$Year,month=hh$Month,day=hh$Day))/yday(make_date(year=hh$Year,month=12,day=31))
ca_flounder <- left_join(ca_flounder,hh %>% select(haul.id,Jday))
ca_flounder <- ca_flounder %>% filter(!(is.na(Jday) | is.na(Age)))
#plot(ca_flounder$Age+ca_flounder$Jday,ca_flounder$LngtCm)
ca_flounder <- ca_flounder %>% filter(Age<40)
ca_flounder <- ca_flounder %>% filter(!(Age<1 & Quarter==1))
table(ca_flounder$cohort)
ca_flounder$Species <- "flounder"

cohorts <- 1988:2022#sort(unique(ca_flounder$cohort))[-c(1:6)] # cohort with sufficient observations to fit
ca_flounder <- ca_flounder %>% filter(cohort %in% cohorts)
ca_flounder$logLength <- log(ca_flounder$LngtCm)
ca_flounder$aq <- paste(ca_flounder$Age,ca_flounder$Quarter)
ca_flounder$yq <- paste(ca_flounder$Year,ca_flounder$Quarter)
t <- aggregate(Jday~yq,data=ca_flounder,FUN=mean)
ca_flounder$qday <- t$Jday[match(ca_flounder$yq,t$yq)]
ca_flounder$Length <- ca_flounder$LngtCm
ca_flounder$t <- ca_flounder$Age+ca_flounder$qday

df <- ca_flounder %>% select(Length,t,cohort,HLNoAtLngt,haul.id,Year,Species,Age,Quarter)

ca_flounder$n <- 1
df.n <- aggregate(n~haul.id+Length,data=ca_flounder,FUN=sum)

df <- left_join(df,df.n,)
df$mul <- df$HLNoAtLngt/df$n

dat <- as.list(df)
rm(list=setdiff(ls(),c('dat')))

source('fitVBGR_ca.R')

sdr <- fitVBGR_ca(dat,weights=c(3,3)) # weights are for mean and sd, respectively
sdr
#####

# plaice growth functions
#####
d <- rbind(readRDS("N_plaice1.RData"),readRDS("N_plaice4.RData"))
ca_plaice <- d %>% filter(N!=0)
ca_plaice$cohort <- ca_plaice$Year-ca_plaice$Age
ca_plaice$t <- ca_plaice$Age+ca_plaice$qday

names(ca_plaice)[names(ca_plaice)=="LngtCm"] <- "Length"

df <- ca_plaice %>% select(Length,t,cohort,N,Year,Species,Quarter,Age)
table(df$cohort)
cohortX <- 1984:2022

tplot <- df %>% filter(cohort %in% cohortX) %>%
  group_by(t,Age,Quarter,Year) %>%
  summarise(
    w_mean = sum(Length * N) / sum(N)
  )

dag <- aggregate(w_mean~Age+Quarter,data=tplot,FUN=mean,ylim=c(0,60))
plot(tplot$t,tplot$w_mean,xlim=c(0,20))
axis(1,
     at = seq(0, 20, by = 1),
     tck = 1, lty = 2, col = "gray")

lines(c(10,10),c(-10,100),lwd=2,col="red")
lines(c(2,2),c(-10,100),lwd=2,col="red")
lines(c(12/2,12/2),c(-10,100),lwd=2,col="red")

plot(dag$Age+dag$Quarter/4-1/8,dag$w_mean,xlim=c(0,20),ylim=c(0,60))
axis(1,
     at = seq(0, 20, by = 1),
     tck = 1, lty = 2, col = "gray")
lines(c(10,10),c(-10,100),lwd=2,col="red")
lines(c(2,2),c(-10,100),lwd=2,col="red")
lines(c(12/2,12/2),c(-10,100),lwd=2,col="red")
dat <- as.list(df %>% filter(cohort %in% cohortX & t>1))

rm(list=setdiff(ls(),c('dat','ca_plaice')))
source('fitVBGR.R')

sdr <- fitVBGR(dat,weights=c(2,2))
sdr
#####

#plaice growth functions CA
#####
# working directory to coefficient-file
library(sf)
library(DATRAS)
wd_coef <- "C:/Users/mavast/Documents/GitHub/SOS data/"

#dAll <- readRDS(paste(wd_coef,"DatrasExchange.R",sep=""))
dAll <- readRDS("datras.R")
WB <- st_read(paste(wd_coef,"shapefiles/ICES_areas.shp",sep=""))
#plot(WB %>% filter(ICES_SUB %in% 22:24))
WB <- as(WB,"Spatial")

dAll <- addSpatialData(dAll,WB)

# plaice
d<-subset(dAll, Species=="Pleuronectes platessa",ICES_SUB %in% 22:24)
d<-addSpectrum(d,by=1)
options(warn = -1)
ca <- d[['CA']]
hh <- d[['HH']]
hl <- d[['HL']]

hl <- aggregate(HLNoAtLngt~haul.id+LngtCm,data=hl,FUN=sum)

ca$lengthID <- paste(ca$haul.id,ca$LngtCm)
hl$lengthID <- paste(hl$haul.id,hl$LngtCm)
ca_plaice <- left_join(ca,hl %>% select(lengthID,HLNoAtLngt))
ca_plaice$cohort <- as.numeric(as.character(ca_plaice$Year))-ca_plaice$Age
#ca_plaice <- ca_plaice %>% filter(Age>0 & !is.na(HLNoAtLngt))
ca_plaice <- ca_plaice %>% filter(!is.na(HLNoAtLngt))
hh$Jday <- yday(make_date(year=hh$Year,month=hh$Month,day=hh$Day))/yday(make_date(year=hh$Year,month=12,day=31))
ca_plaice <- left_join(ca_plaice,hh %>% select(haul.id,Jday))
ca_plaice <- ca_plaice %>% filter(!(is.na(Jday) | is.na(Age)))
#plot(ca_plaice$Age+ca_plaice$Jday,ca_plaice$LngtCm)
ca_plaice <- ca_plaice %>% filter(!(Age>10 & LngtCm<20))
ca_plaice <- ca_plaice %>% filter(!(Age<1 & Quarter==1))
table(ca_plaice$cohort)
ca_plaice$Species <- "plaice"

cohorts <- 1990:2023#sort(unique(ca_plaice$cohort))[-c(1:6)] # cohort with sufficient observations to fit
ca_plaice <- ca_plaice %>% filter(cohort %in% cohorts)
ca_plaice$logLength <- log(ca_plaice$LngtCm)
ca_plaice$aq <- paste(ca_plaice$Age,ca_plaice$Quarter)
ca_plaice$yq <- paste(ca_plaice$Year,ca_plaice$Quarter)
t <- aggregate(Jday~yq,data=ca_plaice,FUN=mean)
ca_plaice$qday <- t$Jday[match(ca_plaice$yq,t$yq)]
ca_plaice$Length <- ca_plaice$LngtCm
ca_plaice$t <- ca_plaice$Age+ca_plaice$qday
ca_plaice$Species <- "plaice"
df <- ca_plaice %>% select(Length,t,cohort,HLNoAtLngt,haul.id,Year,Species,Age,Quarter)

ca_plaice$n <- 1
df.n <- aggregate(n~haul.id+Length,data=ca_plaice,FUN=sum)

df <- left_join(df,df.n,)
df$mul <- df$HLNoAtLngt/df$n

dat <- as.list(df)
rm(list=setdiff(ls(),c('dat')))

source('fitVBGR_ca.R')

sdr <- fitVBGR_ca(dat,weights=c(3,3)) # weights are for mean and sd, respectively
sdr
#####

# dab growth functions
#####
d <- rbind(readRDS("N_dab1.RData"),readRDS("N_dab4.RData"))
ca_dab <- d %>% filter(N!=0)
ca_dab$cohort <- ca_dab$Year-ca_dab$Age
ca_dab$t <- ca_dab$Age+ca_dab$qday

names(ca_dab)[names(ca_dab)=="LngtCm"] <- "Length"

df <- ca_dab %>% select(Length,t,cohort,N,Year,Species,Quarter,Age)
table(df$cohort)
cohortX <- c(1989:1997,2000:2023)

tplot <- df %>% filter(cohort %in% cohortX) %>%
  group_by(t,Age,Quarter,Year) %>%
  summarise(
    w_mean = sum(Length * N) / sum(N)
  )

dag <- aggregate(w_mean~Age+Quarter,data=tplot,FUN=mean,ylim=c(0,60))
plot(tplot$t,tplot$w_mean,xlim=c(0,20))
axis(1,
     at = seq(0, 20, by = 1),
     tck = 1, lty = 2, col = "gray")

lines(c(10,10),c(-10,100),lwd=2,col="red")
lines(c(2,2),c(-10,100),lwd=2,col="red")
lines(c(12/2,12/2),c(-10,100),lwd=2,col="red")

plot(dag$Age+dag$Quarter/4-1/8,dag$w_mean,xlim=c(0,20),ylim=c(0,60))
axis(1,
     at = seq(0, 20, by = 1),
     tck = 1, lty = 2, col = "gray")
lines(c(10,10),c(-10,100),lwd=2,col="red")
lines(c(2,2),c(-10,100),lwd=2,col="red")
lines(c(12/2,12/2),c(-10,100),lwd=2,col="red")
dat <- as.list(df %>% filter(cohort %in% cohortX & t>1.5))

rm(list=setdiff(ls(),c('dat','ca_dab')))
source('fitVBGR.R')

sdr <- fitVBGR(dat,weights=c(2,2))
sdr
#####

#dab growth functions CA
#####
# working directory to coefficient-file
library(sf)
library(DATRAS)
wd_coef <- "C:/Users/mavast/Documents/GitHub/SOS data/"

#dAll <- readRDS(paste(wd_coef,"DatrasExchange.R",sep=""))
dAll <- readRDS("datras.R")
WB <- st_read(paste(wd_coef,"shapefiles/ICES_areas.shp",sep=""))
#plot(WB %>% filter(ICES_SUB %in% 22:24))
WB <- as(WB,"Spatial")

dAll <- addSpatialData(dAll,WB)

# dab
d<-subset(dAll, Species=="Limanda limanda",ICES_SUB %in% 22:24)
d<-addSpectrum(d,by=1)
options(warn = -1)
ca <- d[['CA']]
hh <- d[['HH']]
hl <- d[['HL']]

hl <- aggregate(HLNoAtLngt~haul.id+LngtCm,data=hl,FUN=sum)

ca$lengthID <- paste(ca$haul.id,ca$LngtCm)
hl$lengthID <- paste(hl$haul.id,hl$LngtCm)
ca_dab <- left_join(ca,hl %>% select(lengthID,HLNoAtLngt))
ca_dab$cohort <- as.numeric(as.character(ca_dab$Year))-ca_dab$Age
#ca_dab <- ca_dab %>% filter(Age>0 & !is.na(HLNoAtLngt))
ca_dab <- ca_dab %>% filter(!is.na(HLNoAtLngt))
hh$Jday <- yday(make_date(year=hh$Year,month=hh$Month,day=hh$Day))/yday(make_date(year=hh$Year,month=12,day=31))
ca_dab <- left_join(ca_dab,hh %>% select(haul.id,Jday))
ca_dab <- ca_dab %>% filter(!(is.na(Jday) | is.na(Age)))
#plot(ca_dab$Age+ca_dab$Jday,ca_dab$LngtCm)

ca_dab <- ca_dab %>% filter(!(Age<1 & Quarter==1))
table(ca_dab$cohort)
ca_dab$Species <- "dab"

cohorts <- c(1989:1997,2000:2023) # cohort with sufficient observations to fit
ca_dab <- ca_dab %>% filter(cohort %in% cohorts)
ca_dab$logLength <- log(ca_dab$LngtCm)
ca_dab$aq <- paste(ca_dab$Age,ca_dab$Quarter)
ca_dab$yq <- paste(ca_dab$Year,ca_dab$Quarter)
t <- aggregate(Jday~yq,data=ca_dab,FUN=mean)
ca_dab$qday <- t$Jday[match(ca_dab$yq,t$yq)]
ca_dab$Length <- ca_dab$LngtCm
ca_dab$t <- ca_dab$Age+ca_dab$qday
ca_dab$Speceis <- "dab"
df <- ca_dab %>% select(Length,t,cohort,HLNoAtLngt,haul.id,Year,Species,Age,Quarter)

ca_dab$n <- 1
df.n <- aggregate(n~haul.id+Length,data=ca_dab,FUN=sum)

df <- left_join(df,df.n,)
df$mul <- df$HLNoAtLngt/df$n

dat <- as.list(df)
rm(list=setdiff(ls(),c('dat')))

source('fitVBGR_ca.R')

sdr <- fitVBGR_ca(dat,weights=c(3,3)) # weights are for mean and sd, respectively
sdr
#####

########################
# Plot growth functions
########################

source("plotVBGR.R")
plotVBGR(dat,sdr)
plotVBGR(dat,sdr,version="together")
plotVBGR(dat,sdr,version="data")
plotVBGR(dat,sdr,version="coefficients")

