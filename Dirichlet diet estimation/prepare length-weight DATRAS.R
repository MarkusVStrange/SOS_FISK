
data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="") # data working directory

# Read data from DATRAS
#####
hldata <- read.table(paste(data_wd,"hldata.csv",sep=""),sep=',',header=TRUE)
cadata <- read.table(paste(data_wd,"cadata.csv",sep=""),sep=',',header=TRUE)
hhdata <- read.table(paste(data_wd,"hhdata.csv",sep=""),sep=',',header=TRUE)
hl.cpue <- read.table(paste(data_wd,"hl_cpue.csv",sep=""),sep=',',header=TRUE)
ca.cpue <- read.table(paste(data_wd,"ca_cpue.csv",sep=""),sep=',',header=TRUE)
#####

#Prepare data
#####

aphia.cod <- 126436
aphia.dab <- 127139
aphia.plaice <- 127143
aphia.flounder <- 127141
aphia.herring <- 126417
specs <- data.frame(Valid_Aphia = c(aphia.flounder,aphia.plaice,aphia.dab,aphia.cod,aphia.herring),
                    species = c('flounder','plaice','dab','cod','herring'))

hhdata$haulID <- paste(hhdata$Year,hhdata$Quarter,hhdata$Ship,hhdata$Gear,hhdata$HaulNo)
hh  <- hhdata %>% filter(haulID %in% hl.cpue$haulID)
table(table(hh$haulID))
doubles <- as.data.frame(table(hh$haulID))
table(doubles$Freq)
hh <- hh %>% filter(haulID %in% doubles$Var1[doubles$Freq==1])
rm(list = c('doubles'))


hldata$haulID <- paste(hldata$Year,hldata$Quarter,hldata$Ship,hldata$Gear,hldata$HaulNo)
#cadata$haulID <- paste(cadata$Year,cadata$Quarter,cadata$Ship,cadata$HaulNo)
hl  <- hldata %>% filter(Valid_Aphia %in% specs$Valid_Aphia & haulID %in% hh$haulID)
hl <- hl %>% left_join(specs, by='Valid_Aphia')
hl <- hl %>% left_join(hhdata %>% select(HaulDur,haulID),by='haulID')


# remove NA lengths
hl <- hl %>% filter(!(is.na(LngtClass)))
# Fix that some lengths are in mm and some in cm
hl$LngtClass[hl$LngtCode==0] <- round(hl$LngtClass[hl$LngtCode==0]/10)
hl$LngtClass[hl$LngtCode=='.'] <- round(hl$LngtClass[hl$LngtCode=='.']/10)

#hl$CPUE <- hl$HLNoAtLngt/(hl$HaulDur/60) # CPUE in fish pr. hour

hl <- hl %>% dplyr::select(Year,Quarter,LngtClass,HLNoAtLngt,species,haulID)
hl <- hl %>% filter(HLNoAtLngt>0)
#hl <- hl %>% select(Year,Quarter,LngtClass,CPUE,species,haulID)
#hl <- hl %>% filter(CPUE>0)

cadata$haulID <- paste(cadata$Year,cadata$Quarter,cadata$Ship,cadata$Gear,cadata$HaulNo)
cadata <- cadata %>% filter(IndWgt>0)

ca  <- cadata %>% filter(Valid_Aphia %in% specs$Valid_Aphia & haulID %in% hh$haulID)
#ca  <- cadata %>% filter(Valid_Aphia %in% hl.cpue$AphiaID  & haulID %in% hl.cpue$haulID)
ca <- ca %>% left_join(specs, by='Valid_Aphia')

# Fix that some lengths are in mm and some in cm
ca$LngtClass[ca$LngtCode==0] <- round(ca$LngtClass[ca$LngtCode==0]/10)
ca$LngtClass[ca$LngtCode=='.'] <- round(ca$LngtClass[ca$LngtCode=='.']/10)


ca <- ca %>% dplyr::select(CANoAtLngt,Quarter,Year,LngtClass,species,Age,haulID,Country,IndWgt)
ca <- as.data.frame(na.omit(ca))

#hl$HLNoAtLngt <- round(hl$HLNoAtLngt)# some numbers have decimals - ask Casper Berg why

# Make some sound cuts in the age-length file
hhWB <- hhdata %>% filter(haulID %in% ca$haulID)

#hhWB <- hhdata %>% filter(StatRec %in% WB_rectangles)
hhWB$jday <- yday(paste(hhWB$Year,'-',hhWB$Month,'-',hhWB$Day,sep=""))/365
survey.day <- hhWB %>% select(haulID, jday)
#survey.day <- aggregate(jday~Year+Quarter,data=hhWB,FUN=mean)
#survey.day$jday <- survey.day$jday/365
#quarter.day <- aggregate(jday~Quarter,data=survey.day,FUN=mean)
ca <- ca %>% left_join(survey.day,by=c('haulID'))



ca <- ca %>% filter(Age<30)
ca$LngtClass[ca$LngtClass>120] <- ca$LngtClass[ca$LngtClass>120]/10 # wrong unit in haulID "2000 4 26D4 50" 
ca <- ca[-which(ca$LngtClass>30 & ca$Age==0),]
ca <- ca[-which(ca$LngtClass>80 & ca$species!="cod"),]

dat <- left_join(hl,ca %>% dplyr::select(haulID,LngtClass,species,IndWgt),by=c('haulID','LngtClass','species'))
dat <- dat %>% filter(!is.na(IndWgt) & !is.na(HLNoAtLngt))


lw.cod <- dat %>% filter(species=="cod")
lw.herring <- dat %>% filter(species=="herring")
lw.plaice <- dat %>% filter(species==c("plaice"))
lw.flounder <- dat %>% filter(species==c("flounder"))
lw.dab <- dat %>% filter(species==c("dab"))



exp_fit <- function(L,a,b){
  a*L^b
}
#plot(lw.cod$LngtClass,lw.cod$IndWgt)
#####

# cod
#############
lw.cod <- lw.cod %>% filter(Quarter %in% c(1,4))# only consider data from quarter 1 and 4
lw.cod$yq <- paste(lw.cod$Year,lw.cod$Quarter) # redefine year as a variable
years <- sort(unique(lw.cod$Year))

# calculate mean and sd length for each year and age
n_years <- length(years) # number of years

# define empty matrices
cod.a <- rep(0,n_years)
cod.b <- rep(0,n_years)

# normal
for (i in 1:n_years){
  dada <- lw.cod %>% filter(Year==years[i])
  fit <- nls(IndWgt~exp_fit(LngtClass,a,b),data=dada,start=c(a=0.2,b=2.3),
             weights = HLNoAtLngt)
  cod.a[i] <- coef(fit)[1]
  cod.b[i] <- coef(fit)[2]
  
}
par(mfrow=c(1,2))
plot(dada$LngtClass,dada$IndWgt,xlab="length [cm]",ylab="weight [g]")
lines(1:100,exp_fit(1:100,cod.a[i],cod.b[i]),col="red",lwd=2)

colors <- colorRampPalette(c("red", "blue"))(n_years)

plot(1:100,cod.a[1]*(1:100)^cod.b[1],type = 'l',lwd=2,col=colors[1]
     ,xlab="length [cm]",ylab="weight [g]",main="cod",xlim = c(0,80),ylim=c(0,8000))
for(i in 2:n_years){
  lines(1:100,cod.a[i]*(1:100)^cod.b[i],type = 'l',lwd=2,col=colors[i])
}
cod_LW <- data.frame(year=years,a=cod.a,b=cod.b,
                     species=rep("cod",n_years))
################

# herring
#############
lw.herring <- lw.herring %>% filter(Quarter %in% c(1,4))# only consider data from quarter 1 and 4
lw.herring$yq <- paste(lw.herring$Year,lw.herring$Quarter) # redefine year as a variable
years <- sort(unique(lw.herring$Year))

# calculate mean and sd length for each year and age
n_years <- length(years) # number of years

# define empty matrices
herring.a <- rep(0,n_years)
herring.b <- rep(0,n_years)
herring.RMSE <- rep(0,n_years)

# normal
for (i in 1:n_years){
  dada <- lw.herring %>% filter(Year==years[i])
  fit <- nls(IndWgt~exp_fit(LngtClass,a,b),data=dada,start=c(a=0.001,b=3.5),
             weights=HLNoAtLngt)
  herring.a[i] <- coef(fit)[1]
  herring.b[i] <- coef(fit)[2]
  herring.RMSE[i] <- mean(sqrt((predict(fit,newdata = dada)-dada$IndWgt)^2))
  
}
colors <- colorRampPalette(c("red", "blue"))(n_years)
plot(dada$LngtClass,dada$IndWgt)
lines(1:40,coef(fit)[1]*(1:40)^coef(fit)[2],col="red",lwd=2)

par(mfrow=c(1,2))
plot(1:40,herring.a[1]*(1:40)^herring.b[1],type = 'l',lwd=2,col=colors[1]
     ,xlab="length [cm]",ylab="weight [g]",main="Herring",xlim=c(0,35))
for(i in 2:n_years){
  lines(1:40,herring.a[i]*(1:40)^herring.b[i],type = 'l',lwd=2,col=colors[i])
}
herring_LW <- data.frame(a=herring.a,b=herring.b,species=rep("herring",n_years))
plot(lw.herring$LngtClass,lw.herring$IndWgt,xlab="length [cm]",ylab="weight [g]")
fit <- nls(IndWgt~exp_fit(LngtClass,a,b),data=lw.herring,start=c(a=0.001,b=3.5),
           weights=HLNoAtLngt)
lines(1:40,coef(fit)[1]*(1:40)^coef(fit)[2],col="red",lwd=2)

herring_LW <- data.frame(year="generic",a=coef(fit)[1],b=coef(fit)[2],
                         species="herring")


################

# flounder
#############
lw.flounder <- lw.flounder %>% filter(Quarter %in% c(1,4))# only consider data from quarter 1 and 4
lw.flounder$yq <- paste(lw.flounder$Year,lw.flounder$Quarter) # redefine year as a variable
years <- sort(unique(lw.flounder$Year))

# calculate mean and sd length for each year and age
n_years <- length(years) # number of years

# define empty matrices
flounder.a <- rep(0,n_years)
flounder.b <- rep(0,n_years)
flounder.RMSE <- rep(0,n_years)

# normal
for (i in 1:n_years){
  dada <- lw.flounder %>% filter(Year==years[i])
  fit <- nls(IndWgt~exp_fit(LngtClass,a,b),data=dada,start=c(a=0.01,b=3),
             weights=HLNoAtLngt)
  flounder.a[i] <- coef(fit)[1]
  flounder.b[i] <- coef(fit)[2]
  flounder.RMSE[i] <- mean(sqrt((predict(fit,newdata = dada)-dada$IndWgt)^2))
  
}
colors <- colorRampPalette(c("red", "blue"))(n_years)
par(mfrow=c(1,2))
plot(dada$LngtClass,dada$IndWgt,xlab="length [cm]",ylab="weight [g]")
lines(1:50,coef(fit)[1]*(1:50)^coef(fit)[2],col="red",lwd=2)

plot(1:50,flounder.a[1]*(1:50)^flounder.b[1],type = 'l',lwd=2,col=colors[1]
     ,xlab="length [cm]",ylab="weight [g]",main="Flounder",xlim=c(0,45))
for(i in 2:n_years){
  lines(1:50,flounder.a[i]*(1:50)^flounder.b[i],type = 'l',lwd=2,col=colors[i])
}
flounder_LW <- data.frame(year=years,a=flounder.a,b=flounder.b,
                          species=rep("flounder",n_years))
################

# plaice
#############
lw.plaice <- lw.plaice %>% filter(Quarter %in% c(1,4))# only consider data from quarter 1 and 4
lw.plaice$yq <- paste(lw.plaice$Year,lw.plaice$Quarter) # redefine year as a variable
years <- sort(unique(lw.plaice$Year))

# calculate mean and sd length for each year and age
n_years <- length(years) # number of years

# define empty matrices
plaice.a <- rep(0,n_years)
plaice.b <- rep(0,n_years)
plaice.RMSE <- rep(0,n_years)

# normal
for (i in 1:n_years){
  dada <- lw.plaice %>% filter(Year==years[i])
  fit <- nls(IndWgt~exp_fit(LngtClass,a,b),data=dada,start=c(a=0.005,b=3),
             weights=HLNoAtLngt)
  plaice.a[i] <- coef(fit)[1]
  plaice.b[i] <- coef(fit)[2]
  plaice.RMSE[i] <- mean(sqrt((predict(fit,newdata = dada)-dada$IndWgt)^2))
  
}
colors <- colorRampPalette(c("red", "blue"))(n_years)
par(mfrow=c(1,2))
plot(dada$LngtClass,dada$IndWgt,xlab="length [cm]",ylab="weight [g]")
lines(1:50,coef(fit)[1]*(1:50)^coef(fit)[2],col="red",lwd=2)

plot(1:50,plaice.a[1]*(1:50)^plaice.b[1],type = 'l',lwd=2,col=colors[1]
     ,xlab="length [cm]",ylab="weight [g]",main="Plaice",xlim=c(0,45))
for(i in 2:n_years){
  lines(1:50,plaice.a[i]*(1:50)^plaice.b[i],type = 'l',lwd=2,col=colors[i])
}
plaice_LW <- data.frame(year=years,a=plaice.a,b=plaice.b,
                        species=rep("plaice",n_years))
################

# dab
#############
lw.dab <- lw.dab %>% filter(Quarter %in% c(1,4))# only consider data from quarter 1 and 4
lw.dab$yq <- paste(lw.dab$Year,lw.dab$Quarter) # redefine year as a variable
years <- sort(unique(lw.dab$Year))

# calculate mean and sd length for each year and age
n_years <- length(years) # number of years

# define empty matrices
dab.a <- rep(0,n_years)
dab.b <- rep(0,n_years)
dab.RMSE <- rep(0,n_years)

# normal
for (i in 1:n_years){
  dada <- lw.dab %>% filter(Year==years[i])
  fit <- nls(IndWgt~exp_fit(LngtClass,a,b),data=dada,start=c(a=0.02,b=3),
             weights=HLNoAtLngt)
  dab.a[i] <- coef(fit)[1]
  dab.b[i] <- coef(fit)[2]
  dab.RMSE[i] <- mean(sqrt((predict(fit,newdata = dada)-dada$IndWgt)^2))
  
}
colors <- colorRampPalette(c("red", "blue"))(n_years)
par(mfrow=c(1,2))
plot(dada$LngtClass,dada$IndWgt,xlab="length [cm]",ylab="weight [g]")
lines(1:50,coef(fit)[1]*(1:50)^coef(fit)[2],col="red",lwd=2)

plot(1:40,dab.a[1]*(1:40)^dab.b[1],type = 'l',lwd=2,col=colors[1],main="Dab",
     xlab="length [cm]",ylab="weight [g]",xlim=c(0,35))
for(i in 2:n_years){
  lines(1:40,dab.a[i]*(1:40)^dab.b[i],type = 'l',lwd=2,col=colors[i])
}
dab_LW <- data.frame(year=years,a=dab.a,b=dab.b,
                     species=rep("dab",n_years))
################

df.LW <- rbind(cod_LW,herring_LW,flounder_LW,plaice_LW,dab_LW)
write.table(df.LW,paste(data_wd,"length-weight.csv",sep=""),sep=';',row.names=FALSE)
#rm(list=setdiff(ls(),c('df.LW')))







