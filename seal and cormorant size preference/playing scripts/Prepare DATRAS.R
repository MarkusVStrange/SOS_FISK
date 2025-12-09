library(tidyverse)

data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="")
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
WB_statRec <- c("37G0","37G1","37G2","37G3","37G4",
                "38F9","38G0","38G1","38G2","38G3","38G4",
                "39F9","39G0","39G1","39G2","39G3","39G4",
                "40F9","40G0","40G1","40G2","41G0","41G2")


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
hh <- hh %>% filter(StatRec %in% WB_statRec)

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

hl <- hl %>% dplyr::select(Year,Quarter,LngtClass,HLNoAtLngt,species,haulID,HaulDur)
hl <- hl %>% filter(HLNoAtLngt>0)
#hl <- hl %>% select(Year,Quarter,LngtClass,CPUE,species,haulID)
#hl <- hl %>% filter(CPUE>0)
hl <- hl %>% left_join(hh %>% select(haulID,HaulDur))



cadata$haulID <- paste(cadata$Year,cadata$Quarter,cadata$Ship,cadata$Gear,cadata$HaulNo)
#cadata <- left_join(cadata,hh_unique)

ca  <- cadata %>% filter(Valid_Aphia %in% specs$Valid_Aphia & haulID %in% hh$haulID)
#ca  <- cadata %>% filter(Valid_Aphia %in% hl.cpue$AphiaID  & haulID %in% hl.cpue$haulID)
ca <- ca %>% left_join(specs, by='Valid_Aphia')

# Fix that some lengths are in mm and some in cm
ca$LngtClass[ca$LngtCode==0] <- round(ca$LngtClass[ca$LngtCode==0]/10)
ca$LngtClass[ca$LngtCode=='.'] <- round(ca$LngtClass[ca$LngtCode=='.']/10)

ca <- ca %>% dplyr::select(CANoAtLngt,Quarter,Year,LngtClass,species,Age,haulID,Country)
ca <- as.data.frame(na.omit(ca))

#hl$HLNoAtLngt <- round(hl$HLNoAtLngt)# some numbers have decimals - ask Casper Berg why

# Make some sound cuts in the age-length file
hhWB <- hhdata %>% filter(haulID %in% ca$haulID)

hhWB$jday <- yday(paste(hhWB$Year,'-',hhWB$Month,'-',hhWB$Day,sep=""))
survey.day <- aggregate(jday~Year+Quarter,data=hhWB,FUN=mean)

ca <- ca %>% left_join(survey.day)

ca <- ca %>% filter(Age<30)
ca$LngtClass[ca$LngtClass>120] <- ca$LngtClass[ca$LngtClass>120]/10 # wrong unit in haulID "2000 4 26D4 50" 
ca <- ca[-which(ca$LngtClass>30 & ca$Age==0),]
ca <- ca[-which(ca$LngtClass>80 & ca$species!="cod"),]
ca$LngtClass <- round(ca$LngtClass)
ca.cod <- ca %>% filter(species=="cod")
ca.herring <- ca %>% filter(species=="herring")
ca.plaice <- ca %>% filter(species==c("plaice"))
ca.flounder <- ca %>% filter(species==c("flounder"))
ca.dab <- ca %>% filter(species==c("dab"))


plot(ca.cod$Age+ca.cod$jday/365,ca.cod$LngtClass)
#plot(ca.herring$Age+ca.herring$jday,ca.herring$LngtClass)
#plot(ca.plaice$Age+ca.plaice$jday,ca.plaice$LngtClass)
#plot(ca.flounder$Age+ca.flounder$jday,ca.flounder$LngtClass)
#plot(ca.dab$Age+ca.dab$jday,ca.dab$LngtClass,xlim=c(0,12))
hl <- hl %>% filter(LngtClass<200)
hl.cod <- hl %>% filter(species=="cod")
hl.herring <- hl %>% filter(species=="herring")
hl.plaice <- hl %>% filter(species==c("plaice"))
hl.flounder <- hl %>% filter(species==c("flounder"))
hl.dab <- hl %>% filter(species==c("dab"))

#hist(hl.cod$LngtClass)
#hist(hl.herring$LngtClass)
#hist(hl.plaice$LngtClass)
#hist(hl.flounder$LngtClass)
#hist(hl.dab$LngtClass)


#names(hl.cpue)[names(hl.cpue) %in% c("AphiaID","LngtClas")] <- c("Valid_Aphia","LngtClass")
#hl.cpue <- hl.cpue %>% left_join(specs, by='Valid_Aphia')
#hl.cpue$LngtClass <- round(hl.cpue$LngtClass/10)
#table(hl$HLNoAtLngt>200)
#hist(hl$HLNoAtLngt)
#####
#rm(list=setdiff(ls(),c('hl','hl.cpue','ca','ca.cpue')))
#write.table(hl,"hl_markus.csv",row.names = FALSE,sep=';')

# Combine CA an HL to get age-structure population  - ALK by quarter only. OBS: currently preferred
#####
ca.ag <- aggregate(CANoAtLngt~Quarter+LngtClass+species+Age,data = ca,FUN=sum)

sum(ca.ag$CANoAtLngt)
sum(ca$CANoAtLngt)
# Construct age-length-key (ALK) as probabilities of age given length over years and quarters
alk <- ca.ag %>%
  group_by(Quarter,LngtClass,species, Age) %>%
  summarise(n_age = sum(CANoAtLngt), .groups = "drop") %>%
  group_by(Quarter,species,LngtClass) %>%
  mutate(p_age = n_age / sum(n_age)) 

alk <- as.data.frame(alk)

alk <- alk %>% filter(paste(Quarter,LngtClass,species) %in% unique(paste(hl$Quarter,hl$LngtClass,hl$species)))
#cpue.hl <- aggregate(CPUE~Quarter+Year+LngtClass+species,data = hl,FUN=sum)
# Join ALK to HL using Haul keys and LngtClass
hl_age <- hl %>%
  left_join(alk %>% dplyr::select(-n_age),
            by = c("Quarter","LngtClass","species"),relationship = "many-to-many") %>%
  mutate(N_age = HLNoAtLngt * p_age)

#hl_N <- rbind(hl_age %>% filter(!(is.na(N_age))),hl_nas)
hl_N <- hl_age %>% filter(!(is.na(N_age)))

hhWB <- hhdata
hhWB$jday <- yday(paste(hhWB$Year,'-',hhWB$Month,'-',hhWB$Day,sep=""))
survey.day <- aggregate(jday~Year+Quarter,data=hhWB,FUN=mean)
survey.day$jday <- survey.day$jday/365
quarter.day <- aggregate(jday~Quarter,data=survey.day,FUN=mean)

hl_N <- hl_N %>% left_join(quarter.day,by=c('Quarter'))

hl_N$jday[is.na(hl_N$jday)] <- 365/4+365/8
#####

rm(list=setdiff(ls(),c('hl_N')))











# Combine CA an HL to get age-structure population - ALK by year and quarter
#####
#ca.ag <- aggregate(CANoAtLngt~Quarter+Year+LngtClass+species+Age+haulID,data = ca,FUN=sum)

#sum(ca.ag$CANoAtLngt)
#sum(ca$CANoAtLngt)
# Construct age-length-key (ALK) as probabilities of age given length over years and quarters
#alk <- ca.ag %>%
#  group_by(Quarter,Year,LngtClass,species, Age) %>%
#  summarise(n_age = sum(CANoAtLngt), .groups = "drop") %>%
#  group_by(Quarter,Year,species,LngtClass) %>%
#  mutate(p_age = n_age / sum(n_age)) 

#alk <- as.data.frame(alk)
#no.hl <- aggregate(HLNoAtLngt~haulID+Year+Quarter+LngtClass+species,data = hl,FUN=sum)
#alk <- alk %>% filter(paste(Year,Quarter,LngtClass,species) %in% unique(paste(no.hl$Year,no.hl$Quarter,no.hl$LngtClass,no.hl$species)))
#cpue.hl <- aggregate(CPUE~Quarter+Year+LngtClass+species,data = hl,FUN=sum)
# Join ALK to HL using Haul keys and LngtClass
#hl_age <- no.hl %>%
#  left_join(alk %>% dplyr::select(-n_age),
#            by = c("Year","Quarter","LngtClass","species"),relationship = "many-to-many") %>%
#  mutate(N_age = HLNoAtLngt * p_age)


#hl_age <- cpue.hl %>%
#  left_join(alk %>% dplyr::select(-n_age),
#            by = c("Quarter", "Year","LngtClass","species")) %>%
#  mutate(N_age = CPUE * p_age)

#ca.ag.mean <- aggregate(CANoAtLngt~LngtClass+species+Age,data = ca,FUN=sum)

# Construct ALK as probabilities of age given length as mean of all years
#mean_alk <- ca.ag.mean %>%
#  group_by(LngtClass,species, Age) %>%
#  summarise(n_age = sum(CANoAtLngt), .groups = "drop") %>%
#  group_by(species,LngtClass) %>%
#  mutate(p_age = n_age / sum(n_age)) 
# remove lengths with no age key
#mean_alk$al <- paste(mean_alk$LngtClass,mean_alk$species)
#hl_age$al <- paste(hl_age$LngtClass,hl_age$species)
#hl_age <- hl_age %>% filter(al %in% mean_alk$al)
#mean_alk <- mean_alk %>% dplyr::select(-al)
#hl_age <- hl_age %>% dplyr::select(-al)



#hl_nas <- hl_age %>% filter(is.na(Age)) %>%
#  select(-c("Age","p_age","N_age"))


#hl_est.age <- hl_nas %>%
#  left_join(mean_alk %>% dplyr::select(-n_age),
#            by = c("LngtClass","species"),relationship = "many-to-many") %>%
#  mutate(N_age = HLNoAtLngt * p_age)


#nas <- hl_age %>% filter(is.na(N_age))
#hl_nas <- hl_age[1,]
#count <- 1
#for(i in (1:length(hl_age$N_age))[is.na(hl_age$N_age)]){
#  alk_i <- mean_alk %>% filter(LngtClass==hl_age$LngtClass[i],species==hl_age$species[i])
#  n <- length(alk_i$LngtClass)
#  hl_i <- hl_age[i,]
#  hl_new <- data.frame(haulID=rep(hl_i$haulID,n),Quarter=rep(hl_i$Quarter,n),Year=rep(hl_i$Year,n),
#                       LngtClass=rep(hl_i$LngtClass,n),species=rep(hl_i$species,n),
#                       HLNoAtLngt=rep(hl_i$HLNoAtLngt,n),Age=alk_i$Age,p_age=alk_i$p_age,N_age=NA)
#  hl_new$N_age <- hl_new$HLNoAtLngt*alk_i$p_age
#  hl_nas <- rbind(hl_nas,hl_new)
#  print(paste(round(count/length(is.na(hl_age$N_age)==TRUE),1),"% done"))
#  count <- count+1
#}
#hl_nas <- hl_nas[-1,]
#head(hl_nas)

#hl_N <- rbind(hl_age %>% filter(!(is.na(N_age))),hl_nas)
#hl_N <- rbind(hl_age %>% filter(!(is.na(N_age))),hl_est.age)

# a little test for conservation of fish
#ag1 <- aggregate(HLNoAtLngt~Quarter+Year+LngtClass+species,data=hl_age,FUN=mean)
#ag2 <- aggregate(HLNoAtLngt~Quarter+Year+LngtClass+species,data=hl_N,FUN=mean)
#sum(ag1$HLNoAtLngt)
#sum(ag2$HLNoAtLngt)

#hhdata$haulID <- paste(hhdata$Year,hhdata$Quarter,hhdata$Ship,hhdata$Gear,hhdata$HaulNo)
#hhWB <- hhdata %>% filter(haulID %in% ca$haulID)
#hhWB$jday <- yday(paste(hhWB$Year,'-',hhWB$Month,'-',hhWB$Day,sep=""))
#survey.day <- aggregate(jday~Year+Quarter,data=hhWB,FUN=mean)
#survey.day$jday <- survey.day$jday/365
#quarter.day <- aggregate(jday~Quarter,data=survey.day,FUN=mean)

#hl_N <- hl_N %>% left_join(quarter.day,by=c('Quarter'))

#hl_N$jday[is.na(hl_N$jday)] <- 365/4+365/8
#####

