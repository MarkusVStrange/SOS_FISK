########################
# a function to handle  diet data.
# 

#' @return 
#' a dataframe containing: year, quarter, site, prey species, predator species, weighted and unweighted diets, and number of excrements
#########################
library(dplyr)

handleData <- function(years,species){
  dietData <- data.frame(year = NA,quarter=NA,site=NA,prey=NA,unweighted.diet=NA,
                         weighted.diet=NA,n=NA,predator=NA) # n denote the number of pellets/scats
  
  ###############
  # Cormorants
  ###############
  # 1992-1994
  #####
  if('cormorant' %in% species & TRUE %in% (1992:1994 %in% years)){
    d<- read.table("data/90s_cormorant_pellets.csv",header=TRUE,
                   sep=";",as.is=TRUE)
    info <- read.table("data/90s_sampling_information.csv",header=TRUE,
                       sep=";",as.is=TRUE)
    
    month <- c("January","February","March","April","May","June","July","August",
               "September","October","November","December")
    # Fix species names
    d$ART[d$ART=="SAN" | d$ART=="SOR"] <- 'goby'
    d$ART[d$ART=="TOR"] <- 'cod'
    d$ART[d$ART=="ISI" | d$ART=="SKB" | d$ART=="R\xd8S"] <- 'flatfish'
    d$ART[d$ART=="\xc5LA"] <- 'eelpout'
    d$ART[d$ART=="ULK"] <- 'sculpin'
    d$ART[d$ART=="TBI"] <- 'sandeel'
    d$ART[d$ART=="SIL"] <- 'herring'
    d$ART[-which(d$ART=="goby" | d$ART=="cod" | d$ART=="flatfish" | d$ART=="eelpout" |
                   d$ART=="sculpin" | d$ART=="sandeel" | d$ART=="herring")] <- 'other'
    
    
    d$day <- rep(0,length(d$GYLP))
    d$month <- rep(0,length(d$GYLP))
    d$year <- rep(0,length(d$GYLP))
    d$dmyl <- rep(0,length(d$GYLP))
    d$site <- rep(0,length(d$GYLP))
    
    d <- d[-which(d$FISKLGD==0 & d$FISKVGT==0),]
    d <- d[-which(is.na(d$FISKLGD) & is.na(d$FISKVGT)),]
    d$size <- d$FISKLGD
    d$size <- floor(d$size/10)*10
    
    for (i in 1:length(d$GYLP)){
      gylp <- d$GYLP[i]
      inf <- subset(info,info$GYLP==gylp)
      d$day[i] <- inf$DAG
      d$month[i] <- month[inf$MON]
      d$year[i] <- inf$YEAR+1900
      d$site[i] <- inf$LOKAL
      d$dmyl[i] <- paste(as.character(inf$DAG),as.character(inf$MON),as.character(inf$YEAR),inf$LOKAL)
      
    }
    
    
    d <- d[-which(is.na(d$size)),]
   
    samp <- aggregate(FISKVGT~dmyl+ART+day+month+year+site,data=d,FUN=sum)
    
    
    # get the number of pellets pr. sampling
    samp$n_pellets <- rep(0,length(samp$year))
    for (i in 1:length(samp$n_pellets)){
      dada <- d %>% filter(dmyl==samp$dmyl[i])
      samp$n_pellets[i] <- length(unique(dada$GYLP))
    }
    
    # define quarter
    samp$quarter <- rep('O',length(samp$year))
    samp$quarter[which(samp$month=="January" | samp$month=="February" | 
                         samp$month=="March")] <- 'Q1'
    samp$quarter[which(samp$month=="April" | samp$month=="May" | 
                         samp$month=="June")] <- 'Q2'
    samp$quarter[which(samp$month=="July" | samp$month=="August" | 
                         samp$month=="September")] <- 'Q3'
    samp$quarter[which(samp$month=="October" | samp$month=="November" | 
                         samp$month=="December")] <- 'Q4'
    
    samp$rel_biom <- rep(0,length(samp$dmyl))
    # total biomass eaten pr. sampling
    samp_tot <- aggregate(FISKVGT~n_pellets+dmyl,data=samp,FUN = sum)
    # diet at each sampling
    for (i in 1:length(samp$dmyl)){
      idx <- samp$dmyl[i]
      samp$rel_biom[i] <- samp$FISKVGT[i]/samp_tot$FISKVGT[samp_tot$dmyl==idx]
    }
    
    # define quarter ID's
    samp$Q.ID <- paste(samp$year,samp$quarter,samp$site)
    
    # find the mean quarterly diet composition weighted by number of pellets in each sampling
    # multiply by number of pellets
    samp$weighted.diet <- samp$n_pellets*samp$rel_biom
    # the quarterly sums of relative biomass in the diet
    diet.92.94 <- aggregate(weighted.diet~year+quarter+site+ART+Q.ID,data=samp,FUN=sum)
    # the number of pellets in each sampling
    n_pellets <- aggregate(n_pellets~Q.ID,data=aggregate(n_pellets~dmyl+Q.ID,data=samp,FUN=mean),FUN=sum)
    # put number of pellets in each sampling into the diet dataframe
    diet.92.94 <- diet.92.94 %>%
      left_join(n_pellets %>% select(Q.ID, n_pellets), by = "Q.ID")
    # divide the quarterly sums of relative biomass in the diet by the number of pellets in respective quarter
    diet.92.94$weighted.diet <- diet.92.94$weighted.diet/diet.92.94$n_pellets
    
    # find the unweighted quarterly diet. Mean of all samplings within a quarter
    # the sum of relative biomass in diet within a quarter
    unweighted <- aggregate(rel_biom~ART+Q.ID,data=samp,FUN=sum)
    # sampling ID's
    samplings <- unique(samp$Q.ID)
    # divide the sum of relative biomass in diet within a quarter with the number of samplings within said quarter
    for (i in 1:length(samplings)){
      dada <- samp %>% filter(Q.ID==samplings[i])
      dmyl <- length(unique(dada$dmyl))
      unweighted$rel_biom[unweighted$Q.ID==samplings[i]] <- unweighted$rel_biom[unweighted$Q.ID==samplings[i]]/dmyl
    }
    
    # put the unweighted diets into the diet dataframe. Trim the dataframe
    diet.92.94$unweighted.diet <- unweighted$rel_biom
    diet.92.94$Q.ID <- NULL
    diet.92.94 <- diet.92.94[c(1:4,7,5,6)]
    names(diet.92.94)[names(diet.92.94)=="ART"] <- "prey"
    names(diet.92.94)[names(diet.92.94)=="n_pellets"] <- "n"
    diet.92.94$predator <- rep('cormorant',length(diet.92.94$year))
    # combine with the "mother" dataframe of all diets
    dietData <- rbind(dietData,diet.92.94)
    rm(list=setdiff(ls(),c('dietData','species','years')))
    
  }  
  #####
  
  ###############
  # Grey seals
  ###############
  # 2015-2017
  #####
  if('grey seal' %in% species & TRUE %in% (2015:2017 %in% years)){
    # Load in the datafile
    D <- read.table("data/scat_otoliths.csv",header=TRUE,
                    sep=";",as.is=TRUE)
    # Convert important numbers tro numeric. Warning because some numbers are empty 
    suppressWarnings({D$FL.with.SCF <- as.numeric(D$FL.with.SCF)
                     D$FW.with.SCF <- as.numeric(D$FW.with.SCF)})
    
    # adjust the length unit of cod to match the other species
    D$FL.with.SCF[which(D$Species=="cod" | D$Species=="Cod" |
                          D$Species=="cod ")] <- D$FL.with.SCF[which(D$Species=="cod" | D$Species=="Cod" |
                                                                       D$Species=="cod ")]/10
    # convert length all units from cm to mm
    D$FL.with.SCF <- D$FL.with.SCF*10
    # implement 10 mm length classes. 
    D$size <- floor(D$FL.with.SCF/10)*10
    # remove NA and length = 0 observations
    D <- D[-which(is.na(D$size) | D$size==0),]

    # Clean up the species names
    D$Species[which(D$Species=="cod" | D$Species=="Cod" |
              D$Species=="cod ")] <- "cod"
    D$Species[which(D$Species=="flatfish" |
                      D$Species=="Flatfish")] <- "flatfish"
    D$Species[which(D$Species=="GOBY" | 
                      D$Species=="goby")] <- "goby"
    D$Species[which(D$Species=="herring" | 
                      D$Species=="Herring")] <- "herring"
    D$Species[which(D$Species=="sprat" | 
                      D$Species=="Sprat")] <- "sprat"
    D$Species[which(D$Species=="whiting" | 
                      D$Species=="Whiting")] <- "whiting"
    D$Species[which(D$Species=="four-bearded rockling" | 
                      D$Species=="Unidentified")] <- "other"
    # define a sampling ID 
    D$ymp <- paste(D$Year,D$Month,D$Site)
    
    # these are the species considered
    fish <- c("cod","herring","flatfish","goby",
              "sandeel","sprat","whiting","other")
    
    # NCF's from Lundström et al., 2007 
    # Implement number correction factors, i.e., a fraction of eaten otoliths are never found. This differ between species.
    # mean between plaice and flounder used for flatfish NCF - assuming flatfish are largely plaice, dab and flounder, and that NCF's are similar between these species.
    # Cod NCF used for whiting
    # NFC for "other" is a mean value of the named species
    # Species are cod, herring, flatfish, goby, sandeel, sprat,whiting, others
    NCF <- data.frame(species = fish, ncf = c(1.2,3.7,1.65,6.3,4.5,5.9,1.2,mean(c(1.2,3.7,1.65,6.3,4.5,5.9,1.2)))) 
    # aggregate all eaten fish by samplings based on sampling year, month, and location
    samp <- aggregate(FW.with.SCF~Year+Month+Site+Species+ymp,data=D,FUN = sum)
    # get the number of scats pr. sampling
    samp$n_scats <- rep(0,length(samp$Year))
    for (i in 1:length(samp$n_scats)){
      d <- D %>% filter(ymp==samp$ymp[i])
      samp$n_scats[i] <- length(unique(d$Scat))
    }
    
    # define a column for the NCF corrected biomass
    samp$corB <- rep(0,length(samp$ymp))
    # define a column for the relative diet pr. sampling
    samp$rel_biom <- rep(0,length(samp$ymp))
    
    # correction with NCF's
    for (i in 1:length(samp$ymp)){
      idx <- samp$ymp[i]
      sp <- samp$Species[i]
      samp$corB[i] <- samp$FW.with.SCF[i]*NCF$ncf[NCF$species==sp]
    }
    
    # total biomass eaten pr. sampling
    samp_tot <- aggregate(corB~Year+Month+Site+ymp,data=samp,FUN = sum)
    for (i in 1:length(samp$ymp)){
      idx <- samp$ymp[i]
      samp$rel_biom[i] <- samp$corB[i]/samp_tot$corB[samp_tot$ymp==idx]
      
    }
    # remove data from Utklippan
    samp <- samp %>% filter(!(Site=="Utklippan"))
    
    # define quarter
    samp$quarter <- rep('O',length(samp$Year))
    samp$quarter[which(samp$Month=="January" | samp$Month=="February" | 
                         samp$Month=="March")] <- 'Q1'
    samp$quarter[which(samp$Month=="April" | samp$Month=="May" | 
                         samp$Month=="June")] <- 'Q2'
    samp$quarter[which(samp$Month=="July" | samp$Month=="August" | 
                         samp$Month=="September")] <- 'Q3'
    samp$quarter[which(samp$Month=="October" | samp$Month=="November" | 
                         samp$Month=="December")] <- 'Q4'
    # define quarter ID's
    samp$Q.ID <- paste(samp$Year,samp$quarter,samp$Site)
    
    # find the mean quarterly diet composition weighted by number of scats in each sampling
    # multiply by number of scats
    samp$weighted.diet <- samp$n_scats*samp$rel_biom
    # the quarterly sums of relative biomass in the diet
    diet.15.17 <- aggregate(weighted.diet~Year+quarter+Site+Species+Q.ID,data=samp,FUN=sum)
    # the number of scats in each sampling
    n_scats <- aggregate(n_scats~Q.ID,data=aggregate(n_scats~Month+Q.ID,data=samp,FUN=mean),FUN=sum)
    # put number of scats in each sampling into the diet dataframe
    diet.15.17 <- diet.15.17 %>%
      left_join(n_scats %>% select(Q.ID, n_scats), by = "Q.ID")
    # divide the quarterly sums of relative biomass in the diet by the number of scats in respective quarter
    diet.15.17$weighted.diet <- diet.15.17$weighted.diet/diet.15.17$n_scats
    
    # find the unweighted quarterly diet. Mean of all samplings within a quarter
    # the sum of relative biomass in diet within a quarter
    unweighted <- aggregate(rel_biom~Species+Q.ID,data=samp,FUN=sum)
    # sampling ID's
    samplings <- unique(samp$Q.ID)
    # divide the sum of relative biomass in diet within a quarter with the number of samplings within said quarter
    for (i in 1:length(samplings)){
      d <- samp %>% filter(Q.ID==samplings[i])
      months <- length(unique(d$Month))
      unweighted$rel_biom[unweighted$Q.ID==samplings[i]] <- unweighted$rel_biom[unweighted$Q.ID==samplings[i]]/months
    }
    
    # put the unweighted diets into the diet dataframe. Trim the dataframe
    diet.15.17$unweighted.diet <- unweighted$rel_biom
    diet.15.17$Q.ID <- NULL
    diet.15.17 <- diet.15.17[c(1:4,7,5,6)]
    names(diet.15.17)[names(diet.15.17)=="Year"] <- "year"
    names(diet.15.17)[names(diet.15.17)=="Site"] <- "site"
    names(diet.15.17)[names(diet.15.17)=="Species"] <- "prey"
    names(diet.15.17)[names(diet.15.17)=="n_scats"] <- "n"
    diet.15.17$predator <- rep('grey seal',length(diet.15.17$year))
    # combine with the "mother" dataframe of all diets
    dietData <- rbind(dietData,diet.15.17)
    rm(list=setdiff(ls(),c('dietData','species','years')))
  }
  #####
  # remove the first row of NA's
  dietData <- dietData[-1,]
  dietData
}