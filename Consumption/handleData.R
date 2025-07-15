########################
# a function to handle  diet data.
# 

#' @return 
#' a dataframe containing: year, quarter, site, prey species, predator species, weighted and unweighted diets, and number of excrements
#########################
library(tidyverse)

handleData <- function(years,species){
  dietData <- data.frame(year = NA,quarter=NA,site=NA,prey=NA,unweighted.diet=NA,
                         weighted.diet=NA,n=NA,predator=NA) # n denote the number of pellets/scats
  
  data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="") # data working directory
  ###############
  # Cormorants
  ###############
  # 1992-1994
  #####
  if('cormorant' %in% species & TRUE %in% (1992:1994 %in% years)){
  
    d<- read.table(paste(data_wd,"corm_diet_prop.csv",sep=""),header=TRUE,
                   sep=";",as.is=TRUE)
    d$site <- str_sub(d$dmyl,-3)
    samp <- d
    names(samp)[names(samp)=="n"] <- "n_pellets"

    # define quarter
    samp$quarter <- rep('O',length(samp$year))
    samp$quarter[which(samp$month=="Jan" | samp$month=="Feb" | 
                         samp$month=="Mar")] <- 'Q1'
    samp$quarter[which(samp$month=="Apr" | samp$month=="May" | 
                         samp$month=="Jun")] <- 'Q2'
    samp$quarter[which(samp$month=="Jul" | samp$month=="Aug" | 
                         samp$month=="Sep")] <- 'Q3'
    samp$quarter[which(samp$month=="Oct" | samp$month=="Nov" | 
                         samp$month=="Dec")] <- 'Q4'
    
    samp$rel_biom <- rep(0,length(samp$dmyl))
    # total biomass eaten pr. sampling
    samp_tot <- aggregate(B~n_pellets+dmyl,data=samp,FUN = sum)
    # diet at each sampling
    for (i in 1:length(samp$dmyl)){
      idx <- samp$dmyl[i]
      samp$rel_biom[i] <- samp$B[i]/samp_tot$B[samp_tot$dmyl==idx]
    }
    
    # define quarter ID's
    samp$Q.ID <- paste(samp$year,samp$quarter,samp$site)
    
    # find the mean quarterly diet composition weighted by number of pellets in each sampling
    # multiply by number of pellets
    samp$weighted.diet <- samp$n_pellets*samp$rel_biom
    # the quarterly sums of relative biomass in the diet
    diet.92.94 <- aggregate(weighted.diet~year+quarter+site+species+Q.ID,data=samp,FUN=sum)
    # the number of pellets in each sampling
    n_pellets <- aggregate(n_pellets~Q.ID,data=aggregate(n_pellets~dmyl+Q.ID,data=samp,FUN=mean),FUN=sum)
    # put number of pellets in each sampling into the diet dataframe
    diet.92.94 <- diet.92.94 %>%
      left_join(n_pellets %>% select(Q.ID, n_pellets), by = "Q.ID")
    # divide the quarterly sums of relative biomass in the diet by the number of pellets in respective quarter
    diet.92.94$weighted.diet <- diet.92.94$weighted.diet/diet.92.94$n_pellets
    
    # find the unweighted quarterly diet. Mean of all samplings within a quarter
    # the sum of relative biomass in diet within a quarter
    unweighted <- aggregate(rel_biom~species+Q.ID,data=samp,FUN=sum)
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
    names(diet.92.94)[names(diet.92.94)=="species"] <- "prey"
    names(diet.92.94)[names(diet.92.94)=="n_pellets"] <- "n"
    diet.92.94$predator <- 'cormorant'
    # combine with the "mother" dataframe of all diets
    dietData <- rbind(dietData,diet.92.94)
    rm(list=setdiff(ls(),c('dietData','species','years','data_wd')))
    
  }  
  #####
  
  ###############
  # Grey seals
  ###############
  # 2015-2017
  #####
  if('grey seal' %in% species & TRUE %in% (2015:2017 %in% years)){
    # Load in the datafile
    D <- read.table(paste(data_wd,"gSeal_diet_prop.csv",sep=""),header=TRUE,
                    sep=";",as.is=TRUE)
    names(D)[names(D)=="n"] <- "n_scats"
    # remove data from Utklippan
    samp <- D %>% filter(!(Site=="Utklippan"))
    
    # define quarter
    samp$quarter <- rep('O',length(samp$Year))
    samp$quarter[which(samp$Month=="Jan" | samp$Month=="Feb" | 
                         samp$Month=="Mar")] <- 'Q1'
    samp$quarter[which(samp$Month=="Apr" | samp$Month=="May" | 
                         samp$Month=="Jun")] <- 'Q2'
    samp$quarter[which(samp$Month=="Jul" | samp$Month=="Aug" | 
                         samp$Month=="Sep")] <- 'Q3'
    samp$quarter[which(samp$Month=="Oct" | samp$Month=="Nov" | 
                         samp$Month=="Dec")] <- 'Q4'
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