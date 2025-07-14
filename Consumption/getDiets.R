########################
# a function to get the number of predators hunting in a region of interest
# 
#' @param Pradators = a dataframe containing information of numbers of species in years and quarters
#' @param method = a string containing either 'model' or 'data', which emphasize which method to use
#' 'model' gets time varying diet based on prey stocks, whereas 'data' gets a time constant diet based on the mean of available years
#' 
#' 
#' @return 
#' a dataframe containing diets of predators in specified years
#########################

getDiets <- function(Predators,method){
  diet <- data.frame(year = NA,quarter=NA,prey=NA,diet=NA,n=NA,predator=NA)
  
  source('handleData.R')
  d <- handleData(years = unique(Predators$year), species = unique(Predators$species))
  
  # Cormorants
  #####
  if(method=='data' & 'cormorant' %in% unique(Predators$species)){
    # take out cormorant data
    cormorants <- d %>% filter(predator=='cormorant')
    
    # find the weighted diet in each quarter by multiplying all observation by their n, and then dividing by the total n of the quarter
    # the multiplication of n
    cormorants$prelim <- cormorants$weighted.diet*cormorants$n
    # the summized diets in each quarter
    corm.diet <- aggregate(prelim~prey+quarter,data = cormorants,FUN=sum)
    # the number of pellets in each quarter
    Q.pellets <- aggregate(n~quarter,data=aggregate(n~quarter+site+year,data = cormorants,FUN=mean),FUN=sum)
    # put number of pellets in each sampling into the diet dataframe
    corm.diet <- corm.diet %>%
      left_join(Q.pellets %>% dplyr::select(quarter, n), by = "quarter")
    # divide by the total number of pellets in each quarter
    corm.diet$diet <- corm.diet$prelim/corm.diet$n
    # set the consumption of sculpin, herring and sandeel in Q4 to 0 manually
    corm.diet <- rbind(corm.diet,corm.diet[29,],corm.diet[28:29,])
    corm.diet$prey[c(29,31,32)] <- c('herring','sandeel','sculpin')
    corm.diet[c(29,31,32),c(3,5)] <- rep(0,6)
    
    # put it into a correct dataframe
    corm.diet.fin <- data.frame(year = rep(unique(Predators$year),each=length(corm.diet$quarter)),
                                quarter = rep(corm.diet$quarter,length(unique(Predators$year))),
                                prey = rep(corm.diet$prey,length(unique(Predators$year))),
                                diet = rep(corm.diet$diet,length(unique(Predators$year))),
                                predator = rep('cormorant',length(unique(Predators$year))*length(corm.diet$quarter)))
    # put the quarterly number of pellets into the dataframe
    corm.diet.fin <- corm.diet.fin %>%
      left_join(Q.pellets %>% dplyr::select(quarter, n), by = "quarter")
    # combine with the "mother" dataframe of all diets
    diet <- rbind(diet,corm.diet.fin)
    rm(list=setdiff(ls(),c('diet','Predators','method','d')))
  }
  
  
  if(method=='model' & 'cormorant' %in% unique(Predators$species)){
    
  }
  #####
  
  
  
  # Grey seals
  #####
  if(method=='data' & 'grey seal' %in% unique(Predators$species)){
    # take out grey seal data
    g.seals <- d %>% filter(predator=='grey seal')
    # find the weighted diet in each quarter by multiplying all observation by their n, and then dividing by the total n of the quarter
    # the multiplication of n
    g.seals$prelim <- g.seals$weighted.diet*g.seals$n
    # the summized diets in each quarter
    g.seal.diet <- aggregate(prelim~prey+quarter,data = g.seals,FUN=sum)
    # the number of scats in each quarter
    Q.scats <- aggregate(n~quarter,data=aggregate(n~quarter+site+year,data = g.seals,FUN=mean),FUN=sum)
    # put number of pellets in each sampling into the diet dataframe
    g.seal.diet <- g.seal.diet %>%
      left_join(Q.scats %>% dplyr::select(quarter, n), by = "quarter")
    #divide by the total number of pellets in each quarter
    g.seal.diet$diet <- g.seal.diet$prelim/g.seal.diet$n
    # set the consumption of "other species" in Q3 to 0 manually
    g.seal.diet <- rbind(g.seal.diet[1:21,],g.seal.diet[21:31,])
    g.seal.diet$prey[21] <- 'other'
    g.seal.diet[21,c(3,5)] <- c(0,0)
    
    # put it into a correct dataframe
    g.seal.diet.fin <- data.frame(year = rep(unique(Predators$year),each=length(g.seal.diet$quarter)),
                                quarter = rep(g.seal.diet$quarter,length(unique(Predators$year))),
                                prey = rep(g.seal.diet$prey,length(unique(Predators$year))),
                                diet = rep(g.seal.diet$diet,length(unique(Predators$year))),
                                predator = rep('grey seal',length(unique(Predators$year))*length(g.seal.diet$quarter)))
    # put the quarterly number of scats into the dataframe
    g.seal.diet.fin <- g.seal.diet.fin %>%
      left_join(Q.scats %>% dplyr::select(quarter, n), by = "quarter")
    # combine with the "mother" dataframe of all diets
    diet <- rbind(diet,g.seal.diet.fin)
    rm(list=setdiff(ls(),c('diet','Predators','method','d')))
  }
  
  
  if(method=='model' & 'grey seal' %in% unique(Predators$species)){
    
  }
  #####
  
  diet <- diet[-1,]
  diet
}