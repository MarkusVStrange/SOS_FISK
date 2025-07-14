########################
# a function to get biomass fish eaten by a predator, based on bioenergetics
# 
#' @param Diets = a dataframe containing information of predators and their diets
#' 
#' 
#' @return 
#' a dataframe containing consumption in tonnes/year of prey for the different predators
#########################

getConsumption <- function(g_eaten,Predators){
  Consumption <- data.frame(year = NA,quarter=NA,prey=NA,consumption_t=NA,predator=NA)
  # create a year and quarter variable for merging of dataframes
  g_eaten$yq <- paste(g_eaten$year,g_eaten$quarter)
  Predators$yq <- paste(Predators$year,Predators$quarter)
  # Cormorants
  #####
  if('cormorant' %in% unique(Predators$species)){
    # get cormorant data
    corm.eat <- g_eaten %>% filter(predator=='cormorant')
    cormorants <- Predators %>% filter(species=='cormorant')
    
    # merge number of birds into the eaten dataframe
    corm.eat <- corm.eat %>%
      left_join(cormorants %>% select(yq, N_preds), by = "yq")
    
    # calculate overall consumption (in tonnes)
    corm.eat$consumption_t <- corm.eat$N_preds*corm.eat$g_eaten/10^6*(365/4) # calculate quarterly consumption - convert from grams to tonnes. Extrapolate to the whole quarter, assume the are all equally long
    # trim the dataframe
    corm.eat <- corm.eat[,names(Consumption)]
    # put into mother dataframe
    Consumption <- rbind(Consumption,corm.eat)
    
  } 
  #####
  
  
  
  # Grey seals
  #####
  if('grey seal' %in% unique(Predators$species)){
    # get grey seal data
    g.seal.eat <- g_eaten %>% filter(predator=='grey seal')
    greyseals <- Predators %>% filter(species=='grey seal')
    # only take years when there is grey seal counts
    g.seal.eat <- g.seal.eat %>% filter(year %in% greyseals$year)
    
    
    # merge number of birds into the eaten dataframe
    g.seal.eat <- g.seal.eat %>%
      left_join(greyseals %>% select(yq, N_preds), by = "yq")
    
    # calculate overall consumption (in tonnes)
    g.seal.eat$consumption_t <- g.seal.eat$N_preds*g.seal.eat$g_eaten/10^6*(365/4) # calculate quarterly consumption - convert from grams to tonnes. Extrapolate to the whole quarter, assume the are all equally long
    # trim the dataframe
    g.seal.eat <- g.seal.eat[,names(Consumption)]
    # put into mother dataframe
    Consumption <- rbind(Consumption,g.seal.eat)
  }
  #####
  
  Consumption <- Consumption[-1,]
  Consumption
}