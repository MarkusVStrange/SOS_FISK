########################
# a function to get individuals consumed pr. age group, based on biomass consumption
# 
#' @param consumption = a dataframe containing information of predator'a consumption of prey fish
#' 
#' 
#' @return 
#' a dataframe containing consumption in ind. pr. age pr. year of prey for the different predators
#########################

getConsInd <- function(consumption){
  wd <- "C:/Users/mavast/Documents/GitHub/SOS data/"
  LW <- read.table(paste(wd,"length-weight.csv",sep=""),header=TRUE,sep=';')
  
  source(paste(wd,"getFunction.R",sep=""))
  
  q_days <- c(1/8,3/8,5/8,7/8)
  
  cons.ind <- data.frame(N = NA,quarter=NA,ages=NA,years=NA,prey=NA,predator=NA)
  
  # Cormorants
  #####
  if('cormorant' %in% unique(consumption$predator)){
    Age.frac_corm <- read.table(paste(wd,"CormorantfoodAges_sim.csv",sep=""),header=TRUE,sep=';')
    cons.B_corm <- consumption %>% filter(predator=="cormorant")
    
    if("cod" %in% consumption$prey){
      cons_cod <- cons.B_corm %>% filter(prey=="cod")
      cod_frac <- Age.frac_corm %>% filter(species=="cod")
      cod_frac <- cod_frac[1:10]
      yrs_constant <- 1985:1990
      frac_constant <- rbind(cod_frac[1:4,],cod_frac[1:4,],cod_frac[1:4,],cod_frac[1:4,],cod_frac[1:4,],cod_frac[1:4,])
      frac_constant$sampling <- paste(rep(yrs_constant,each=4),1:4)
      
      cod_frac <- rbind(frac_constant,cod_frac)
      
      
      cons_frac <- cons_cod$consumption_t*(cod_frac[3:10])
      
      cod.growth <- getFunction(species="cod",type="growth")
      t_hatch <- ((16+197)/2)/365 # Cod time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
      length_matrix <- matrix(c(cod.growth(0:7+q_days[1]-t_hatch),
                                cod.growth(0:7+q_days[2]-t_hatch),
                                cod.growth(0:7+q_days[3]-t_hatch),
                                cod.growth(0:7+q_days[4]-t_hatch)),nrow=4,byrow=TRUE)
      cod.coef <- LW %>% filter(species=="cod")
      coef_constant <- data.frame(year=1985:1990,a=rep(cod.coef$a[1],6),
                                  b=rep(cod.coef$b[1],6),species="cod")
      cod.coef <- rbind(coef_constant,cod.coef[-length(cod.coef$year),])
      
      fish_weights <- rep(cod.coef$a,each=4)*
        length_matrix[rep(1:4,length.out = nrow(cod_frac)),]^(
          rep(cod.coef$b,each=4))
      
      cod_ind_cons <- ((cons_frac*10^6)/fish_weights)/1000 # unit is 1000 ind.
      cod_ind_cons[is.na(cod_ind_cons)] <- 0
    
      suppressMessages({df_cod.corm <- melt(cod_ind_cons,value.name="N")})
      df_cod.corm$quarter <- 1:4
      df_cod.corm$ages <- rep(0:7,each=length(cod_frac$sampling))
      df_cod.corm$years <- rep(1985:2023,each=4) 
      df_cod.corm$prey <- "cod"
      df_cod.corm$predator <- "cormorant"
      
      df_cod.corm <- df_cod.corm %>% select(-variable)
      
      cons.ind <- rbind(cons.ind,df_cod.corm)
    }
    
  } 
  #####
  
  
  
  # Grey seals
  #####
  if('grey seal' %in% unique(Predators$species)){
    Age.frac_seal <- read.table(paste(wd,"SealfoodAges_sim.csv",sep=""),header=TRUE,sep=';')
    cons.B_seal <- consumption %>% filter(predator=="grey seal")
    
    if("cod" %in% consumption$prey){
      cons_cod <- cons.B_seal %>% filter(prey=="cod")
      cod_frac <- Age.frac_seal %>% filter(species=="cod")
      cod_frac <- cod_frac[1:10]
      yrs_constant <- 1985:1990
      frac_constant <- rbind(cod_frac[1:4,],cod_frac[1:4,],cod_frac[1:4,],cod_frac[1:4,],cod_frac[1:4,],cod_frac[1:4,])
      frac_constant$sampling <- paste(rep(yrs_constant,each=4),1:4)
      
      cod_frac <- cod_frac %>% filter(as.numeric(substr(cod_frac$sampling,1,4)) %in% cons_cod$year)
      cons_frac <- cons_cod$consumption_t*(cod_frac[3:10])
      
      cod.growth <- getFunction(species="cod",type="growth")
      t_hatch <- ((16+197)/2)/365 # Cod time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
      length_matrix <- matrix(c(cod.growth(0:7+q_days[1]-t_hatch),
                                cod.growth(0:7+q_days[2]-t_hatch),
                                cod.growth(0:7+q_days[3]-t_hatch),
                                cod.growth(0:7+q_days[4]-t_hatch)),nrow=4,byrow=TRUE)
      cod.coef <- LW %>% filter(species=="cod")
      cod.coef <- cod.coef %>% filter(year %in% cons_cod$year)
      
      fish_weights <- rep(cod.coef$a,each=4)*
        length_matrix[rep(1:4,length.out = nrow(cod_frac)),]^(
          rep(cod.coef$b,each=4))
      
      cod_ind_cons <- ((cons_frac*10^6)/fish_weights)/1000 # unit is 1000 ind.
      cod_ind_cons[is.na(cod_ind_cons)] <- 0
      
      suppressMessages({df_cod.seal <- melt(cod_ind_cons,value.name="N")})
      df_cod.seal$quarter <- 1:4
      df_cod.seal$ages <- rep(0:7,each=length(cod_frac$sampling))
      df_cod.seal$years <- rep(unique(cons_cod$year),each=4) 
      df_cod.seal$prey <- "cod"
      df_cod.seal$predator <- "grey seal"
      
      df_cod.seal <- df_cod.seal %>% select(-variable)
      
      cons.ind <- rbind(cons.ind,df_cod.seal)
    }
  }
  #####
  
  cons.ind <- cons.ind[-1,]
  cons.ind
}