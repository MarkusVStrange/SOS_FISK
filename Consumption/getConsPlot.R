########################
# a function to plot the consumption of predators and fishery on a fish species
# 
#' @param prey = one prey species to plot. E.g., 'cod' 
#' @param Consumption = a dataframe with consumption data
#' @param type = a string specifying your wanted plot: "total", "biomass", "proportions", "recent"
#' 
#' 
#' @return 
#' a plot of consumption of a specific species
#########################



  

getConsPlot <- function(prey,consumption,type){
  data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="") # data working directory
  
  #####
  if(prey=='cod'){
    #from sam
    #TSB <- data.frame(B=c(48480,39312,39560,36266,25835,19730,20639,26697,35936,47799,56040,67807,71018,64829,57342,49628,47490,
    #                      43846,40447,41435,46248,40313,36314,24138,18603,19222,18800,19677,17450,21476,21812,16169,12598,14102,
    #                      13898,9421,4816,3583,3037,3512,4689),years=1985:2025)
    
    #catch_SAM <- c(31466,23831,27224,25648,16732,11946,11950,15877,23919,33087,34031,45760,42416,
    #               36452,39590,32038,27795,24079,23779,20485,25849,25682,20872,14279,9056,10959,
    #               13573,11755,9713,12597,13371,9759,6886,8306,7454,5168,2887,2234,1281,1825)
    
    # observations
    #wd.dat <- "C:/Users/mavast/Desktop/Github/SOS-FISK/data/user3-WBCod2025/WBCod2025/data/"
    #dat.sam <- read.sam.data(wd.dat)
    #catch.obs.nr <- dat.sam$Catchobs[,,1] # catch observations, ind.
    #weca <- dat.sam$mtrx$weca[,,1] # weight-at-age in the catch
    #catch.obs <- colSums(catch.obs.nr*weca)
    
    # Plot WB cod
    yrs <- c(1985:2023)
    
    # get early consumptions
    cod.cons <- aggregate(consumption_t~year+predator,data=consumption %>% filter(prey=='cod' & year %in% yrs),FUN=sum)
    
    sms_nis <- readRDS(paste(data_wd,"sms_mortality.RDS",sep=""))
    # read in WB cod population data
    #M_age <- read.table("data/M_sms.csv",header=TRUE,sep=';')
    # read in WB cod catch data
    #catch <- read.table("data/C_sms.csv",header=TRUE,sep=';')
    catch <- sms_nis$Yield %>% filter(years %in% yrs)
    # put catches into the consumption df
    cod.cons <- rbind(cod.cons,data.frame(year=catch$years,predator=rep('fishery',length(catch$years)),
                                          consumption_t=catch$Yield))
    
    
    # total WB cod biomass
    #wb_cod <- aggregate(B~years,data = M_age,FUN=sum)
    #wb_cod <- TSB[-length(TSB$B),]
    wb_cod <- data.frame(year=catch$years,cod_B=sms_nis$TSB$TSB[sms_nis$TSB$years %in% yrs])
    # rename the year and cod biomass column
    #names(wb_cod)[names(wb_cod)=='years'] <- 'year'
    #names(wb_cod)[names(wb_cod)=='B'] <- 'cod_B'
    # put total cod biomass into the consumption df
    cod.cons <- cod.cons %>%
      left_join(wb_cod %>% select(year, cod_B), by = "year")
    # find the biomass specific proportion eaten
    cod.cons$proportion <- cod.cons$consumption_t/cod.cons$cod_B
    #rename the predators
    cod.cons$predator[cod.cons$predator=="cormorant"] <- "danish cormorants"
    cod.cons$predator[cod.cons$predator=="grey seal"] <- "grey seals"
    cod.cons$predator[cod.cons$predator=="fishery"] <- "fisheries"
    cod.cons$predator[cod.cons$predator=="cod biomass"] <- "cod stock"
    
    # set the order of the stacked bars
    cod.cons$predator <- factor(cod.cons$predator,levels = c("grey seals",
                                                             "danish cormorants",
                                                             "fisheries"))
    
    rel_cons <- data.frame(year=rep(yrs,3),predator=rep(unique(cod.cons$predator),each=length(yrs)),
                           rel_B = rep(0,length(yrs))*3)
    for(i in 1:length(rel_cons$year)){
      dada <- cod.cons %>% filter(year==rel_cons$year[i])
      cons_Pi <- dada$consumption_t[dada$predator==rel_cons$predator[i]]
      cons_yi <- sum(dada$consumption_t)
      if(rel_cons$predator[i] %in% dada$predator) rel_cons$rel_B[i] <- cons_Pi/cons_yi
    }
    # plot consumption of WB cod
    ptot <- ggplot(data=cod.cons %>% filter(year>1990)) +geom_bar(aes(x = year, y = proportion*120000,fill=predator),
                                stat = "identity",position="stack",alpha=0.25,width=0.7,show.legend = FALSE)+   
    geom_line(aes(x = year, y = cod_B,color='cod stock'),linewidth=2)+
    geom_line(aes(x = year, y = consumption_t,color=predator),linewidth=2)+
    ggtitle("Cod in the western Baltic")+
    scale_fill_manual(values = c("danish cormorants" = "darkblue",
                                 "grey seals"="darkorange",
                                 "fisheries"="darkred"),name="")+
    scale_color_manual(values = c("cod stock" = "black","danish cormorants" = "darkblue",
                                  "grey seals"="darkorange","fisheries"="darkred"),name="")+
    xlab("year")+
    scale_y_continuous(name = 'biomass [tonnes]',limits = c(0,120000),expand=c(0,0),
                       sec.axis = sec_axis(~.*1/120000, name = "proportion removed")) +
      scale_x_continuous(expand = c(0, 0)) +
    theme_bw()+ theme(axis.line = element_line(color='black'),
                      plot.background = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank(),
                      text=element_text(color="black", size=18,family="serif"),
                      legend.position = c(0.8,0.85),
                      legend.title = element_blank())
    
    # plot consumption of WB cod in proportions
    pprop <- ggplot(data=cod.cons %>% filter(year>1985)) +geom_bar(aes(x = year, y = proportion*100,fill=predator),
                                                                   stat = "identity",position="stack",alpha=0.35,width=0.7,show.legend = FALSE)+   
      ggtitle("Cod in the western Baltic")+
      scale_fill_manual(values = c("danish cormorants" = "darkblue",
                                   "grey seals"="darkorange",
                                   "fisheries"="darkred"),name="")+
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(limits=c(0,200),expand = c(0, 0))+
      xlab("year")+ylab("stock biomass removed [%]")+
      theme_bw()+ theme(axis.line = element_line(color='black'),
                        plot.background = element_blank(),
                        text=element_text(color="black", size=18,family="serif"),
                        legend.position = c(0.8,0.85),
                        legend.title = element_blank())
    
    # plot consumption of WB cod 
    pbiom <- ggplot(data=cod.cons %>% filter(year>1985)) +
      geom_line(aes(x = year, y = cod_B/1000,color='cod stock'),linewidth=2)+
      geom_line(aes(x = year, y = consumption_t/1000,color=predator),linewidth=2)+
      ggtitle("Cod in the western Baltic")+
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(limits = c(0,120),expand = c(0, 0))+
      scale_color_manual(values = c("cod stock" = "black","danish cormorants" = "darkblue",
                                    "grey seals"="darkorange","fisheries"="darkred"),name="")+
      xlab("year")+ylab("cod biomass [1000 tonnes]")+
      theme_bw()+ theme(axis.line = element_line(color='black'),
                        plot.background = element_blank(),
                        text=element_text(color="black", size=18,family="serif"),
                        legend.position = c(0.8,0.85),
                        legend.title = element_blank())
    
    # plot consumption of WB cod 
    pbiom2015 <- ggplot(data=cod.cons %>% filter(year>2010)) +
      geom_line(aes(x = year, y = cod_B/1000,color='cod stock'),linewidth=2)+
      geom_line(aes(x = year, y = consumption_t/1000,color=predator),linewidth=2)+
      ggtitle("Cod in the western Baltic")+
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(limits = c(0,35),expand = c(0, 0))+
      scale_color_manual(values = c("cod stock" = "black","danish cormorants" = "darkblue",
                                    "grey seals"="darkorange","fisheries"="darkred"),name="")+
      xlab("year")+ylab("cod biomass [1000 tonnes]")+
      theme_bw()+ theme(axis.line = element_line(color='black'),
                        plot.background = element_blank(),
                        text=element_text(color="black", size=18,family="serif"),
                        legend.position = c(0.8,0.85),
                        legend.title = element_blank())
    
    df_tsb <- data.frame(year = yrs,tsb=sms_nis$TSB$TSB[sms_nis$TSB$years %in% yrs])
    
    # plot consumption of WB cod in proportions
    prelprop <- ggplot(data=rel_cons) +geom_bar(aes(x = year, y = rel_B*100,fill=predator),
                                                                   stat = "identity",position="stack",alpha=0.35,width=0.7)+   
      ggtitle("Cod in the western Baltic")+geom_line(data=df_tsb,aes(x=year,y=tsb/max(tsb)*100),linewidth=1.5)+
      scale_fill_manual(values = c("danish cormorants" = "darkblue",
                                   "grey seals"="darkorange",
                                   "fisheries"="darkred"),name="")+
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(limits=c(0,100.01),expand = c(0, 0))+
      xlab("year")+ylab("relative biomass removed [%]")+
      theme_bw()+ theme(axis.line = element_line(color='black'),
                        plot.background = element_blank(),
                        text=element_text(color="black", size=18,family="serif"),
                        legend.position = c(0.4,0.4),
                        legend.title = element_blank())
    
    ggsave("WB cod consumption.png",prelprop,width=9,height=6)
    
    # Create dummy data for manual legend
    legend_df <- data.frame(
      x = c(2000, 2000, 2000, 2000),
      y = c(-1, -1, -1, -1),  # or small values like 1000, 2000, etc.
      group = c("Cod stock", "Danish cormorants", "Grey seals", "Fisheries")
    )
    
    ptot <- ggplot(data=cod.cons %>% filter(year>1990)) +geom_bar(aes(x = year, y = proportion*90000,fill=predator),
                                                                  stat = "identity",position="stack",alpha=0.40,width=0.7,show.legend = FALSE)+   
      geom_line(aes(x = year, y = cod_B),color='black',linewidth=2)+
      ggtitle("Cod in the western Baltic")+
      scale_fill_manual(values = c("danish cormorants" = "darkblue",
                                   "grey seals"="darkorange",
                                   "fisheries"="darkred"),name="")+
      xlab("year")+
      scale_y_continuous(name = 'Total stock biomass [tonnes]',limits = c(0,90000),expand=c(0,0),
                         sec.axis = sec_axis(~.*1/90000, name = "proportion removed")) +
      scale_x_continuous(expand = c(0, 0))+ 
      theme_bw()
    
    
    ptot <- ptot+ geom_line(data = legend_df, aes(x = x, y = y, color = group),linewidth=2,show.legend = FALSE) +
      scale_color_manual(values = c(
        "Grey seals" = "darkorange",
        "Danish cormorants" = "darkblue",
        "Fisheries" = "darkred",
        "Cod stock" = "black"
      )) +
      labs(color = "Predators") +
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            text=element_text(color="black", size=14,family="serif"),
            legend.position = c(0.27,0.15),
            legend.title = element_blank())
    
  }
  if(type=="proportions") p <- pprop
  if(type=="biomass") p <- pbiom
  if(type=="recent") p <- pbiom2015
  if(type=="total") p <- ptot
  p
}