data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="") # data working directory
corm.diet <-  read.table(paste(data_wd,"corm_diet_22_24.csv",sep=""),header=TRUE,sep=';')
fishes <- c("Gadus morhua" ="cod","Platichthys flesus"="flatfish",
            "Pleuronectes platessa"="flatfish","Limanda limanda"="flatfish",
            "other"="other")

corm.diet$Species <- as.character(fishes[corm.diet$Species])


# Put in 0-values
Mal_flat <- corm.diet %>% filter(Colony_1=="Malurtholm" & Species=="cod")
Mal_flat[c(1,5)] <- c("flatfish",0)
Vor_cod <- (corm.diet %>% filter(Colony_1=="Vorsø" & Species=="flatfish" & Month==8))[1,]
Vor_cod[c(1,5)] <- c("cod",0)
corm.diet <- rbind(corm.diet,Mal_flat,Vor_cod)
corm.diet$weight <- as.numeric(corm.diet$weight)

d <- aggregate(weight~Species+Colony_1+Year+Month+total_fish_weight+n,data=corm.diet,FUN=sum)
rm(list=setdiff(ls(),c('d')))

library(vegan)


df <- data.frame(Site=d$Colony_1[d$Species=="cod"],
                 month=d$Month[d$Species=="cod"],
                 year=d$Year[d$Species=="cod"],
                 cod=d$weight[d$Species=="cod"],
                 flatfish=d$weight[d$Species=="flatfish"],
                 other=d$weight[d$Species=="other"])

df[,4:6] <- df[,4:6]/rowSums(df[,4:6])



df$Site <- as.numeric(as.factor(df$Site))

dat <- df[,4:6]

x <- metaMDS(dat,k=2)
x
s <- numeric(10)
for(i in 1:10){
  x <- metaMDS(df,k=i,silent=TRUE)
  s[i] <- (x$stress)
}
plot(x, display = c("sites", "species"))



plot(x)               # Quick ordination plot
ordiplot(x, type="text")   # Labels sites directly
orditorp(x, display="species", col="blue") # Add species labels nicely













nmds <- metaMDS(dat,k=2)
# Check stress
nmds$stress
# Lower stress (~<0.2 is okay; <0.1 is good for 2D) — interpret cautiously

# Extract site scores (points) and species scores (optional)
site_scores <- as.data.frame(scores(nmds, display = "sites"))
# Add a grouping variable from dune.env, e.g. Management
site_scores$colony <- df$Site
site_scores$year <- factor(df$year)

# Optional: species scores
species_scores <- as.data.frame(scores(nmds, display = "species"))
species_scores$Species <- rownames(species_scores)

# Plot NMDS ordination — sites colored by Management, species as text
ggplot(site_scores, aes(x = NMDS1, y = NMDS2, color = year)) +
  geom_point(size = 3, stroke = 0.5) +
  stat_ellipse(aes(group = year), linetype = "dashed", show.legend = FALSE) +
  geom_text(data = species_scores, aes(x = NMDS1, y = NMDS2, label = Species),
            inherit.aes = FALSE, size = 3, alpha = 0.6) +
  labs(title = "NMDS (metaMDS) of dune species",
       subtitle = paste0("Stress = ", round(nmds$stress, 3)),
       x = "NMDS1", y = "NMDS2") +
  theme_minimal(base_size = 14)
