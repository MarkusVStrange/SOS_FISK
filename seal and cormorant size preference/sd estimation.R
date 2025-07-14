source('prepare DATRAS.R')
##############
# cod
#############
cod <- hl_N %>% filter(species=='cod' & Age>1) # Look at cod age 2+
cod <- cod %>% filter(!(LngtClass<20 & Age>5)) # remove a 10 cm 8 year old cod. Must be a mistake
cod <- cod %>% filter(Quarter %in% c(1,4))# only consider data from quarter 1 and 4
cod$year <- as.numeric(substr(cod$haulID,1,4)) # redefine year as a variable
cod$yd <- cod$Age+cod$jday # define a julian day pseudo age
ages <- sort(unique(cod$yd))
years <- sort(unique(cod$year))

# Group by age only
byAge <- aggregate(N_age~yd+LngtClass,data=cod,FUN=sum) # group by age
nAge  <- aggregate(N_age~yd,data=cod,FUN=sum) # Count observations
X <- byAge %>% filter(yd==ages[1]) # choose the age to plot
# normal
m <-  sum(X$LngtClass*X$N_age)/sum(X$N_age) # weighted mean
s <- sqrt(sum(X$N_age*(X$LngtClass-m)^2)/sum(X$N_age)) # weighted sd
r <- round(rnorm(sum(X$N_age),mean=m,sd=s)) # simulated data
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r))) #make simulation a df

ggplot(X, aes(x = (LngtClass), y = N_age)) + # plot the observed and simulated length-frequencies together
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# log-normal
m <-  sum(log(X$LngtClass)*X$N_age)/sum(X$N_age) # weighted mean of log(observations)
s <- sqrt(sum(X$N_age*(log(X$LngtClass)-m)^2)/sum(X$N_age))# weighted sd of log(observations)
r <- round(exp(rnorm(sum(X$N_age),mean = m,sd = s)))-1 # simulated data
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r))) #make simulation a df

ggplot(X, aes(x = (LngtClass), y = N_age)) +# plot the observed and simulated length-frequencies together
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# Group by age and year
byYear <- aggregate(N_age~yd+LngtClass+year,data=cod,FUN=sum) # group by age and year
X <- byYear %>% filter(yd==ages[14] & year==2001) # select age and year to plot
# normal
m <-  sum(X$LngtClass*X$N_age)/sum(X$N_age) # can be inactive if a comparison to the only age-based simulation
s <- sqrt(sum(X$N_age*(X$LngtClass-m)^2)/sum(X$N_age)) # can be inactive if a comparison to the only age-based simulation
r <- round(rnorm(sum(X$N_age),mean=m,sd=s))
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r)))

ggplot(X, aes(x = (LngtClass), y = N_age)) +
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# log-normal
m <-  sum(log(X$LngtClass)*X$N_age)/sum(X$N_age) # can be inactive if a comparison to the only age-based simulation
s <- sqrt(sum(X$N_age*(log(X$LngtClass)-m)^2)/sum(X$N_age)) # can be inactive if a comparison to the only age-based simulation
r <- round(exp(rnorm(sum(X$N_age),mean = m,sd = s)))-1
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r)))

ggplot(X, aes(x = (LngtClass), y = N_age)) +
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()


# calculate mean and sd length for each year and age
n_years <- length(years) # number of years
n_ages <- length(ages) # number of ages
byYear <- aggregate(N_age~yd+LngtClass+year,data=cod,FUN=sum) # data by length, year, and age
# define empty matrices
mat_sd <- matrix(NA,ncol=n_years,nrow = n_ages)
mat_mean <- matrix(NA,ncol=n_years,nrow = n_ages)
mat_n <- matrix(NA,ncol=n_years,nrow = n_ages)
# normal
for (i in 1:n_years){
  for (j in 1:n_ages){
    dada <- byYear %>% filter(yd==ages[j] & year==years[i])
    obs <- dada$LngtClass
    mat_mean[j,i] <- sum(obs*dada$N_age)/sum(dada$N_age)
    mat_sd[j,i] <- sqrt(sum(dada$N_age*(obs-mat_mean[j,i])^2)/sum(dada$N_age))
    mat_n[j,i] <- sum(dada$N_age)
  }
}

mat_sd[mat_sd==0] <- NA # do not consider sd's of 0
sd_red <- mat_sd[-c(12,14:20),] # remove ages with low number of observations
mean_red <- mat_mean[-c(12,14:20),]
n_red <- mat_n[-c(12,14:20),]
# Convert SD to data frame with row and column indices for plot purposes
dsd <- as.data.frame(sd_red)
dsd$row <- 1:nrow(dsd)
# Pivot to long format
df_sd <- dsd %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_sd, aes(x = col+min(years)-1, y = row, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "SD [cm]") +
  theme_minimal()


# Convert Mean to data frame with row and column indices for plot purposes
dmean <- as.data.frame(mean_red)
dmean$row <- 1:nrow(dmean)
# Pivot to long format
df_mean <- dmean %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_mean, aes(x = col+min(years)-1, y = row, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "Mean [cm]") +
  theme_minimal()

# Convert n to data frame with row and column indices for plot purposes
dn <- as.data.frame(n_red)
dn$row <- 1:nrow(dn)
# Pivot to long format
df_n <- dn %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_mean, aes(x = col+min(years)-1, y = row, fill = log10(value))) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "log10(ind.)") +
  theme_minimal()

mean_sds <- rowMeans(sd_red,na.rm = TRUE)
mean_means <- rowMeans(mean_red,na.rm = TRUE)
n_sds <- apply(n_red,1,FUN=mean,na.rm = TRUE)

plot(mat_n[1,],mat_sd[1,])
plot(ages[-c(12,14:20)],mean_sds)
plot(mean_means,mean_sds)

df.cod <- data.frame(yd = ages[-c(12,14:20)],
                     mean=mean_means,sd=mean_sds,
                     species=rep("cod",length(mean_means)),
                     n=n_sds) # mean catch pr. quarter
# log-normal if necessary at some point
for (i in 1:n_years){
  for (j in 1:length(ages)){
    dada <- byYear %>% filter(yd==ages[j] & year==years[i])
    Lobs <- log(dada$LngtClass)
    mat_mean[j,i] <- sum(Lobs*dada$N_age)/sum(dada$N_age)
    mat_sd[j,i] <- sqrt(sum(dada$N_age*(Lobs-mat_mean[j,i])^2)/sum(dada$N_age))
    mat_n[j,i] <- sum(dada$N_age)
  }
}

rm(list=setdiff(ls(),c('hl_N','df.cod')))

################

##############
# herring
#############
herring <- hl_N %>% filter(species=='herring' & Age>1) # Look at cod age 2+
herring$year <- as.numeric(substr(herring$haulID,1,4)) # redefine year as a variable
herring$yd <- herring$Age+herring$jday # define a julian day pseudo age
ages <- sort(unique(herring$yd))
years <- sort(unique(herring$year))

# Group by age only
byAge <- aggregate(N_age~yd+LngtClass,data=herring,FUN=sum) # group by age
nAge  <- aggregate(N_age~yd,data=herring,FUN=sum) # Count observations
X <- byAge %>% filter(yd==ages[12]) # choose the age to plot
# normal
m <-  sum(X$LngtClass*X$N_age)/sum(X$N_age) # weighted mean
s <- sqrt(sum(X$N_age*(X$LngtClass-m)^2)/sum(X$N_age)) # weighted sd
r <- round(rnorm(sum(X$N_age),mean=m,sd=s)) # simulated data
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r))) #make simulation a df

ggplot(X, aes(x = (LngtClass), y = N_age)) + # plot the observed and simulated length-frequencies together
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# log-normal
m <-  sum(log(X$LngtClass)*X$N_age)/sum(X$N_age) # weighted mean of log(observations)
s <- sqrt(sum(X$N_age*(log(X$LngtClass)-m)^2)/sum(X$N_age))# weighted sd of log(observations)
r <- round(exp(rnorm(sum(X$N_age),mean = m,sd = s)))-1 # simulated data
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r))) #make simulation a df

ggplot(X, aes(x = (LngtClass), y = N_age)) +# plot the observed and simulated length-frequencies together
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# Group by age and year
byYear <- aggregate(N_age~yd+LngtClass+year,data=herring,FUN=sum) # group by age and year
X <- byYear %>% filter(yd==ages[1] & year==2001) # select age and year to plot
# normal
m <-  sum(X$LngtClass*X$N_age)/sum(X$N_age) # can be inactive if a comparison to the only age-based simulation
s <- sqrt(sum(X$N_age*(X$LngtClass-m)^2)/sum(X$N_age)) # can be inactive if a comparison to the only age-based simulation
r <- round(rnorm(sum(X$N_age),mean=m,sd=s))
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r)))

ggplot(X, aes(x = (LngtClass), y = N_age)) +
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# log-normal
m <-  sum(log(X$LngtClass)*X$N_age)/sum(X$N_age) # can be inactive if a comparison to the only age-based simulation
s <- sqrt(sum(X$N_age*(log(X$LngtClass)-m)^2)/sum(X$N_age)) # can be inactive if a comparison to the only age-based simulation
r <- round(exp(rnorm(sum(X$N_age),mean = m,sd = s)))-1
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r)))

ggplot(X, aes(x = (LngtClass), y = N_age)) +
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()


# calculate mean and sd length for each year and age
n_years <- length(years) # number of years
n_ages <- length(ages) # number of ages
byYear <- aggregate(N_age~yd+LngtClass+year,data=herring,FUN=sum) # data by length, year, and age
# define empty matrices
mat_sd <- matrix(NA,ncol=n_years,nrow = length(ages))
mat_mean <- matrix(NA,ncol=n_years,nrow = length(ages))
mat_n <- matrix(NA,ncol=n_years,nrow = length(ages))
# normal
for (i in 1:n_years){
  for (j in 1:n_ages){
    dada <- byYear %>% filter(yd==ages[j] & year==years[i])
    obs <- dada$LngtClass
    mat_mean[j,i] <- sum(obs*dada$N_age)/sum(dada$N_age)
    mat_sd[j,i] <- sqrt(sum(dada$N_age*(obs-mat_mean[j,i])^2)/sum(dada$N_age))
    mat_n[j,i] <- sum(dada$N_age)
  }
}

mat_sd[mat_sd==0] <- NA # do not consider sd's of 0
sd_red <- mat_sd[-c(18:20),] # remove ages with low number of observations
mean_red <- mat_mean[-c(18:20),]
n_red <- mat_n[-c(18:20),]
# Convert to data frame with row and column indices for plot purposes
dsd <- as.data.frame(sd_red)
dsd$row <- 1:nrow(dsd)
# Pivot to long format
df_sd <- dsd %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_sd, aes(x = col+min(years)-1, y = row, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "SD [cm]") +
  theme_minimal()


# Convert to data frame with row and column indices for plot purposes
dmean <- as.data.frame(mean_red)
dmean$row <- 1:nrow(dmean)
# Pivot to long format
df_mean <- dmean %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_mean, aes(x = col+min(years)-1, y = row, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "Mean [cm]") +
  theme_minimal()
# Convert n to data frame with row and column indices for plot purposes
dn <- as.data.frame(n_red)
dn$row <- 1:nrow(dn)
# Pivot to long format
df_n <- dn %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_mean, aes(x = col+min(years)-1, y = row, fill = log10(value))) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "log10(ind.)") +
  theme_minimal()

mean_sds <- rowMeans(sd_red,na.rm = TRUE)
mean_means <- rowMeans(mean_red,na.rm = TRUE)
n_sds <- apply(n_red,1,FUN=mean,na.rm = TRUE)

plot(mat_n[1,],mat_sd[1,])
plot(ages[-c(18:20)],mean_sds)
plot(mean_means,mean_sds)

df.herring <- data.frame(yd = ages[-c(18:20)],
                     mean=mean_means,sd=mean_sds,
                     species=rep("herring",length(mean_means)),
                     n=n_sds) # mean catch pr. quarter
# log-normal if necessary at some point
for (i in 1:n_years){
  for (j in 1:length(ages)){
    dada <- byYear %>% filter(yd==ages[j] & year==years[i])
    Lobs <- log(dada$LngtClass)
    mat_mean[j,i] <- sum(Lobs*dada$N_age)/sum(dada$N_age)
    mat_sd[j,i] <- sqrt(sum(dada$N_age*(Lobs-mat_mean[j,i])^2)/sum(dada$N_age))
    mat_n[j,i] <- sum(dada$N_age)
  }
}
rm(list=setdiff(ls(),c('hl_N','df.cod','df.herring')))

#################

##############
# flounder
#############
flounder <- hl_N %>% filter(species=='flounder' & Age>1) # Look at cod age 2+
flounder$year <- as.numeric(substr(flounder$haulID,1,4)) # redefine year as a variable
flounder$yd <- flounder$Age+flounder$jday # define a julian day pseudo age
ages <- sort(unique(flounder$yd))
years <- sort(unique(flounder$year))

# Group by age only
byAge <- aggregate(N_age~yd+LngtClass,data=flounder,FUN=sum) # group by age
nAge  <- aggregate(N_age~yd,data=flounder,FUN=sum) # Count observations
X <- byAge %>% filter(yd==ages[2]) # choose the age to plot
# normal
m <-  sum(X$LngtClass*X$N_age)/sum(X$N_age) # weighted mean
s <- sqrt(sum(X$N_age*(X$LngtClass-m)^2)/sum(X$N_age)) # weighted sd
r <- round(rnorm(sum(X$N_age),mean=m,sd=s)) # simulated data
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r))) #make simulation a df

ggplot(X, aes(x = (LngtClass), y = N_age)) + # plot the observed and simulated length-frequencies together
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# log-normal
m <-  sum(log(X$LngtClass)*X$N_age)/sum(X$N_age) # weighted mean of log(observations)
s <- sqrt(sum(X$N_age*(log(X$LngtClass)-m)^2)/sum(X$N_age))# weighted sd of log(observations)
r <- round(exp(rnorm(sum(X$N_age),mean = m,sd = s)))-1 # simulated data
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r))) #make simulation a df

ggplot(X, aes(x = (LngtClass), y = N_age)) +# plot the observed and simulated length-frequencies together
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# Group by age and year
byYear <- aggregate(N_age~yd+LngtClass+year,data=flounder,FUN=sum) # group by age and year
X <- byYear %>% filter(yd==ages[1] & year==1997) # select age and year to plot
# normal
m <-  sum(X$LngtClass*X$N_age)/sum(X$N_age) # can be inactive if a comparison to the only age-based simulation
s <- sqrt(sum(X$N_age*(X$LngtClass-m)^2)/sum(X$N_age)) # can be inactive if a comparison to the only age-based simulation
r <- round(rnorm(sum(X$N_age),mean=m,sd=s))
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r)))

ggplot(X, aes(x = (LngtClass), y = N_age)) +
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# log-normal
m <-  sum(log(X$LngtClass)*X$N_age)/sum(X$N_age) # can be inactive if a comparison to the only age-based simulation
s <- sqrt(sum(X$N_age*(log(X$LngtClass)-m)^2)/sum(X$N_age)) # can be inactive if a comparison to the only age-based simulation
r <- round(exp(rnorm(sum(X$N_age),mean = m,sd = s)))-1
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r)))

ggplot(X, aes(x = (LngtClass), y = N_age)) +
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()


# calculate mean and sd length for each year and age
n_years <- length(years) # number of years
n_ages <- length(ages)
byYear <- aggregate(N_age~yd+LngtClass+year,data=flounder,FUN=sum) # data by length, year, and age
# define empty matrices
mat_sd <- matrix(NA,ncol=n_years,nrow = length(ages))
mat_mean <- matrix(NA,ncol=n_years,nrow = length(ages))
mat_n <- matrix(NA,ncol=n_years,nrow = length(ages))
# normal
for (i in 1:n_years){
  for (j in 1:n_ages){
    dada <- byYear %>% filter(yd==ages[j] & year==years[i])
    obs <- dada$LngtClass
    mat_mean[j,i] <- sum(obs*dada$N_age)/sum(dada$N_age)
    mat_sd[j,i] <- sqrt(sum(dada$N_age*(obs-mat_mean[j,i])^2)/sum(dada$N_age))
    mat_n[j,i] <- sum(dada$N_age)
  }
}

mat_sd[mat_sd==0] <- NA # do not consider sd's of 0
sd_red <- mat_sd[-c(24:33),] # remove ages with low number of observations
mean_red <- mat_mean[-c(24:33),]
n_red <- mat_n[-c(24:33),]
# Convert to data frame with row and column indices for plot purposes
dsd <- as.data.frame(sd_red)
dsd$row <- 1:nrow(dsd)
# Pivot to long format
df_sd <- dsd %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_sd, aes(x = col+min(years)-1, y = row, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "SD [cm]") +
  theme_minimal()


# Convert to data frame with row and column indices for plot purposes
dmean <- as.data.frame(mean_red)
dmean$row <- 1:nrow(dmean)
# Pivot to long format
df_mean <- dmean %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_mean, aes(x = col+min(years)-1, y = row, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "Mean [cm]") +
  theme_minimal()
# Convert n to data frame with row and column indices for plot purposes
dn <- as.data.frame(n_red)
dn$row <- 1:nrow(dn)
# Pivot to long format
df_n <- dn %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_mean, aes(x = col+min(years)-1, y = row, fill = log10(value))) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "log10(ind.)") +
  theme_minimal()


mean_sds <- rowMeans(sd_red,na.rm = TRUE)
mean_means <- rowMeans(mean_red,na.rm = TRUE)
n_sds <- apply(n_red,1,FUN=mean,na.rm = TRUE)
#par(mfrow=c(5,5))
#for(i in 1:25) plot(mat_n[i,],mat_sd[i,])
par(mfrow=c(1,1))

plot(ages[-c(24:33)],mean_sds)
plot(mean_means,mean_sds)

df.flounder <- data.frame(yd = ages[-c(24:33)],
                         mean=mean_means,sd=mean_sds,
                         species=rep("flounder",length(mean_means)),
                         n=n_sds) # mean catch pr. quarter
# log-normal if necessary at some point
for (i in 1:n_years){
  for (j in 1:length(ages)){
    dada <- byYear %>% filter(yd==ages[j] & year==years[i])
    Lobs <- log(dada$LngtClass)
    mat_mean[j,i] <- sum(Lobs*dada$N_age)/sum(dada$N_age)
    mat_sd[j,i] <- sqrt(sum(dada$N_age*(Lobs-mat_mean[j,i])^2)/sum(dada$N_age))
    mat_n[j,i] <- sum(dada$N_age)
  }
}
rm(list=setdiff(ls(),c('hl_N','df.cod','df.herring','df.flounder')))
#################

##############
# plaice
#############
plaice <- hl_N %>% filter(species=='plaice' & Age>1) # Look at cod age 2+
plaice$year <- as.numeric(substr(plaice$haulID,1,4)) # redefine year as a variable
plaice$yd <- plaice$Age+plaice$jday # define a julian day pseudo age
ages <- sort(unique(plaice$yd))
years <- sort(unique(plaice$year))

# Group by age only
byAge <- aggregate(N_age~yd+LngtClass,data=plaice,FUN=sum) # group by age
nAge  <- aggregate(N_age~yd,data=plaice,FUN=sum) # Count observations
X <- byAge %>% filter(yd==ages[2]) # choose the age to plot
# normal
m <-  sum(X$LngtClass*X$N_age)/sum(X$N_age) # weighted mean
s <- sqrt(sum(X$N_age*(X$LngtClass-m)^2)/sum(X$N_age)) # weighted sd
r <- round(rnorm(sum(X$N_age),mean=m,sd=s)) # simulated data
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r))) #make simulation a df

ggplot(X, aes(x = (LngtClass), y = N_age)) + # plot the observed and simulated length-frequencies together
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# log-normal
m <-  sum(log(X$LngtClass)*X$N_age)/sum(X$N_age) # weighted mean of log(observations)
s <- sqrt(sum(X$N_age*(log(X$LngtClass)-m)^2)/sum(X$N_age))# weighted sd of log(observations)
r <- round(exp(rnorm(sum(X$N_age),mean = m,sd = s)))-1 # simulated data
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r))) #make simulation a df

ggplot(X, aes(x = (LngtClass), y = N_age)) +# plot the observed and simulated length-frequencies together
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# Group by age and year
byYear <- aggregate(N_age~yd+LngtClass+year,data=plaice,FUN=sum) # group by age and year
X <- byYear %>% filter(yd==ages[1] & year==1997) # select age and year to plot
# normal
m <-  sum(X$LngtClass*X$N_age)/sum(X$N_age) # can be inactive if a comparison to the only age-based simulation
s <- sqrt(sum(X$N_age*(X$LngtClass-m)^2)/sum(X$N_age)) # can be inactive if a comparison to the only age-based simulation
r <- round(rnorm(sum(X$N_age),mean=m,sd=s))
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r)))

ggplot(X, aes(x = (LngtClass), y = N_age)) +
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# log-normal
m <-  sum(log(X$LngtClass)*X$N_age)/sum(X$N_age) # can be inactive if a comparison to the only age-based simulation
s <- sqrt(sum(X$N_age*(log(X$LngtClass)-m)^2)/sum(X$N_age)) # can be inactive if a comparison to the only age-based simulation
r <- round(exp(rnorm(sum(X$N_age),mean = m,sd = s)))-1
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r)))

ggplot(X, aes(x = (LngtClass), y = N_age)) +
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()


# calculate mean and sd length for each year and age
n_years <- length(years) # number of years
n_ages <- length(ages)
byYear <- aggregate(N_age~yd+LngtClass+year,data=plaice,FUN=sum) # data by length, year, and age
# define empty matrices
mat_sd <- matrix(NA,ncol=n_years,nrow = length(ages))
mat_mean <- matrix(NA,ncol=n_years,nrow = length(ages))
mat_n <- matrix(NA,ncol=n_years,nrow = length(ages))
# normal
for (i in 1:n_years){
  for (j in 1:n_ages){
    dada <- byYear %>% filter(yd==ages[j] & year==years[i])
    obs <- dada$LngtClass
    mat_mean[j,i] <- sum(obs*dada$N_age)/sum(dada$N_age)
    mat_sd[j,i] <- sqrt(sum(dada$N_age*(obs-mat_mean[j,i])^2)/sum(dada$N_age))
    mat_n[j,i] <- sum(dada$N_age)
  }
}

mat_sd[mat_sd==0] <- NA # do not consider sd's of 0
sd_red <- mat_sd[-c(27:34),] # remove ages with low number of observations
mean_red <- mat_mean[-c(27:34),]
n_red <- mat_n[-c(27:34),]
# Convert to data frame with row and column indices for plot purposes
dsd <- as.data.frame(sd_red)
dsd$row <- 1:nrow(dsd)
# Pivot to long format
df_sd <- dsd %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_sd, aes(x = col+min(years)-1, y = row, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "SD [cm]") +
  theme_minimal()


# Convert to data frame with row and column indices for plot purposes
dmean <- as.data.frame(mean_red)
dmean$row <- 1:nrow(dmean)
# Pivot to long format
df_mean <- dmean %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_mean, aes(x = col+min(years)-1, y = row, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "Mean [cm]") +
  theme_minimal()
# Convert n to data frame with row and column indices for plot purposes
dn <- as.data.frame(n_red)
dn$row <- 1:nrow(dn)
# Pivot to long format
df_n <- dn %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_mean, aes(x = col+min(years)-1, y = row, fill = log10(value))) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "log10(ind.)") +
  theme_minimal()


mean_sds <- rowMeans(sd_red,na.rm = TRUE)
mean_means <- rowMeans(mean_red,na.rm = TRUE)
n_sds <- apply(n_red,1,FUN=mean,na.rm = TRUE)
#par(mfrow=c(6,5))
#for(i in 1:25) plot(mat_n[i,],mat_sd[i,])
par(mfrow=c(1,1))
plot(mat_n[1,],mat_sd[1,])
plot(ages[-c(27:34)],mean_sds)
plot(mean_means,mean_sds)

df.plaice <- data.frame(yd = ages[-c(27:34)],
                          mean=mean_means,sd=mean_sds,
                          species=rep("plaice",length(mean_means)),
                          n=n_sds) # mean catch pr. quarter
# log-normal if necessary at some point
for (i in 1:n_years){
  for (j in 1:length(ages)){
    dada <- byYear %>% filter(yd==ages[j] & year==years[i])
    Lobs <- log(dada$LngtClass)
    mat_mean[j,i] <- sum(Lobs*dada$N_age)/sum(dada$N_age)
    mat_sd[j,i] <- sqrt(sum(dada$N_age*(Lobs-mat_mean[j,i])^2)/sum(dada$N_age))
    mat_n[j,i] <- sum(dada$N_age)
  }
}
rm(list=setdiff(ls(),c('hl_N','df.cod','df.herring','df.flounder','df.plaice')))
#################

##############
# dab
#############
dab <- hl_N %>% filter(species=='dab' & Age>1) # Look at cod age 2+
dab$year <- as.numeric(substr(dab$haulID,1,4)) # redefine year as a variable
dab$yd <- dab$Age+dab$jday # define a julian day pseudo age
ages <- sort(unique(dab$yd))
years <- sort(unique(dab$year))

# Group by age only
byAge <- aggregate(N_age~yd+LngtClass,data=dab,FUN=sum) # group by age
nAge  <- aggregate(N_age~yd,data=dab,FUN=sum) # Count observations
X <- byAge %>% filter(yd==ages[2]) # choose the age to plot
# normal
m <-  sum(X$LngtClass*X$N_age)/sum(X$N_age) # weighted mean
s <- sqrt(sum(X$N_age*(X$LngtClass-m)^2)/sum(X$N_age)) # weighted sd
r <- round(rnorm(sum(X$N_age),mean=m,sd=s)) # simulated data
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r))) #make simulation a df

ggplot(X, aes(x = (LngtClass), y = N_age)) + # plot the observed and simulated length-frequencies together
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# log-normal
m <-  sum(log(X$LngtClass)*X$N_age)/sum(X$N_age) # weighted mean of log(observations)
s <- sqrt(sum(X$N_age*(log(X$LngtClass)-m)^2)/sum(X$N_age))# weighted sd of log(observations)
r <- round(exp(rnorm(sum(X$N_age),mean = m,sd = s)))-1 # simulated data
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r))) #make simulation a df

ggplot(X, aes(x = (LngtClass), y = N_age)) +# plot the observed and simulated length-frequencies together
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# Group by age and year
byYear <- aggregate(N_age~yd+LngtClass+year,data=dab,FUN=sum) # group by age and year
X <- byYear %>% filter(yd==ages[1] & year==1997) # select age and year to plot
# normal
m <-  sum(X$LngtClass*X$N_age)/sum(X$N_age) # can be inactive if a comparison to the only age-based simulation
s <- sqrt(sum(X$N_age*(X$LngtClass-m)^2)/sum(X$N_age)) # can be inactive if a comparison to the only age-based simulation
r <- round(rnorm(sum(X$N_age),mean=m,sd=s))
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r)))

ggplot(X, aes(x = (LngtClass), y = N_age)) +
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()

# log-normal
m <-  sum(log(X$LngtClass)*X$N_age)/sum(X$N_age) # can be inactive if a comparison to the only age-based simulation
s <- sqrt(sum(X$N_age*(log(X$LngtClass)-m)^2)/sum(X$N_age)) # can be inactive if a comparison to the only age-based simulation
r <- round(exp(rnorm(sum(X$N_age),mean = m,sd = s)))-1
dr <- data.frame(Freq = as.numeric(table(r)),l=sort(unique(r)))

ggplot(X, aes(x = (LngtClass), y = N_age)) +
  geom_col(fill = "steelblue") +geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
  labs(x = "Length Class", y = "Summed N_age") +
  theme_minimal()


# calculate mean and sd length for each year and age
n_years <- length(years) # number of years
n_ages <- length(ages)
byYear <- aggregate(N_age~yd+LngtClass+year,data=dab,FUN=sum) # data by length, year, and age
# define empty matrices
mat_sd <- matrix(NA,ncol=n_years,nrow = length(ages))
mat_mean <- matrix(NA,ncol=n_years,nrow = length(ages))
mat_n <- matrix(NA,ncol=n_years,nrow = length(ages))
# normal
for (i in 1:n_years){
  for (j in 1:n_ages){
    dada <- byYear %>% filter(yd==ages[j] & year==years[i])
    obs <- dada$LngtClass
    mat_mean[j,i] <- sum(obs*dada$N_age)/sum(dada$N_age)
    mat_sd[j,i] <- sqrt(sum(dada$N_age*(obs-mat_mean[j,i])^2)/sum(dada$N_age))
    mat_n[j,i] <- sum(dada$N_age)
  }
}

mat_sd[mat_sd==0] <- NA # do not consider sd's of 0
sd_red <- mat_sd[-c(2,5,8,11,14,22:31),] # remove ages with low number of observations
mean_red <- mat_mean[-c(2,5,8,11,14,22:31),]
n_red <- mat_n[-c(2,5,8,11,14,22:31),]
# Convert to data frame with row and column indices for plot purposes
dsd <- as.data.frame(sd_red)
dsd$row <- 1:nrow(dsd)
# Pivot to long format
df_sd <- dsd %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_sd, aes(x = col+min(years)-1, y = row, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "SD [cm]") +
  theme_minimal()


# Convert to data frame with row and column indices for plot purposes
dmean <- as.data.frame(mean_red)
dmean$row <- 1:nrow(dmean)
# Pivot to long format
df_mean <- dmean %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_mean, aes(x = col+min(years)-1, y = row, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "Mean [cm]") +
  theme_minimal()
# Convert n to data frame with row and column indices for plot purposes
dn <- as.data.frame(n_red)
dn$row <- 1:nrow(dn)
# Pivot to long format
df_n <- dn %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "value"
  ) %>%
  mutate(col = as.integer(gsub("V", "", col)))  # Convert column names to numeric

# Plot the beast
ggplot(df_mean, aes(x = col+min(years)-1, y = row, fill = log10(value))) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x = "Year", y = "Age [years]", fill = "log10(ind.)") +
  theme_minimal()


mean_sds <- rowMeans(sd_red,na.rm = TRUE)
mean_means <- rowMeans(mean_red,na.rm = TRUE)
n_sds <- apply(n_red,1,FUN=mean,na.rm = TRUE)
#par(mfrow=c(4,5))
#for(i in 1:19) plot(mat_n[i,],mat_sd[i,])
#par(mfrow=c(1,1))
plot(mat_n[1,],mat_sd[1,])
plot(ages[-c(2,5,8,11,14,22:31)],mean_sds)
plot(mean_means,mean_sds)


df.dab <- data.frame(yd = ages[-c(2,5,8,11,14,22:31)],
                        mean=mean_means,sd=mean_sds,
                        species=rep("dab",length(mean_means)),
                        n=n_sds) # mean catch pr. quarter
# log-normal if necessary at some point
for (i in 1:n_years){
  for (j in 1:length(ages)){
    dada <- byYear %>% filter(yd==ages[j] & year==years[i])
    Lobs <- log(dada$LngtClass)
    mat_mean[j,i] <- sum(Lobs*dada$N_age)/sum(dada$N_age)
    mat_sd[j,i] <- sqrt(sum(dada$N_age*(Lobs-mat_mean[j,i])^2)/sum(dada$N_age))
    mat_n[j,i] <- sum(dada$N_age)
  }
}
rm(list=setdiff(ls(),c('hl_N','df.cod','df.herring','df.flounder','df.plaice','df.dab')))
#################

df <- rbind(df.cod,df.herring,df.flounder,df.plaice,df.dab)
rm(list=setdiff(ls(),c('df')))
