# Tutorial MSOM Binom (lme4 syntaxe)

# requires runjags package
# To install runjags:
# 1. Install the JAGS program in your computer:
# https://sourceforge.net/projects/mcmc-jags/files/
# 2. Run install.packages("runjags")

# requires the reshape2 package
# To install reshape2:
# 1. Run install.packages("reshape2")


# Load packages ####
library(runjags)
library(reshape2)

# Install MSOM functions ####
source("MixedBinom_wrapper.R")

# Import data ####
occurrences<-read.csv("Data/simu1.csv",row.names = 1)

head(occurrences)


# These data represent a simulation in which frogs were sampled in a set of locations using audio recorders. 12 recorders were placed in each of 30 locations (rows). The data has information on the number of times each species (columns) was observed in each location (from zero to 12). The maximum number of observations was 4 (out of 12 possible).

# rows = sites/locations
# columns = species
# matrix fill = values ranging from 0 - 12
## 0 = species not observed in a site
## 1 = species was observed in the site in only 1 recorder
## 12 = species was observed in all recorders in the site

# Step 1. Transform data ####
# Note, J is the number of samples taken in each site. This number must be >1 and must be provided 
bf<-bFrameOccuMultiBinom(y = occurrences,J = 12,siteCovs = data.frame(var=scale(1:30),var2=sample(scale(1:30))))

head(bf)

# Run model (no covariates)
msom<-occuMSOM(~var+var2+(var2|species)~1+(1|species),bf=bf,n.adapt = 500,n.iter = 100,burnin = 50,thin = 10,method = "parallel")



# Warning messages can be ignored

# Extract coefficients from model
coefs<-coefficients(msom)

# Estimated overall psi (logit scale)
coefs$fixed$occupancy
plogis(coefs$fixed$occupancy)

# Estimated overall rho
coefs$fixed$detection
plogis(coefs$fixed$detection[1])

# Estimated psi per species
coefs$spcoefs$occupancy
plogis(coefs$spcoefs$occupancy)

# Estimated rho per species
coefs$spcoefs$detection
plogis(coefs$spcoefs$detection)

# Check HDI/HPD ("confidence intervals")
HDI<-summary.occuMSOM(msom)

plot(rowSums(occurrences>0),HDI$alphaDiv$Median)
abline(0,1)





head(ver$fixed$detection)
coefs$fixed$detection



ver$alphaDiv




sim1df<-as.data.frame(as.matrix(msom$mcmc))

# Logis scale
hist(sim1df$psi0)


# Converted to probabilities
hist(plogis(sim1df$psi0))
abline(v=quantile(plogis(sim1df$psi0),c(0.025,0.5,0.975)))








ver<-summary(msom)
# Convert results to data frame
sim1df<-as.data.frame(as.matrix(msom$mcmc))

msomcopy<-msom

summary(msom,vars="psi0")


# Get variable names
names<-colnames(sim1df)

### Get results for fixed effects
Fdetection<-sim1df[,grepl("rho0",names),drop=FALSE]
Foccupancy<-sim1df[,grepl("psi0",names),drop=FALSE]



ver<-summary(msom)

ver








