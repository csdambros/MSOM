# Script created by CSDambros
# Bayesian inference to estimate species detection, occurrence, and unbiased species richness based on species counts

#--------------- Load R packages ---------------

# requires the jags program to be installed in computer (visit https://sourceforge.net/projects/mcmc-jags/files/)
# requires the rjags package to be installed in R (use install.packages ("rjags"))
library(rjags)
library(vegan)

# library(runjags) # Only used for fast parallel processing (not required; see below)


## The model DorazioMSOM.jags in the JAGS folder is based on Dorazio and Royle 2005 and Dorazio et al. 2006
#Code adapted from http://www.mbr-pwrc.usgs.gov/site/communitymodeling/basic%20model.zip


###############################################

#--------------- Using model with empirical data ---------------

###############################################

# 1. Load data from Dambros et al. (2016)

# Termite data: Each row is a record (vial). Data transformation is required to run model. The following code shows how this tranformation is done.
termiteLongData<-read.csv("https://ndownloader.figshare.com/files/3299687")

# Environmental data: Each entry is a value of soil clay content in a transect
termiteEnvData<-read.csv("https://ndownloader.figshare.com/files/3299696")

# Create object with standardized environmental variable
E<-termiteEnvData$Clay
E.std<-as.vector(decostand(E,"standardize"))

# For each species, create a table with sites in rows and sections in columns
# Jags code was writen to receive an array (3d matrix with sites x sections x species)
# Before running the model, you need to have this array named here as detmat

# However, you will not save your data as a 3d matrix. Instead, you have the option to 1) save each individual record as a row (long table) or 2) save a site x species table (short format) in which the number of presences in each site (row) for each species (columns) fill in the table.

# If the data is organized in the long format (as imported here, prefereable), use the FIRST method to create the detmat table. If your data is organized in the short format (site x species), then use the SECOND method (see below).

# FIRST METHOD
counts<-termiteLongData$Frequency # Extract species counts
N<-nlevels(termiteLongData$Taxon) # Count number of species in dataset

transect<-paste(termiteLongData$Trail,termiteLongData$Plot,sep ="_")
section<-termiteLongData$Distance_beginning
species<-termiteLongData$Taxon

# Create array 
detmat<-tapply(counts,list(transect,section,species),sum)
detmat[is.na(detmat)]<-0

# Add fictitious species that could have occurred (with zero records)
nzero<-300 # number of species to add. Any large number certainly larger than estimated

# Add fictitious species with zero occurrences to array (required)
detmat<-array(detmat,dim=dim(detmat)+c(0,0,nzero))
detmat[,,(N+1):(N+nzero)]<-0

### END OF TRANSFORMATION USING FIRST METHOD

# SECOND METHOD

#################################################################
## If you have a dataframe with sites x species filled with the number of sections in which a species was found, use the following code to create the detmat array

# We will first convert the long format into the site x species matrix to show how it works. This is not required if you have the long format and already have the detmat array

# Create transect x species matrix to show how it works
counts<-termiteLongData$Frequency
N<-nlevels(termiteLongData$Taxon)

transect<-paste(termiteLongData$Trail,termiteLongData$Plot,sep ="_")
section<-termiteLongData$Distance_beginning
species<-termiteLongData$Taxon

spDet<-tapply(counts,list(transect,species),sum) # Create matrix (spDet is the name of your dataset)
spDet[is.na(spDet)]<-0

# If you have the short format, you would only need to import here and name it spDet


### START OF TRANSFORMATION
# The transformation assumes all transects have the same number of sections, but this can be modified
K<-10 # maximum number of sections per transect
N<-ncol(spDet) # Number of observed species
J<-nrow(spDet) # Number of transects surveyed

# Add fictitious species that could have occurred (with zero records)
nzero<-300 #Just any large number

detmat<-array(NA,dim=c(J,K,N+nzero))# Create 3d matrix (array) to receive data from individual sections
detmat[,,(N+1):(N+nzero)]<-0 # Set detection of unnobserved species to zero

# Insert detection history for observed species
for(s in 1:N){
  detmat[,,s]<-(matrix(1:K,nrow=J,ncol=K,byrow = TRUE)<=spDet[,s])*1
}

### END OF TRANSFORMATION USING SECOND METHOD
######################################################################

# 2. Estimate species richness using bayesian inference
# Load all the data into a list
sp.data = list(n=N, nzeroes=nzero, J=J, K=rep(K,J), X=detmat,Env=E.std)

# Define parameters to track during MCMC
sp.params = c('u', 'v','b1.psi','b1.p', 'mu.u', 'mu.v', 'tau.u', 'tau.v', 'omega'
              , 'N','N3','mu.b1.psi','mu.b1.p','tau.b1.psi','tau.b1.p')

# Set initial parameters (makes it easier for the model to start properly but is not required)
sp.inits = function() {
  omegaGuess = runif(1, N/(N+nzero), 1)
  psi.meanGuess = runif(1, .25,1)
  list(omega=omegaGuess,w=c(rep(1, N), rbinom(nzero, size=1, prob=0.5)),
       u=rnorm(N+nzero), v=rnorm(N+nzero),
       Z = cbind(spDet*0+1,spDet[,rep(1,nzero)]*0)
  )
}

#Run the model and call the results

#Normal use of jags (caution: can take over 2 hours to run in a ordinary computer!!)

# Run burn in with 5000 iterations
fit = jags.model(file = "basicmodel.jags",inits = sp.inits,
                 data = sp.data, n.chains = 3,n.adapt = 10000)

# Run model with 10000 iterations x 3 chains. Thinnig reduces amount of stored data
coda.resu<-coda.samples(fit,variable.names = sp.params,thin = 5,n.iter = 20000)


#Running in parallel (multi core), requires runjags package (1 chain = 1 core)
#Uncomment lines to use

## Use detectCores from the parallel package if number of cores is not known
## Default = use 3 cores (quadcore pc)
## Leave one core unused to avoid computer to crash
fit = run.jags(model = "basicmodel.jags",inits = sp.inits,
               data = sp.data, n.chains = 3,monitor = sp.params,
               method = "parallel",thin = 3)
coda.resu<-as.mcmc.list(fit)

# 3. See a summary of parameter estimates
coda.matrix<-as.data.frame(as.matrix(coda.resu))
summary(coda.matrix)

#See estimates of total richness in metacommunity (N)
BN = coda.matrix[,"N"]
summary(BN)
hist(BN,freq = FALSE,main="",ylab="Posterior probability",xlab="Species richness",cex.lab=1.2,col=1,breaks = 30)
box()
abline(v=quantile(BN,c(0.025,0.975)),col=2,lwd=2)
legend("bottom",legend = "95% HDI",bty = "n",text.col = 2,cex=2)

#See estimates of richness per transect (and 95% Highest Density Interval - HDI)
BNt = coda.matrix[,grepl("N3\\[",colnames(coda.matrix),perl=TRUE)]
BNt.median<-apply(BNt,2,quantile,c(0.025,0.5,0.975))
matplot(sort(E),t(BNt.median[,order(E)]),type="l",ylab="Species richness",ylim=c(0,N+nzero),cex.lab=1.2,xlab="Soil clay content",lty=1,col=c(2,1,2))

median(BNt.median)
mean(BNt.median)

# 3.1 Extract overall detection and occurrence data (mean of all species)

# Overall species occurrence (hyperparameter)
Bb0 = coda.matrix[,"mu.u"]
summary(Bb0)
hist(Bb0,freq = FALSE)

# Overall species detection (hyperparameter)
Bc0 = coda.matrix[,"mu.v"]
summary(Bc0)
hist(Bc0,freq = FALSE)

# Overall effect of environment on species occurrence (hyperparameter)
Bb1 = coda.matrix[,"mu.b1.psi"]
summary(Bb1)
hist(Bb1,freq = FALSE)

summary(1/(1+exp(-Bb1)))# Prob scale (logit transformed)
hist(1/(1+exp(-Bb1)),freq = FALSE)

# Overall effect of environment on species detection (hyperparameter)
Bc1 = coda.matrix[,"mu.b1.p"]
summary(Bc1)
hist(Bc1,freq = FALSE)

summary(1/(1+exp(-Bc1)))# Prob scale (logit transformed)
hist(1/(1+exp(-Bc1)),freq = FALSE)

# Variation (Standard deviation) across species in occurrence (hyperparameters)
Bb0sd = sqrt(1/coda.matrix[,"tau.u"])
summary(Bb0)
hist(Bb0,freq = FALSE)

# Variation (Standard deviation) across species in detection (hyperparameters)
Bc0sd = sqrt(1/coda.matrix[,"tau.v"])
summary(Bc0)
hist(Bc0,freq = FALSE)

# Variation (Standard deviation) across species in how occurrence responds to environment (hyperparameters)
Bb1sd = sqrt(1/coda.matrix[,"tau.b1.psi"])
summary(Bb1)
hist(Bb1,freq = FALSE)

# Variation (Standard deviation) across species in how detection responds to environment (hyperparameters)
Bc1sd = sqrt(1/coda.matrix[,"tau.b1.p"])
summary(Bc1)
hist(Bc1,freq = FALSE)

#####
# 3.2 Extract occurrence and detection data for individual species

# Individual species occurrence
Bb0s = coda.matrix[,grepl("u\\[",colnames(coda.matrix),perl=TRUE)]
summary(Bb0s)
summary(1/(1+exp(-Bb0s)))# Prob scale (logit transformation)

# Individual species detection
Bc0s = coda.matrix[,grepl("v\\[",colnames(coda.matrix),perl=TRUE)]
summary(Bc0s)
summary(1/(1+exp(-Bc0s)))# Prob scale (logit transformation)

# Effect of environment on individual species occurrences
Bb1s = coda.matrix[,grepl("b1.psi\\[",colnames(coda.matrix),perl=TRUE)]
summary(Bb1s)
summary(1/(1+exp(-Bb1s)))# Prob scale (logit transformation)

# Effect of environment on individual species detection
Bc1s = coda.matrix[,grepl("b1.p\\[",colnames(coda.matrix),perl=TRUE)]
summary(Bc1s)
summary(1/(1+exp(-Bc1s)))# Prob scale (logit transformation)

## Plot occurrence and detection for first 10 species
# Extract median of parameter estimates for each individual species

Bb0sMed<-apply(Bb0s,2,quantile,0.5)
Bb1sMed<-apply(Bb1s,2,quantile,0.5)

Bc0sMed<-apply(Bc0s,2,quantile,0.5)
Bc1sMed<-apply(Bc1s,2,quantile,0.5)

# Calculate occurrence along enviornmental gradient
Bpsi.logit<-cbind(1,E.std)%*%rbind(Bb0sMed,Bb1sMed)
Bpsi<-1/(1+exp(-Bpsi.logit)) # Logit transformation

# Calculate detection along enviornmental gradient
Btheta.logit<-cbind(1,E.std)%*%rbind(Bc0sMed,Bc1sMed)
Btheta<-1/(1+exp(-Btheta.logit)) # Logit transformation


# Plot curves for original and estimated curves
par(mar=c(5,5,1,1),mfrow=c(2,2))

matplot(sort(E),Bpsi[order(E),1:80],type="l",ylim=c(0,1),ylab=expression(psi[s] (estimated)),cex.lab=1.4,xlab="Soil clay content")
matplot(sort(E),Btheta[order(E),1:80],type="l",ylim=c(0,1),ylab=expression(theta[s] (estimated)),cex.lab=1.4,xlab="Soil clay content")

boxplot(Bpsi[,order(-colMeans(Bpsi[,1:N]))],ylab=expression(psi[s] (estimated)),cex.lab=1.4,xaxt="n",ylim=c(0,1))
axis(1,at=1:N,labels=colnames(spDet)[order(-colMeans(Btheta[,1:N]))],las=2,cex.axis=.3)

boxplot(Btheta[,order(-colMeans(Btheta[,1:N]))],ylab=expression(theta[s] (estimated)),cex.lab=1.4,xaxt="n",ylim=c(0,1))
axis(1,at=1:N,labels=colnames(spDet)[order(-colMeans(Btheta[,1:N]))],las=2,cex.axis=.3)

### Analyze a single species
# Syntermes molestus
SmolNum<-match("Syntermes_molestus",colnames(spDet))
samples<-sample(nrow(Bb0s),100)

Bb0Smol<-Bb0s[samples,match("Syntermes_molestus",colnames(spDet))]
Bb1Smol<-Bb1s[samples,match("Syntermes_molestus",colnames(spDet))]

Bc0Smol<-Bc0s[samples,match("Syntermes_molestus",colnames(spDet))]
Bc1Smol<-Bc1s[samples,match("Syntermes_molestus",colnames(spDet))]

# Calculate occurrence along enviornmental gradient
BpsiSmol.logit<-cbind(1,E.std)%*%rbind(Bb0Smol,Bb1Smol)
BpsiSmol<-1/(1+exp(-BpsiSmol.logit)) # Logit transformation

# Calculate detection along enviornmental gradient
BthetaSmol.logit<-cbind(1,E.std)%*%rbind(Bc0Smol,Bc1Smol)
BthetaSmol<-1/(1+exp(-BthetaSmol.logit)) # Logit transformation

par(mar=c(5,5,1,1),mfrow=c(1,2))
matplot(sort(E),BpsiSmol[order(E),],type="l",ylim=c(0,1),ylab=expression(psi[s] (estimated)),cex.lab=1.4,xlab="Soil clay content",lty=1,col=adjustcolor(2,0.02))
points(sort(E),Bpsi[order(E),match("Syntermes_molestus",colnames(spDet))],type="l",ylim=c(0,1),ylab=expression(psi[s] (estimated)),cex.lab=1.4,xlab="Soil clay content",col=1)

matplot(sort(E),BthetaSmol[order(E),],type="l",ylim=c(0,1),ylab=expression(theta[s] (estimated)),cex.lab=1.4,xlab="Soil clay content",lty=1,col=adjustcolor(2,0.02))
points(sort(E),Btheta[order(E),match("Syntermes_molestus",colnames(spDet))],type="l",ylim=c(0,1),ylab=expression(psi[s] (estimated)),cex.lab=1.4,xlab="Soil clay content",col=1)



###############################################

#------ Using model with simulated data ------#

###############################################

# 1. Simulate data

N <- 100 # number of species to be simulated
J <- 30 # Number of transects to be simulated
K <- 20 # Number of sections per transect to be simulated

# Create environmental gradient
E<-sort(runif(J,min = 0, max =100)) # J uniform values from 0 to 100
E.std<-as.vector(decostand(E,"standardize")) # Standardize variable (makes easier for bayesian computation)

# Simulate occurrence and detection for 100 spp
b0 <- rnorm (N, 5, .1) # Simulate mean occurrence for 100 spp (logit scale)
b1 <- rnorm (N, -0.08,0.01) # Simulate effect of env on occurrence for 100 spp (logit scale)

c0 <- rnorm (N, -3, 0.1) # Simulate mean detection for 100 spp (logit scale)
c1 <- rnorm (N, 0, 0.01) # Simulate effect of env on detection for 100 spp (logit scale)

# Calculate occurrence along enviornmental gradient
psi.logit<-cbind(1,E)%*%rbind(b0,b1)
psi<-1/(1+exp(-psi.logit)) # Logit transformation

# Calculate detection along enviornmental gradient
theta.logit<-cbind(1,E)%*%rbind(c0,c1)
theta<-1/(1+exp(-theta.logit)) # Logit transformation

# Visualize first 10 spp
par(mar=c(5,5,1,1),mfrow=c(1,2))

matplot(E,psi[,1:10],type="l",ylab=expression(psi[s]),cex.lab=1.4,xlab="Gradient")
matplot(E,theta[,1:10],type="l",ylab=expression(theta[s]),cex.lab=1.4,xlab="Gradient")

# Simulate occurrence for all species in all transects
spOcc<-matrix(rbinom(psi,1,psi),nrow(psi),ncol(psi))

# Simulate detection of individual species in n out of K sections (when absent from transect then zero, otherwise depends on theta)
spDet<-matrix(rbinom(theta,K*spOcc,theta),nrow(theta),ncol(theta))

# Create data array (dimension 1 (rows) = transects, dimentsion 2 (colmns) = sections, dimension 3 = species)
nzero<-100 #Number of potential unobserved species (just any number certainly higher than expected)

detmat<-array(NA,dim=c(J,K,N+nzero))# Create 3d matrix (array) to receive data from individual sections
detmat[,,(N+1):(N+nzero)]<-0 # Set detection of unnobserved species to zero

# Insert detection history for observed species
for(s in 1:N){
  detmat[,,s]<-(matrix(1:K,nrow=J,ncol=K,byrow = TRUE)<=spDet[,s])*1
}

# 2. Estimate species richness using bayesian inference
# Load all the data into a list
sp.data = list(n=N, nzeroes=nzero, J=J, K=rep(K,J), X=detmat,Env=E.std)

# Define parameters to track during MCMC
sp.params = c('u', 'v','b1.psi','b1.p', 'mu.u', 'mu.v', 'tau.u', 'tau.v', 'omega', 'N','N3','mu.b1.psi','mu.b1.p','tau.b1.psi','tau.b1.p')

# Set initial parameters (makes it easier for the model to start properly but is not required)
sp.inits = function() {
  omegaGuess = runif(1, N/(N+nzero), 1)
  psi.meanGuess = runif(1, .25,1)
  list(omega=omegaGuess,w=c(rep(1, N), rbinom(nzero, size=1, prob=0.5)),
       u=rnorm(N+nzero), v=rnorm(N+nzero),
       Z = cbind(spDet*0+1,spDet*0)
  )
}


#Run the model and call the results

#Normal use of jags (caution: can take over 2 hours to run in a regular computer!!)

# Run burn in with 5000 iterations
fit = jags.model(file = "basic model/basicmodel.jags",inits = sp.inits,
                 data = sp.data, n.chains = 3,n.adapt = 10000)

# Run model with 10000 iterations x 3 chains. Thinnig reduces amount of stored data
coda.resu<-coda.samples(fit,variable.names = sp.params,thin = 5,n.iter = 20000)


#Running in parallel (multi core), requires runjags package (1 chain = 1 core)
#Uncomment lines to use

## Use detectCores from the parallel package if number of cores is not known
## Leave one core unused to avoid computer to crash
fit = run.jags(model = "basic model/basicmodel.jags",inits = sp.inits,
               data = sp.data, n.chains = 3,monitor = sp.params,
               method = "parallel",thin = 3)
coda.resu<-as.mcmc.list(fit)

# 3. See a summary of parameter estimates
coda.matrix<-as.data.frame(as.matrix(coda.resu))
summary(coda.matrix)

#See estimates of total richness in metacommunity (N)
BN = coda.matrix[,"N"]
summary(BN)
hist(BN,freq = FALSE)

#See estimates of richness per transect (and 95% Highest Density Interval - HDI)
BNt = coda.matrix[,grepl("N3\\[",colnames(coda.matrix),perl=TRUE)]
BNt.median<-apply(BNt,2,quantile,c(0.025,0.5,0.975))
matplot(E,t(BNt.median),type="l",ylab="Species richness")


# 3.1 Extract overall detection and occurrence data (mean of all species)

# Overall species occurrence (hyperparameter)
Bb0 = coda.matrix[,"mu.u"]
summary(Bb0)
hist(Bb0,freq = FALSE)

# Overall species detection (hyperparameter)
Bc0 = coda.matrix[,"mu.v"]
summary(Bc0)
hist(Bc0,freq = FALSE)

# Overall effect of environment on species occurrence (hyperparameter)
Bb1 = coda.matrix[,"mu.b1.psi"]
summary(Bb1)
hist(Bb1,freq = FALSE)

summary(1/(1+exp(-Bb1)))# Prob scale (logit transformed)
hist(1/(1+exp(-Bb1)),freq = FALSE)

# Overall effect of environment on species detection (hyperparameter)
Bc1 = coda.matrix[,"mu.b1.p"]
summary(Bc1)
hist(Bc1,freq = FALSE)

summary(1/(1+exp(-Bc1)))# Prob scale (logit transformed)
hist(1/(1+exp(-Bc1)),freq = FALSE)

# Variation (Standard deviation) across species in occurrence (hyperparameters)
Bb0sd = sqrt(1/coda.matrix[,"tau.u"])
summary(Bb0)
hist(Bb0,freq = FALSE)

# Variation (Standard deviation) across species in detection (hyperparameters)
Bc0sd = sqrt(1/coda.matrix[,"tau.v"])
summary(Bc0)
hist(Bc0,freq = FALSE)

# Variation (Standard deviation) across species in how occurrence responds to environment (hyperparameters)
Bb1sd = sqrt(1/coda.matrix[,"tau.b1.psi"])
summary(Bb1)
hist(Bb1,freq = FALSE)

# Variation (Standard deviation) across species in how detection responds to environment (hyperparameters)
Bc1sd = sqrt(1/coda.matrix[,"tau.b1.p"])
summary(Bc1)
hist(Bc1,freq = FALSE)

#####
# 3.2 Extract occurrence and detection data for individual species

# Individual species occurrence
Bb0s = coda.matrix[,grepl("u\\[",colnames(coda.matrix),perl=TRUE)]
summary(Bb0s)
summary(1/(1+exp(-Bb0s)))# Prob scale (logit transformation)

# Individual species detection
Bc0s = coda.matrix[,grepl("v\\[",colnames(coda.matrix),perl=TRUE)]
summary(Bc0s)
summary(1/(1+exp(-Bc0s)))# Prob scale (logit transformation)

# Effect of environment on individual species occurrences
Bb1s = coda.matrix[,grepl("b1.psi\\[",colnames(coda.matrix),perl=TRUE)]
summary(Bb1s)
summary(1/(1+exp(-Bb1s)))# Prob scale (logit transformation)

# Effect of environment on individual species detection
Bc1s = coda.matrix[,grepl("b1.p\\[",colnames(coda.matrix),perl=TRUE)]
summary(Bc1s)
summary(1/(1+exp(-Bc1s)))# Prob scale (logit transformation)

## Plot occurrence and detection for first 10 species
# Extract median of parameter estimates for each individual species

Bb0sMed<-apply(Bb0s,2,quantile,0.5)
Bb1sMed<-apply(Bb1s,2,quantile,0.5)

Bc0sMed<-apply(Bc0s,2,quantile,0.5)
Bc1sMed<-apply(Bc1s,2,quantile,0.5)

# Calculate occurrence along enviornmental gradient
Bpsi.logit<-cbind(1,E.std)%*%rbind(Bb0sMed,Bb1sMed)
Bpsi<-1/(1+exp(-Bpsi.logit)) # Logit transformation

# Calculate detection along enviornmental gradient
Btheta.logit<-cbind(1,E.std)%*%rbind(Bc0sMed,Bc1sMed)
Btheta<-1/(1+exp(-Btheta.logit)) # Logit transformation


# Plot curves for original and estimated curves
par(mar=c(5,5,1,1),mfrow=c(2,2))

matplot(E,psi[,1:10],type="l",ylim=c(0,1),ylab=expression(psi[s] (true)),cex.lab=1.4,xlab="Gradient")
matplot(E,theta[,1:10],type="l",ylim=c(0,.2),ylab=expression(theta[s] (true)),cex.lab=1.4,xlab="Gradient")

matplot(E,Bpsi[,1:10],type="l",ylim=c(0,1),ylab=expression(psi[s] (estimated)),cex.lab=1.4,xlab="Gradient")
matplot(E,Btheta[,1:10],type="l",ylim=c(0,.2),ylab=expression(theta[s] (estimated)),cex.lab=1.4,xlab="Gradient")


