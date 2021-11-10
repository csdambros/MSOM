## Simulate data for MSOM

source("MixedBinom_wrapper.R")

# Create labels ####
## Create site names ####
lcomb<-expand.grid(LETTERS,LETTERS,LETTERS)
sites<-sample(paste0(lcomb[,1],lcomb[,2],lcomb[,3]))

## Create species names ####
spnames<-read.csv("Data/anura names.csv")
spnames.split<-do.call(rbind,strsplit(spnames$x," "))
spnames3<-gsub("([A-z]{3}).*","\\1",spnames.split)
species<-sample(paste(spnames3[,1],spnames3[,2],sep = "_"))

# 1. No variation and no covariate ####

# Number of species
nsp = 23

# Number of sites
nsites = 30

# Number of replicates per site
J = 12

# Overall probabilities
# Occupancy
psi = 0.3

# Detection
rho = 0.08

# Simulate observation matrix
obs<-matrix(NA,nsites,nsp,dimnames = list(sort(sites[1:nsites]),sort(species[1:nsp])))

occ<-obs
occ[,]<-rbinom(nsites*nsp,1,psi)
obs[,]<-rbinom(occ,occ*J,rho)


# for(j in 1:nsp){
#   # Site occupancies
#   occ <- rbinom(nsites,1,psi)
#   for(i in 1:nsites){
#     # observations (detections)
#     obs[i,j] <- rbinom(1,occ[i]*nrep,rho)
#   }
# }

# Save simulated data
#write.csv(obs,"Data/simu1.csv")

bf<-bFrameOccuMultiBinom(y = obs,J = J)

msom<-occuMSOM(~1+(1|species)~1+(1|species),bf=bf,n.adapt = 500,n.iter = 100,burnin = 50,thin = 10,n.chains = 3,method = "parallel")

coefs<-coefficients(msom)


# Estimated overall psi
plogis(coefs$fixed$occupancy)
psi

# Estimated overall rho
plogis(coefs$fixed$detection)
rho

# Estimated psi per species
plogis(coefs$spcoefs$occupancy)

# Estimated rho per species
plogis(coefs$spcoefs$detection)

# 2. Two species and 4 sites ####
sp1<-c(1,1,1,1)
sp2<-c(1,0,0,1)

occ<-cbind(sp2,sp1)
rownames(occ)<-LETTERS[1:length(sp1)]
J = 30

rho=0.09

# Create observed data
det<-occ
det[,]<-rbinom(occ,occ*J,rho)

# True occurrences
occ

# Observed data
det

# Transform data
bf<-bFrameOccuMultiBinom(y = det,J = J)

# Run model
msom<-occuMSOM(~1+(1|species)~1+(1|species),bf=bf,n.adapt = 500,n.iter = 100,burnin = 50,thin = 10,n.chains = 3,method = "parallel")

# Extract coefficients
coefs<-coefficients(msom)
coefs

# Compare estimates and truth
# Overall occupancy
# estimated
plogis(coefs$fixed$occupancy)
# truth
mean(occ)

# Overall detection
# estimated
plogis(coefs$fixed$detection)
# truth
rho

### Calculate summary statistics (HDI/HDP), Gelman-Rubin statistic
summ<-summary.occuMSOM(msom,returnZ = TRUE)

summ$alphaDiv
summ$gammaDiv

apply(summ$z,c(1,2),min)
apply(summ$z,c(1,2),mean)
(det>0)*1

# 3. Two species along environmental gradient ####

var <- sample(seq(0,10,length.out=30))
var.std<-scale(var)

coefs1<-c(-1,0.7)
coefs2<-c(0,-2)

plot(function(x)plogis(cbind(1,x)%*%coefs1),xlim=range(var.std),ylim=c(0,1))
plot(function(x)plogis(cbind(1,x)%*%coefs2),xlim=range(var.std),ylim=c(0,1),col=2,add=TRUE)


p1<-plogis(cbind(1,var.std)%*%coefs1)
p2<-plogis(cbind(1,var.std)%*%coefs2)


sp1<-rbinom(p1,1,p1)
sp2<-rbinom(p2,1,p2)

occ<-cbind(sp1,sp2)


matplot(var.std,occ,add=TRUE,pch=21,bg=1:2)


J = 30
rho=0.2

# Create observed data
det<-occ
det[,]<-rbinom(occ,occ*J,rho)

# True occurrences
occ

# Observed data
det

# Transform data
bf<-bFrameOccuMultiBinom(y = det,J = J,siteCovs = data.frame(X=var.std))

# Run model
msom<-occuMSOM(~X+(X|species)~X+(X|species),bf=bf,n.adapt = 500,n.iter = 100,burnin = 50,thin = 10,n.chains = 3,method = "parallel")

# Extract coefficients
coefs<-coefficients(msom)
coefs


# Compare estimates and truth
# Occupancy

plot(function(x)plogis(cbind(1,x)%*%coefs1),xlim=range(var.std),ylim=c(0,1),main="Occupancy",ylab=expression(psi))
plot(function(x)plogis(cbind(1,x)%*%coefs2),xlim=range(var.std),ylim=c(0,1),col=2,add=TRUE)


plot(function(x)plogis(cbind(1,x)%*%coefs$spcoefs$occupancy[1,]),xlim=range(var.std),ylim=c(0,1),lty=2,add=TRUE)
plot(function(x)plogis(cbind(1,x)%*%coefs$spcoefs$occupancy[2,]),xlim=range(var.std),ylim=c(0,1),col=2,lty=2,add=TRUE)

legend("topright",
       c("sp1 (true)","sp2 (true)","sp1 (est)","sp2 (est)"),col=1:2,lty=c(1,1,2,2),bty="n")


# Detection
plot(function(x)plogis(cbind(1,x)%*%c(qlogis(rho),0)),xlim=range(var.std),ylim=c(0,1),main="Detection",ylab=expression(rho))
plot(function(x)plogis(cbind(1,x)%*%c(qlogis(rho),0)),xlim=range(var.std),ylim=c(0,1),col=2,add=TRUE)


plot(function(x)plogis(cbind(1,x)%*%coefs$spcoefs$detection[1,]),xlim=range(var.std),ylim=c(0,1),lty=2,add=TRUE)
plot(function(x)plogis(cbind(1,x)%*%coefs$spcoefs$detection[2,]),xlim=range(var.std),ylim=c(0,1),col=2,lty=2,add=TRUE)

legend("topright",
       c("sp1 (true)","sp2 (true)","sp1 (est)","sp2 (est)"),col=1:2,lty=c(1,1,2,2),bty="n")



# estimated
plogis(coefs$fixed$occupancy)
# truth
mean(occ)

# Overall detection
# estimated
plogis(coefs$fixed$detection)
# truth
rho

### Calculate summary statistics (HDI/HDP), Gelman-Rubin statistic
summ<-summary.occuMSOM(msom,returnZ = TRUE)

summ$fixed$occupancy
summ$spdeviations$occupancy

summ$alphaDiv
summ$gammaDiv

# 4. Species richness estimator ####

nsp<-20

var <- sample(seq(0,10,length.out=30))
var.std<-scale(var)

coefs1<-c(-1,0.7)
sd<-c(0.4,1)
J = 30
rho=0.04

coefs.all<-matrix(NA,2,nsp)

for(i in 1:nsp){
  coefs.all[,i]<-rnorm(coefs1,coefs1,sd)
}

p.all<-plogis(cbind(1,var.std)%*%coefs.all)

occ<-pall
occ[,]<-rbinom(p.all,1,p.all)


plot(function(x)plogis(cbind(1,x)%*%coefs1),xlim=range(var.std),ylim=c(0,1))


# Create observed data
det<-occ
det[,]<-rbinom(occ,occ*J,rho)

matplot(sort(var.std),p.all[order(var.std),],type="l",lty=1)




# True occurrences
occ

# Observed data
det

# Number of species (true)
sum(colSums(occ)>0)

# Number of species observed
sum(colSums(det)>0)

# Augment data
det.aug<-data.frame(det,0,0,0,0,0,0,0)

# Transform data
bf<-bFrameOccuMultiBinom(y = det.aug,J = J,siteCovs = data.frame(X=var.std))

# Run model
msom<-occuMSOM(~X+(X|species)~X+(X|species),bf=bf,n.adapt = 500,n.iter = 100,burnin = 50,thin = 10,n.chains = 3,method = "parallel")

# Extract coefficients
coefs<-coefficients(msom)
coefs


# Compare estimates and truth
# Occupancy
plot(function(x)plogis(cbind(1,x)%*%coefs1),xlim=range(var.std),ylim=c(0,1),main="Occupancy",ylab=expression(psi))

plot(function(x)plogis(cbind(1,x)%*%coefs$fixed$occupancy[1,]),xlim=range(var.std),ylim=c(0,1),lty=2,add=TRUE)


####
p.all.est<-plogis(cbind(1,var.std)%*%t(coefs$spcoefs$occupancy))

matplot(sort(var.std),p.all[order(var.std),],type="l",lty=1)

matplot(sort(var.std),p.all.est[order(var.std),],type="l",lty=2,add=TRUE)


matplot(sort(var.std),p.all[order(var.std),1:5],type="l",lty=1)

matplot(sort(var.std),p.all.est[order(var.std),1:5],type="l",lty=2,add=TRUE)

# Detection
plot(function(x)plogis(cbind(1,x)%*%c(qlogis(rho),0)),xlim=range(var.std),ylim=c(0,1),main="Detection",ylab=expression(rho))

plot(function(x)plogis(cbind(1,x)%*%coefs$fixed$detection[1,]),xlim=range(var.std),ylim=c(0,1),lty=2,add=TRUE)


# estimated
coefs$fixed$occupancy
coefs1

# Overall detection
# estimated
coefs$fixed$detection
# truth
qlogis(rho)

### Calculate summary statistics (HDI/HDP), Gelman-Rubin statistic
summ<-summary.occuMSOM(msom,returnZ = TRUE)

summ$fixed$occupancy
summ$spdeviations$occupancy

coefs$spcoefs$occupancy

summ$alphaDiv

# True vs. observed
plot(rowSums(occ),rowSums(det>0))
abline(0,1)

# True vs. estimated
plot(rowSums(occ),summ$alphaDiv$Median,col=2)
abline(0,1)


summ$gammaDiv


## Change in species composition ####
library(vegan)

jac.true<-1-vegdist(occ>0,method="jac")
jac.obs<-1-vegdist(det>0,method="jac")
jac<-1-apply(summ$z,3,vegdist,method="jac")

jac.est<-apply(jac,1,median,na.rm=TRUE)

plot(jac.true,jac.obs,pch=21,bg="grey80",col="grey80")
points(jac.true,jac.est,pch=21,bg=1)
abline(0,1)
abline(lm(jac.obs~jac.true),lty=2)
abline(lm(jac.est~jac.true),lty=3)
cor(jac.obs,jac.true,use = "complete")
cor(jac.est,jac.true,use = "complete")

plot(jac.true,jac.obs,pch=21,bg="grey80",col="grey80")
abline(0,1)
arrows(jac.true,jac.obs,jac.true,jac.est,length = 0.03)




apply(summ$z,c(1,2),min)
apply(summ$z,c(1,2),mean)
(det>0)*1
