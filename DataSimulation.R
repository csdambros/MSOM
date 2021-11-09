## Simulate data for MSOM

# Create labels ####
# Create site names
ver<-expand.grid(LETTERS,LETTERS,LETTERS)
sites<-sample(paste0(ver[,1],ver[,2],ver[,3]))

## Create species names
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
nrep = 12

# Overall probabilities
# Occupancy
psi = 0.3

# Detection
rho = 0.08

# Simulate observation matrix
obs<-matrix(NA,nsites,nsp,dimnames = list(sort(sites[1:nsites]),sort(species[1:nsp])))

for(j in 1:nsp){
  # Site occupancies
  occ <- rbinom(nsites,1,psi)
  for(i in 1:nsites){
    # observations (detections)
    obs[i,j] <- rbinom(1,occ[i]*nrep,rho)
  }
}

# Save simulated data
#write.csv(obs,"Data/simu1.csv")


bf<-bFrameOccuMultiBinom(y = obs,J = nrep)

msom<-occuMSOM(~1+(1|species)~1+(1|species),bf=bf,n.adapt = 500,n.iter = 100,burnin = 50,thin = 10,n.chains = 3,method = "parallel")

coefs<-coefficients(msom)


# Estimated overall psi
plogis(coefs$fixed$occupancy)

# Estimated overall rho
plogis(coefs$fixed$detection)

# Estimated psi per species
plogis(coefs$spcoefs$occupancy)

# Estimated rho per species
plogis(coefs$spcoefs$detection)


