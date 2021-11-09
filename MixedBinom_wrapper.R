

bFrameOccuMultiBinom<-function(y=NULL,siteCovs=NULL,obsCovs=NULL,spCovs=NULL,J=NULL){
  require(reshape2)
  
  if(!is.null(obsCovs)){
    warning("Observation covariates not implemented in this model. Variables will not be included in the output.")
  }
  
  
  call<-match.call(expand.dots = TRUE)
  call.ls<-as.list(call)
  
  if(!is.data.frame(siteCovs)&!is.null(siteCovs)){
    #warning("siteCovs converted to data.frame\n")
    siteCovs<-as.data.frame(siteCovs)
    colnames(siteCovs)<-as.character(call.ls$siteCovs)
  }
  if(!is.data.frame(spCovs)&!is.null(spCovs)){
    #warning("spCovs converted to data.frame\n")
    spCovs<-as.data.frame(spCovs)
    colnames(spCovs)<-as.character(call.ls$spCovs)
  }
  
  if(length(dim(y))==2){
    y<-as.matrix(y)
    occMelt<-melt(y,varnames = c("site","species"),value.name = "y")
    
    if(is.null(J)){
      warning("Number of replicates not provided, assuming J = max (y). Consider providing J.")
      J = max(y)
    }
    
  }
  if(length(dim(y))==3){
    occMelt_pre<-melt(y,varnames = c("site","replicate","species"),value.name = "y")
    occMelt<-aggregate(data.frame(y=occMelt_pre$y),list(site=occMelt_pre$site,species=occMelt_pre$species),sum,na.rm=TRUE)
    
    J=as.vector(apply(y,c(1,3),function(x)sum(!is.na(x))))
    
  }
  
  bframe<-list(occMelt,
               siteCovs[occMelt$site,,drop=FALSE],
               spCovs[occMelt$species,,drop=FALSE],J=J)
  bframe<-bframe[lengths(bframe)>0]
  
  bframe<-as.data.frame(bframe)
  
  return(bframe)
}






occuMSOM<-function(
  formula=~1+(1|species)~1+(1|species),
  bf=bf,
  n.adapt=500,
  n.chains=3,
  inits=NULL,
  n.iter=500,
  method="parallel",
  burnin=100,thin=1){
  
  
  require(runjags)
  require(rjags)
  
  # Parse formulas
  detformula <-as.formula(paste("y",c(formula[[2]]),collapse=""))
  psiformula <- as.formula(paste("z","~", formula[3], sep = ""))
  
  detterms <- gsub("[[:space:]]", "", attr(terms(detformula), "term.labels"))
  psiterms <- gsub("[[:space:]]", "", attr(terms(psiformula), "term.labels"))
  
  detresponse <- as.character(detformula)[2]
  psiresponse <- as.character(psiformula)[2]
  
  detintercept <- attr(terms(detformula), "intercept")
  psiintercept <- attr(terms(psiformula), "intercept")
  
  detfixed<-detterms[!grepl("|", detterms, fixed = TRUE)]
  psifixed<-psiterms[!grepl("|", psiterms, fixed = TRUE)]  
  
  detrandoms <- detterms[grepl("|", detterms, fixed = TRUE)]
  psirandoms <- psiterms[grepl("|", psiterms, fixed = TRUE)]
  
  detrandomsterms<-sapply(strsplit(detrandoms, "|", fixed = TRUE), 
                          function(x) return(gsub("[[:space:]]", "", x[1])))
  
  psirandomsterms<-sapply(strsplit(psirandoms, "|", fixed = TRUE), 
                          function(x) return(gsub("[[:space:]]", "", x[1])))
  
  
  detrandomsf<-sapply(strsplit(detrandoms, "|", fixed = TRUE), 
                      function(x) return(gsub("[[:space:]]", "", x[2])))
  
  psirandomsf<-sapply(strsplit(psirandoms, "|", fixed = TRUE), 
                      function(x) return(gsub("[[:space:]]", "", x[2])))
  
  
  
  
  detfformula<-formula(paste0(c(paste0("~",detintercept),detfixed),collapse = "+"))
  psifformula<-formula(paste0(c(paste0("~",psiintercept),psifixed),collapse = "+"))
  
  detrformula<-formula(paste0("~",detrandomsterms))
  psirformula<-formula(paste0("~",psirandomsterms))
  
  ### End of formula arrangements ####
  pred<-model.matrix(detfformula,data=bf)
  pred2<-model.matrix(psifformula,data=bf)
  
  prednames<-colnames(pred)
  prednames2<-colnames(pred2)
  
  npred<-ncol(pred)
  npred2<-ncol(pred2)
  
  rmatrix<-model.matrix(detrformula,data=bf)
  rmatrix2<-model.matrix(psirformula,data=bf)
  
  predrnames<-colnames(rmatrix)
  predrnames2<-colnames(rmatrix2)
  
  rnd<-match(predrnames,prednames)
  rnd2<-match(predrnames2,prednames2)
  
  if(length(rnd)==0){
    rnd<-npred+1
  }
  if(length(rnd2)==0){
    rnd2<-npred2+1
  }
  
  nrnd<-(1:npred)[!(1:npred)%in%rnd]
  nrnd2<-(1:npred2)[!(1:npred2)%in%rnd2]
  
  nrnd<-c(nrnd,npred+1)
  nrnd2<-c(nrnd2,npred2+1)
  
  #cbind((1:npred)%in%rnd,(1:npred)%in%nrnd)
  #cbind((1:npred2)%in%rnd2,(1:npred2)%in%nrnd2)
  bf$species<-as.factor(bf$species)
  bf$site<-as.factor(bf$site)

  data.jags<-list(nsp=nlevels(bf$species),
                  species=as.integer(bf$species),
                  nz=nrow(bf),
                  y=bf$y,
                  pred=pred,
                  npred=npred,
                  pred2=pred2,
                  npred2=npred2,
                  rnd=rnd,
                  rnd2=rnd2,
                  nrnd=nrnd,
                  nrnd2=nrnd2,
                  J=bf$J)
  
  
  parallelSeed<-parallel.seeds("base::BaseRNG", n.chains)
  
  
  inits<-function(chain){
    
    z=ifelse(data.jags$y>0,1,rbinom(length(data.jags$y),1,0.5))
    w_pre<-tapply(z,data.jags$species,max)
    w=ifelse(w_pre>0,1,rbinom(w_pre,1,0.5))
    
    list(
      rho0=rnorm(npred,0,2),
      psi0=rnorm(npred2,0,2),
      z=z,
      w=w,
      .RNG.name=parallelSeed[[chain]]$.RNG.name,
      .RNG.state=parallelSeed[[chain]]$.RNG.state
    )
  }


inits2<-lapply(1:3,inits)
  
  sim1rj<-run.jags("JAGS/MixedBinom.jags",
                   monitor = c("rho0","psi0","tau.rho","tau.psi","spX.rho","spX.psi","z","w"),
                   data = data.jags,
                   n.chains = n.chains,inits = inits2,burnin = burnin,adapt = n.adapt,sample = n.iter,method = method,thin = thin,silent.jags = TRUE,summarise = FALSE)
  
  
  sim1rj$vars<-list(detpred=prednames,psipred=prednames2,detrpred=predrnames,psirpred=predrnames2,spnames=levels(bf$species),sitenames=levels(bf$site),rnd=rnd,rnd2=rnd2,input=bf)
  
  class(sim1rj)<-c("runjags","occuMSOM")
  
  return(sim1rj)
  
}


coef.occuMSOM<-function(x){
  
  mcmc<-x$mcmc
  # Convert results to data frame
  sim1df<-as.data.frame(as.matrix(mcmc))
  
  # Get variable names
  names<-colnames(sim1df)
  
  ### Get results for random effects (standard deviations from fixed effects)
  SDdetection<-sim1df[,grepl("tau.rho",names),drop=FALSE]
  SDoccupancy<-sim1df[,grepl("tau.psi",names),drop=FALSE]
  
  colnames(SDdetection)[x$vars$rnd]<-x$vars$detrpred
  colnames(SDoccupancy)[x$vars$rnd2]<-x$vars$psirpred
  
  
  ### Get results for fixed effects
  Fdetection<-sim1df[,grepl("rho0",names),drop=FALSE]
  Foccupancy<-sim1df[,grepl("psi0",names),drop=FALSE]
  
  colnames(Fdetection)<-x$vars$detpred
  colnames(Foccupancy)<-x$vars$psipred
  
  ### Get the deviation from the fixed coefficient (random effec) for each individual species
  Rdetection<-sim1df[,grepl("spX.rho",names),drop=FALSE]
  Roccupancy<-sim1df[,grepl("spX.psi",names),drop=FALSE]
  
  dim(Roccupancy)
  dim(Rdetection)
  
  dim(Foccupancy)
  dim(Roccupancy)
  #SPoccupancy<-Foccupancy[,r2coeforder]+Roccupancy
  #SPdetection<-Fdetection[,spX.rhoorder]+Rdetection
  
  #### Calculate quantiles ####
  
  # Random effects (standard deviations)
  SDdetectionq<-apply(SDdetection,2,quantile,c(0.5,0.025,0.975))
  SDoccupancyq<-apply(SDoccupancy,2,quantile,c(0.5,0.025,0.975))
  
  Fdetectionq<-apply(Fdetection,2,quantile,c(0.5,0.025,0.975))
  Foccupancyq<-apply(Foccupancy,2,quantile,c(0.5,0.025,0.975))
  
  Rdetectionq<-apply(Rdetection,2,quantile,c(0.5,0.025,0.975))
  Roccupancyq<-apply(Roccupancy,2,quantile,c(0.5,0.025,0.975))
  
  SPdetectioncoefs<-matrix(Rdetectionq[1,],,ncol(SDdetection))
  SPoccupancycoefs<-matrix(Roccupancyq[1,],,ncol(SDoccupancy))
  
  colnames(SPdetectioncoefs)<-1:ncol(SDdetection)
  colnames(SPoccupancycoefs)<-1:ncol(SDoccupancy)
  
  colnames(SPdetectioncoefs)[x$vars$rnd]<-x$vars$detrpred
  colnames(SPoccupancycoefs)[x$vars$rnd2]<-x$vars$psirpred
  
  SPdetfinal<-t(t(SPdetectioncoefs[,1:ncol(Fdetectionq)])+Fdetectionq[1,])
  SPoccfinal<-t(t(SPoccupancycoefs[,1:ncol(Foccupancy)])+Foccupancyq[1,])
  
  rownames(SPdetfinal)<-x$vars$spnames
  rownames(SPoccfinal)<-x$vars$spnames
  
  colnames(SPdetfinal)<-x$vars$detpred
  colnames(SPoccfinal)<-x$vars$psipred
  
  SDdetection<-SDdetection[,x$vars$rnd,drop=FALSE]
  SDoccupancy<-SDoccupancy[,x$vars$rnd2,drop=FALSE]
  
  SDdetectionq<-SDdetectionq[,x$vars$rnd,drop=FALSE]
  SDoccupancyq<-SDoccupancyq[,x$vars$rnd2,drop=FALSE]
  
  (random=list(occupancy=SDoccupancyq[1,,drop=FALSE],detection=SDdetectionq[1,,drop=FALSE]))
  
  (fixed=list(occupancy=Foccupancyq[1,,drop=FALSE],detection=Fdetectionq[1,,drop=FALSE]))
  
  (spcoefs<-list(occupancy=SPoccfinal,detection=SPdetfinal))
  
  resu<-list(random=random,fixed=fixed,spcoefs=spcoefs)
  
  return(resu)
  
}

summary.occuMSOM<-function(x,returnZ=FALSE){
  
  require(runjags)
  
  mcmc<-x$mcmc
  # Convert results to data frame
  sim1df<-as.data.frame(as.matrix(mcmc))
  
  runjags.options(silent.runjags=TRUE)
  suppressWarnings(
    summ<-summary(x,vars = c("rho0","psi0","tau.rho","tau.psi","spX.rho","spX.psi"),silent.jags=TRUE)
  )
  runjags.options(silent.runjags=FALSE)

  
  
  # Get variable names
  names<-rownames(summ)
  
  SDdetection<-summ[grepl("tau.rho",names),,drop=FALSE][x$vars$rnd,,drop=FALSE]
  SDoccupancy<-summ[grepl("tau.psi",names),,drop=FALSE][x$vars$rnd2,,drop=FALSE]
  
  rownames(SDdetection)<-x$vars$detrpred
  rownames(SDoccupancy)<-x$vars$psirpred
  
  ### Get results for fixed effects
  Fdetection<-summ[grepl("rho0",names),,drop=FALSE]
  Foccupancy<-summ[grepl("psi0",names),,drop=FALSE]
  
  rownames(Fdetection)<-x$vars$detpred
  rownames(Foccupancy)<-x$vars$psipred
  
  ### Get the deviation from the fixed coefficient (random effect) for each individual species
  
  names.rho<-expand.grid(x$vars$spnames,c(x$vars$detpred,"fake"))
  names.psi<-expand.grid(x$vars$spnames,c(x$vars$psipred,"fake"))
  
  Rdetection<-summ[grepl("spX.rho",names),,drop=FALSE]
  Roccupancy<-summ[grepl("spX.psi",names),,drop=FALSE]
  
  rownames(Rdetection)<-paste(names.rho[,1],names.rho[,2],sep = ":")
  rownames(Roccupancy)<-paste(names.psi[,1],names.psi[,2],sep=":")
  
  Rdetection<-Rdetection[names.rho[,2]%in%x$vars$detrpred,]
  Roccupancy<-Roccupancy[names.psi[,2]%in%x$vars$psirpred,]
  

  #   pattern.rho<-paste0("spX.rho\\[.*,",paste0(x$vars$rnd,collapse = "|"),"\\]")
  #  pattern.rho<-paste0("?spX.rho\\[.*,",paste0(x$vars$rnd,collapse = "|"),"\\]")
  #  Rdetection<-summ[grepl(pattern.rho,names),,drop=FALSE]
  
  # z and w
  names2<-colnames(sim1df)
  
  z<-sim1df[,grepl("z\\[",names2),drop=FALSE]
  w<-sim1df[,grepl("w\\[",names2),drop=FALSE]
  
  zarray<-array(t(z),
                c(length(x$vars$sitenames),
                  length(x$vars$spnames),
                  nrow(z)),dimnames = list(x$vars$sitenames,x$vars$spnames,NULL))
  
  Sest<-apply(zarray,3,rowSums)
  Sestq<-apply(Sest,1,quantile,c(0.025,0.5,0.975))
  
  alpha.summary<-data.frame(t(Sestq),rowMeans(Sest),apply(Sest,1,sd))
  
  colnames(alpha.summary)<-c("Lower95","Median","Upper95", "Mean","SD")
  rownames(alpha.summary)<-x$vars$sitenames
  
  gamma<-rowSums(w)
  
  gamma.summary<-data.frame(t(quantile(gamma,c(0.025,0.5,0.975))),mean(gamma),sd(gamma))
  
  colnames(gamma.summary)<-c("Lower95","Median","Upper95", "Mean","SD")
  rownames(gamma.summary)<-"Regional diversity"
  

  (random=list(occupancy=SDoccupancy,detection=SDdetection))
  (fixed=list(occupancy=Foccupancy,detection=Fdetection))
  (spcoefs<-list(occupancy=Roccupancy,detection=Rdetection))
  
  resu<-list(random=random,fixed=fixed,spdeviations=spcoefs,alphaDiv=alpha.summary,gammaDiv=gamma.summary,z=if(returnZ){zarray})
  
  return(resu)
  
}

