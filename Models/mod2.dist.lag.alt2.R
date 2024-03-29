mod2.func<-function(W, Y, log.offset, nlags=8, nleads=1) {
  mod2 <- "
  model{

for(i in 1:(n.times)){

#Observation error on X
W[i,1] ~ dnorm(X[i], prec.w1)  #Replicate 1, Target 1
W[i,3] ~ dnorm(X[i], prec.w1)  #Replicate 2, Target 1
W[i,2] ~ dnorm(X[i], prec.w1)  #Replicate 1, Target 2
W[i,4] ~ dnorm(X[i], prec.w1)  #Replicate 2, Target 2

#X[i] = phi.x0 + phi.x[i]  #X[i] is a RW1
X[i] ~ dnorm(0.00, 0.0001)


}

phi.x0 ~ dnorm(0.00, 0.0001)

for(i in (nlags+1):(n.times-2)){

Y[i] ~ dpois(lambda[i])

log(lambda[i]) = 
(beta0 + 
 beta1[1]*X[i+1] + 
 beta1[2]*X[i] + 
 beta1[3]*X[i-1] +
 beta1[4]*X[i-2] +
 beta1[5]*X[i-3] +
 beta1[6]*X[i-4] +
 beta1[7]*X[i-5] +
 beta1[8]*X[i-6] +
 beta1[9]*X[i-7] +
 beta1[10]*X[i-8] +
 #beta1[11]*X[i-9] +

 phi.y[i] + #phi.y is an AR1
 log.offset[i]
 ) 
}

beta0 ~ dnorm(0.00, 0.0001)
 
#Distributed lag parameters have a RW structure
beta1[1] ~ dnorm(0.00, 0.0001)
beta1.cum[1] <- beta1[1]

for(k in 2:(nlags + nleads+1)){
  beta1[k] ~ dnorm(beta1[k-1], prec.beta)
  beta1.cum[k] <- sum(beta1[1:k])
}


#precision on observation of viral RNA
prec.w1 <- 1/sdW^2
sdW ~dunif(0,100)

prec.beta <- 1/sd.beta^2
sd.beta ~ dunif(0,100)

#RW for X, AR(1) y
#phi.x[1] ~ dnorm(0.00, tau.rw.x)
phi.y[1] ~ dnorm(0.00, tau.rw.y)

for(g in 2:(n.times)){
#  phi.x[g] ~ dnorm(phi.x[g-1] , tau.rw.x)
  phi.y[g] ~ dnorm(rho.y*phi.y[g-1] , tau.rw.y)
}

rho.y~dunif(0,1)

#tau.rw.x <- 1/sd.rw.x^2
tau.rw.y <- 1/sd.rw.y^2

sd.rw.x ~ dunif(0,100)
sd.rw.y ~ dunif(0,100)


}
"


##############################################################
#Model Fitting
##############################################################
inits1=list(".RNG.seed"=c(123), ".RNG.name"='base::Wichmann-Hill')
inits2=list(".RNG.seed"=c(456), ".RNG.name"='base::Wichmann-Hill')
inits3=list(".RNG.seed"=c(789), ".RNG.name"='base::Wichmann-Hill')


##############################################
#Model Organization
##############################################

model_spec<-textConnection(mod2)
model_jags<-jags.model(model_spec, 
                       inits=list(inits1,inits2, inits3),
                       data=list('Y' = Y,
                                 'W' = W,
                                 'nlags'=nlags,
                                 'nleads'=nleads,
                                 'log.offset'=log.offset,
                                 'n.times'=(length(Y))),
                       n.adapt=20000, 
                       n.chains=3)

params<-c('lambda','phi.x', 'beta1', 'beta1.cum','X', 'phi.y')

##############################################
#Posterior Sampling
##############################################
posterior_samples.mod1<-coda.samples(model_jags, 
                                     params, 
                                     n.iter=10000)

#Process output
posterior_samples.all<-do.call(rbind,posterior_samples.mod1)
post_means<-apply(posterior_samples.all, 2, median)
sample.labs<-names(post_means)
ci<-t(hdi(posterior_samples.all, credMass = 0.90))
row.names(ci)<-sample.labs
names(post_means)<-sample.labs
post.comb <- cbind.data.frame(post_means, ci)

#rw.x.index <- grep('phi.x', sample.labs)
lambda.index <- grep('lambda', sample.labs)
beta1.index <- grep('beta1[', sample.labs, fixed=T)
beta1.cum.index <- grep('beta1.cum', sample.labs)
x.index <- grep('X', sample.labs)


#rw.x <- post.comb[rw.x.index,]
lambda <- post.comb[lambda.index,]
beta1 <- post.comb[beta1.index,]
beta1.cum <- post.comb[beta1.cum.index,]
X <- post.comb[x.index,]

outlist.mod1 <- list('beta1'=beta1,'beta1.cum'=beta1.cum,'lambda'=lambda,'X'=X,'posterior_samples.all'=posterior_samples.all)
return(outlist.mod1)
}

