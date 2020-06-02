mod2.indiv.lag.func<-function(W, Y, log.offset, lag.n=0) {
  mod2 <- "
  model{

for(i in 1:(n.times)){

#Observation error on X
W[i,1] ~ dnorm(X[i], prec.w1)  #Replicate 1, Target 1
W[i,3] ~ dnorm(X[i], prec.w1)  #Replicate 2, Target 1
W[i,2] ~ dnorm(X[i], prec.w1)  #Replicate 1, Target 2
W[i,4] ~ dnorm(X[i], prec.w1)  #Replicate 2, Target 2

X[i] = phi.x0 + phi.x[i]   #X[i] is a RW1

}

phi.x0 ~ dnorm(0.00, 0.0001)

for(i in (7+1):(n.times)){

Y[i] ~ dpois(lambda[i])

log(lambda[i]) = 
(beta0 + 
 beta1*X[i-lag.n] + 
 phi.y[i] + #phi.y is an AR1
 log.offset[i]
 ) 
}

beta0 ~ dnorm(0.00, 0.0001)
 
#Distributed lag parameters have a RW structure
beta1 ~ dnorm(0.00, 0.0001)


#precision on observation of viral RNA
prec.w1 <- 1/sdW^2
sdW ~dunif(0,100)


#RW for X, AR(1) y
phi.x[1] ~ dnorm(0.00, tau.rw.x)
phi.y[1] ~ dnorm(0.00, tau.rw.y)

for(g in 2:(n.times)){
  phi.x[g] ~ dnorm(phi.x[g-1] , tau.rw.x)
  phi.y[g] ~ dnorm(rho.y*phi.y[g-1] , tau.rw.y)
}

rho.y~dunif(0,1)

tau.rw.x <- 1/sd.rw.x^2
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
                                 'lag.n'=lag.n,
                                 'log.offset'=log.offset,
                                 'n.times'=(length(Y))),
                       n.adapt=10000, 
                       n.chains=3)

params<-c('lambda','phi.x', 'beta1','X','delta.dow')

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
ci<-t(hdi(posterior_samples.all, credMass = 0.9875))
row.names(ci)<-sample.labs
names(post_means)<-sample.labs
post.comb <- cbind.data.frame(post_means, ci)

rw.x.index <- grep('phi.x', sample.labs)
lambda.index <- grep('lambda', sample.labs)
beta1.index <- grep('beta1', sample.labs, fixed=T)
x.index <- grep('X', sample.labs)
delta.dow.index <- grep('delta.dow', sample.labs)



rw.x <- post.comb[rw.x.index,]
lambda <- post.comb[lambda.index,]
beta1 <- post.comb[beta1.index,]
X <- post.comb[x.index,]
delta.dow <- post.comb[delta.dow.index,]

outlist.mod1 <- list('rw.x'=rw.x,'beta1'=beta1,'delta.dow'=delta.dow,'lambda'=lambda,'X'=X)
return(outlist.mod1)
}

