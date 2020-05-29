mod1.func<-function(Y, W){
mod1 <- "
model{

for(i in 1:n.times){

Y[i] ~ dpois(lambda[i])

log(lambda[i]) = phi.y0 + phi.y[i]

#Observation error on X
W[i,1] ~ dnorm(X1[i], prec.w1)  #Replicate 1, Target 1
W[i,3] ~ dnorm(X1[i], prec.w1)  #Replicate 2, Target 1
W[i,2] ~ dnorm(X1[i], prec.w1)  #Replicate 1, Target 2
W[i,4] ~ dnorm(X1[i], prec.w1)  #Replicate 2, Target 2

X1[i] =  phi.x0 + phi.x[i]  

}

phi.y0 ~ dnorm(0.00, 0.0001)
phi.x0 ~ dnorm(0.00, 0.0001)

#precision on observation of viral RNA
prec.w1 <- 1/sdW^2
sdW ~ dunif(0,100)

#AR1 for X and Y
phi.x[1] ~ dnorm(0.00, tau.rw.x)
phi.y[1] ~ dnorm(0.00, tau.rw.y)

for(g in 2:n.times){
  phi.x[g]~ dnorm(rho.x*phi.x[g-1], tau.rw.x)
  phi.y[g]~ dnorm(rho.y*phi.y[g-1], tau.rw.y)
}

rho.x~dunif(0,1)
rho.y~dunif(0,1)
tau.rw.x <- 1/sd.rw.x^2
tau.rw.y <- 1/sd.rw.y^2

sd.rw.x ~ dunif(0, 100)
sd.rw.y ~ dunif(0, 100)

deriv.x[1] <-0
deriv.y[1] <-0

for(i in 2:n.times){
  deriv.x[i] <- phi.x[i] - phi.x[i-1]
  deriv.y[i] <- phi.y[i] - phi.y[i-1]

}


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

model_spec<-textConnection(mod1)
model_jags<-jags.model(model_spec, 
                       inits=list(inits1,inits2, inits3),
                       data=list('Y' = Y,
                                 'W' = W,
                                 'n.times'=(length(Y))),
                       n.adapt=10000, 
                       n.chains=3)

params<-c('phi.y','phi.x', 'deriv.x','deriv.y','X1','lambda')

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
ci<-t(hdi(posterior_samples.all, credMass = 0.95))
row.names(ci)<-sample.labs
names(post_means)<-sample.labs
post.comb <- cbind.data.frame(post_means, ci)

rw.x.index <- grep('phi.x', sample.labs)
deriv.x.index <- grep('deriv.x',sample.labs)
rw.y.index <- grep('phi.y', sample.labs)
deriv.y.index <- grep('deriv.y', sample.labs)
x1.index <- grep('X1', sample.labs)
lambda.index <- grep('lambda', sample.labs)

rw.x <- post.comb[rw.x.index,]
rw.y <- post.comb[rw.y.index,]
X1 <- post.comb[x1.index,]
lambda <- post.comb[lambda.index,]

deriv.x <- post.comb[deriv.x.index,]
deriv.y <- post.comb[deriv.y.index,]
outlist.mod1 <- list('rw.x'=rw.x,'rw.y'=rw.y,'deriv.x'=deriv.x,'deriv.y'=deriv.y,'X1'=X1,'lambda'=lambda)
return(outlist.mod1)
}

