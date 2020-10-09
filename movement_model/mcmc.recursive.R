#######################################################
### MCMC Algorithms for fitting two-stage procedure ###
#######################################################
# Author H McCaslin
# April 2020

# Inverse gamma density #
ldIG <-function(x,q,r){
  -(q+1)*log(x)-1/(r*x)-q*log(r)-lgamma(q)  
}

##### Stage one mcmc #####
s1.mcmc <- function(times,S,n.mcmc,q0,r0,mu.0,s2.0, phi.tune){
  
  # Set up #
  S <- as.matrix(S)   #telemetry observations
  t <- length(times)  #timestamps of observations
  dt <- diff(times)
  if(!is.numeric(times))  units(dt) <- "hours"  #when times, dt posixct objects
  dt <- as.numeric(dt)
  dim(dt) <- c((t-1),1) #accommodate non-constant dts
  
  # Prepare to save output #
  beta.save <- rep(0,n.mcmc)
  s2.save <- rep(0,n.mcmc)
  phi.save <- rep(0,n.mcmc)
  
  # Starting values #
  beta <- 0
  phi <- pi/4
  rot <- matrix(c(sin(phi), cos(phi)),nrow=2,ncol=1)
  accept <- 0
  
  for(k in 1:n.mcmc){
    
    # update s2 #
    tmp.r0 <- 1/(1/r0+0.5*sum(apply((S[-1,]-S[-t,]+dt%*%(beta*t(rot)))^2,1,sum)%*% dt^(-1)))
    tmp.q0 <- t-1+q0
    s2 <- 1/rgamma(1,shape=tmp.q0,scale=tmp.r0) 
    
    # update beta #
    tmp.var0 <- 1/( sum(dt)/s2 + 1/s2.0 ) 
    tmp.mn0 <- tmp.var0*(mu.0/s2.0 + sum(-(S[-1,]-S[-t,])%*%rot)/s2)
    beta <- rnorm(1,tmp.mn0,sqrt(tmp.var0))
    
    # update phi #
    phi.star <- rnorm(1,phi, phi.tune)

    if(phi.star>0 && phi.star<pi){
      rot.star <- matrix(c(sin(phi.star), cos(phi.star)),nrow=2,ncol=1)

      mh1 <- sum(dnorm(S[-1,],S[-t,]-dt%*%t(beta*rot.star),sqrt(s2*dt),log=T))
      mh2 <- sum(dnorm(S[-1,],S[-t,]-dt%*%t(beta*rot),sqrt(s2*dt),log=T))
      mh <- exp(mh1-mh2)
      if(mh>runif(1)){
        phi <- phi.star
        rot <- rot.star
        accept <- accept+1
      }
    }

    # save samples #
    beta.save[k] <- beta
    s2.save[k] <- s2
    phi.save[k] <- phi
  }
  # write output #
  list(beta.save=beta.save,s2.save=s2.save,phi.save=phi.save,n.mcmc=n.mcmc, phi.accept=accept/n.mcmc)
}


##### Stage two mcmc #####
s2.mcmc <- function(out.arr, q0,r0, mu.0,s2.0, mu.1,s2.1, q1,r1, q2,r2, mu.2,s2.2, n.mcmc){
  
  z.arr <- out.arr  #out.arr from stage one; dim JxKx3 (beta,s2,phi)
  J <-dim(z.arr)[1]
  Z <- matrix(nrow=J, ncol=3) 
  
  # Prepare to save output #
  beta.save <- matrix(0,nrow=J,ncol=n.mcmc)
  s2.save <- matrix(0,nrow=J,ncol=n.mcmc)
  log.sig.save <- matrix(0,nrow=J,ncol=n.mcmc)
  phi.save <- matrix(0,nrow=J, ncol=n.mcmc)
  mu.beta.save <- numeric(n.mcmc)
  s2.beta.save<- numeric(n.mcmc)
  mu.sigma.save <-numeric(n.mcmc)
  s2.sigma.save <-numeric(n.mcmc)
  accept <- numeric(J)
  
  # Starting values #
  post.length <- dim(z.arr)[2]
  
  beta <- rowMeans(z.arr[,,1])  #use first-stage means as initial values
  s2 <- rowMeans(z.arr[,,2])
  log.sig <- log(sqrt(s2))
  phi <- rowMeans(z.arr[,,3])
  Z <- cbind(beta, s2, phi) 
  
  mu.beta <- 1
  s2.beta <- 0.2
  mu.sigma <- -1.5
  s2.sigma <- 0.2
  
  for(k in 1:n.mcmc){
    
    # jointly update beta, s2, phi #
    for(j in 1:J){
      idx <- sample(1:post.length,1) 
      z.star <- z.arr[j,idx,]  #draws from first stage posterior
      log.sig.star <- log(sqrt(z.star[2]))
      
      mh1 <- dnorm(z.star[1],mu.beta,sqrt(s2.beta),log=T) + 
        dnorm(log.sig.star,mu.sigma,sqrt(s2.sigma),log=T) + 
        dnorm(Z[j,1],mu.0,sqrt(s2.0),log=T) + 
        ldIG(Z[j,2],q0,r0) + log(Z[j,2])
      mh2 <- dnorm(Z[j,1],mu.beta,sqrt(s2.beta),log=T) + 
        dnorm(log.sig[j],mu.sigma,sqrt(s2.sigma),log=T) +
        dnorm(z.star[1],mu.0,sqrt(s2.0),log=T) + 
        ldIG(z.star[2],q0,r0) + log(z.star[2])
      
      mh <- exp(mh1-mh2)
      if(mh > runif(1)){
        Z[j,] <- z.star
        accept[j] <- accept[j] + 1
      }
    }
    
    beta <- Z[,1]
    s2 <- Z[,2]
    log.sig <- log(sqrt(s2))
    phi <- Z[,3]
    
    # update s2.beta #
    tmp.r1 <- 1/(sum((beta-mu.beta)^2)/2 + 1/r1)
    tmp.q1 <- J/2 + q1
    s2.beta <- 1/rgamma(1,shape=tmp.q1,scale=tmp.r1)
    
    # update mu.beta #
    tmp.var <- 1/(J/s2.beta + 1/s2.1)
    tmp.mn <- tmp.var*(sum(beta)/s2.beta + mu.1/s2.1)
    mu.beta <- rnorm(1,tmp.mn,sqrt(tmp.var))
    
    # update s2.sigma #
    tmp.r2 <- 1/(sum((log.sig - mu.sigma)^2)/2 + 1/r2) 
    tmp.q2 <- J/2 + q2
    s2.sigma <- 1/rgamma(1,shape=tmp.q2,scale=tmp.r2)
    
    # update mu.sigma #
    tmp.v.sig <- 1/(J/s2.sigma + 1/s2.2)
    tmp.m.sig <- tmp.v.sig*(sum(log.sig)/s2.sigma + mu.2/s2.2)
    mu.sigma <- rnorm(1,tmp.m.sig, sqrt(tmp.v.sig))
    
    # save samples #
    beta.save[,k] <- beta
    s2.save[,k] <- s2
    log.sig.save[,k] <- log.sig
    phi.save[,k] <- phi
    mu.beta.save[k] <- mu.beta
    s2.beta.save[k] <- s2.beta
    mu.sigma.save[k] <- mu.sigma
    s2.sigma.save[k] <- s2.sigma
  }
  
  # write output #
  list(beta.save=beta.save,s2.save=s2.save,log.sig.save=log.sig.save,phi.save=phi.save,mu.beta.save=mu.beta.save,s2.beta.save=s2.beta.save,mu.sigma.save=mu.sigma.save,s2.sigma.save=s2.sigma.save, accept=accept/n.mcmc, n.mcmc)
  
}


