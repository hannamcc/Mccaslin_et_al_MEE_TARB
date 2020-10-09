#####################################################
### MCMC Algorithm for fitting hierarchical model ###
#####################################################
# Author H McCaslin
# April 2020

# Inverse gamma density #
ldIG <-function(x,q,r){
    -(q+1)*log(x)-1/(r*x)-q*log(r)-lgamma(q)  
}

##### Single stage MCMC #####
hier.mcmc <- function(y.list, mu.1,s2.1, q1,r1, q2,r2, mu.2,s2.2, n.mcmc, logsig.tune, phi.tune){
    
    # Number of individuals #
    J <- length(y.list)  #each item in list is matrix of one individual: time, lat, lon
    
    # Prepare to save output #
    beta.save <- matrix(0,nrow=J,ncol=n.mcmc)
    log.sig.save <- matrix(0,nrow=J,ncol=n.mcmc)
    s2.save <- matrix(0,nrow=J,ncol=n.mcmc)
    phi.save <- matrix(0,nrow=J,ncol=n.mcmc)
    mu.beta.save <- numeric(n.mcmc)
    s2.beta.save<- numeric(n.mcmc)
    mu.sigma.save <-numeric(n.mcmc)
    s2.sigma.save <-numeric(n.mcmc)
    
    accept.sig <- numeric(J)
    accept.phi <- numeric(J)
    
    # Starting values #
    beta <- rep(1,times=J)
    log.sig <- rep(-1,times=J)
    s2 <- exp(log.sig)^2
    phi <- rep(pi/4,times=J)
    rot <- matrix(c(sin(phi),cos(phi)),2,J) 
    
    mu.beta <- 1
    s2.beta <- 0.2
    mu.sigma <- -1.5
    s2.sigma <- 0.2
    
    for(k in 1:n.mcmc){
        
        ## update indiv-level pars ##
        for(j in 1:J){
            S <- as.matrix(y.list[[j]][,2:3])
            t <- dim(S)[1]
            
            dt <- diff(y.list[[j]][,1]) 
            if(!is.numeric(y.list[[j]][,1])) units(dt) <- "hours"    #dt posixct object
            dt <- as.numeric(dt)
            dim(dt) <- c((t-1),1)   #accommodate non-constant dts
            
            # update beta #
            tmp.var0 <- 1/( sum(dt)/s2[j] + 1/s2.beta )
            tmp.mn0 <- tmp.var0*(mu.beta/s2.beta + sum(-(S[-1,]-S[-t,])%*%rot[,j])/s2[j])
            beta[j] <- rnorm(1,tmp.mn0,sqrt(tmp.var0))

            # update log(sigma) #
            log.sig.star <- rnorm(1, log.sig[j], logsig.tune[j])
            
            mh <- sum(dnorm(S[-1,],S[-t,]-dt%*%t(beta[j]*rot[,j]),sqrt(dt)*exp(log.sig.star),log=T)) + 
                dnorm(log.sig.star,mu.sigma,sqrt(s2.sigma),log=T) -
                sum(dnorm(S[-1,],S[-t,]-dt%*%t(beta[j]*rot[,j]),sqrt(dt)*exp(log.sig[j]),log=T)) -
                dnorm(log.sig[j],mu.sigma,sqrt(s2.sigma),log=T)

            mh <- exp(mh)
            
            if(mh > runif(1)){
                log.sig[j] <- log.sig.star
                s2[j] <- exp(log.sig[j])^2
                accept.sig[j] <- accept.sig[j] + 1
            }
            
            # update phi #
            phi.star <- rnorm(1,phi[j], phi.tune)

            if(phi.star>0 && phi.star<pi){
                rot.star <- matrix(c(sin(phi.star), cos(phi.star)),nrow=2,ncol=1)

                mh1 <- sum(dnorm(S[-1,],S[-t,]-dt%*%t(beta[j]*rot.star),sqrt(s2[j]*dt),log=T))
                mh2 <- sum(dnorm(S[-1,],S[-t,]-dt%*%t(beta[j]*rot[,j]),sqrt(s2[j]*dt),log=T))
                mh <- exp(mh1-mh2)
                if(mh>runif(1)){
                    phi[j] <- phi.star
                    rot[,j] <- rot.star
                    accept.phi[j] <- accept.phi[j]+1
                }
            }
        }
        
        ## population-level parameters ##
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
        log.sig.save[,k] <- log.sig
        s2.save[,k] <- s2
        phi.save[,k] <- phi
        mu.beta.save[k] <- mu.beta
        s2.beta.save[k] <- s2.beta
        mu.sigma.save[k] <- mu.sigma
        s2.sigma.save[k] <- s2.sigma
    }
    # write output #
    list(beta.save=beta.save,log.sig.save=log.sig.save,s2.save=s2.save,phi.save=phi.save, mu.beta.save=mu.beta.save,s2.beta.save=s2.beta.save,mu.sigma.save=mu.sigma.save,s2.sigma.save=s2.sigma.save, accept.sig=accept.sig/n.mcmc,accept.phi=accept.phi/n.mcmc, n.mcmc)
    
}

