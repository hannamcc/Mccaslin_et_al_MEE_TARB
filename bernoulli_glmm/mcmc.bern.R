######################################################
### Stage 2 MCMC algorithm for bernoulli glmm with ###
###                  logit link                    ###
######################################################

# Logit & inverse logit functions #
logit <- function(phi){
     log(phi/(1-phi))
}

logit.inv <- function(phi){
        exp(phi)/(1+exp(phi))
}

# Stage 2 algorithm #
bern.s2 <- function(mat, alpha, beta, n.mcmc){
        J <- dim(mat)[1]
        s1_length <- dim(mat)[2]
        
        # prepare to save output #
        mu.save <- numeric(n.mcmc)
        s2.save <- numeric(n.mcmc)
        p.save <- matrix(nrow=J,ncol=n.mcmc)
        logit.p.save <- matrix(nrow=J,ncol=n.mcmc)
        
        # hyperpriors and starting values #
        mu_0 <- 0.7
        s2_0 <- 1
        q <- .01
        r <- 100
        
        mu <- 0.7
        s2 <- 1
        p <- rowMeans(mat) 
        logit.p <- logit(p)
        
        for(k in 1:n.mcmc){
                
                ## update logit(p) for each j##
                for(j in 1:J){
                     idx <- sample(1:s1_length,1,replace=T) 
                     p.star <- mat[j,idx]  
                     logit.p.star <- logit(p.star)
                     
                     mh1 <- dnorm(logit.p.star,mu,sqrt(s2),log=T) + dbeta(p[j],alpha,beta,log=T) + 
                           log(p[j]-p[j]^2)
                     mh2 <- dnorm(logit.p[j],mu,sqrt(s2),log=T) + dbeta(p.star,alpha,beta,log=T) + 
                           log(p.star-p.star^2)
                     
                     mh <- exp(mh1-mh2) 
                     if(mh > runif(1)){
                        p[j] <- p.star
                        logit.p[j] <- logit.p.star
                     }
                }
                
                ## sigma update ##
                tmp_r <- 1/(sum((logit.p-mu)^2)/2 + 1/r) 
                tmp_q <- J/2 + q
                s2 <- 1/rgamma(1,shape=tmp_q,scale=tmp_r)
                
                ## mu update ##
                tmp_s2 <- 1/(J/s2 + 1/s2_0)
                tmp_mu <- tmp_s2*(sum(logit.p)/s2 + mu_0/s2_0)
                mu <- rnorm(1,tmp_mu, sqrt(tmp_s2))
                
                ## save samples ##
                mu.save[k] <- mu
                s2.save[k] <- s2
                p.save[,k] <- p
                logit.p.save[,k] <- logit.p
        }
        # write output #
        list(mu.save=mu.save, s2.save=s2.save, p.save=p.save, logit.p.save=logit.p.save)
}

######################################################
### MCMC algorithm to fit bernoulli GLMM with logit 
### link conventionally, with MH updates for logit(p)
######################################################
bern.full <- function(y.list, tune, n.mcmc){
        
        J=length(y.list)
        
        # prepare to save output #
        mu.save <- numeric(n.mcmc)
        s2.save <- numeric(n.mcmc)
        logit.p.save <- matrix(nrow=J,ncol=n.mcmc)
        accept <- numeric(J)
        
        # hyperpriors and starting values # 
        mu_0 <- 0.7
        s2_0 <- 1
        q <- .01
        r <- 100
        
        mu <- 0.7
        s2 <- 1
        logit.p <- logit(rep(0.5,J))
        
        for(k in 1:n.mcmc){
                
                ## update logit(p) for each j##
                for(j in 1:J){
                        y <- y.list[[j]][,3]
                        
                        logit.p.star <- rnorm(1, logit.p[j], tune[j])
                        
                        mh1 <- sum(dbinom(y,1,logit.inv(logit.p.star),log=T)) + dnorm(logit.p.star,mu,sqrt(s2),log=T)
                        mh2 <- sum(dbinom(y,1,logit.inv(logit.p[j]),log=T)) + dnorm(logit.p[j],mu,sqrt(s2),log=T)       
                        mh <- exp(mh1-mh2) 
                        if(mh > runif(1)){
                                logit.p[j] <- logit.p.star
                                accept[j] <- accept[j] + 1
                        }
                    
                }
                
                ## sigma update ##
                tmp_r <- 1/(sum((logit.p-mu)^2)/2 + 1/r) 
                tmp_q <- J/2 + q
                s2 <- 1/rgamma(1,shape=tmp_q,scale=tmp_r)
                
                ## mu update ##
                tmp_s2 <- 1/(J/s2 + 1/s2_0)
                tmp_mu <- tmp_s2*(sum(logit.p)/s2 + mu_0/s2_0)
                mu <- rnorm(1,tmp_mu, sqrt(tmp_s2))
                
                ## save samples ##
                mu.save[k] <- mu
                s2.save[k] <- s2
                logit.p.save[,k] <- logit.p
        }
        # write output #
        list(mu.save=mu.save, s2.save=s2.save, logit.p.save=logit.p.save, accept=accept/n.mcmc)
}

