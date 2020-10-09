library(scales)
library(parallel)
numCores <- detectCores()

#### Create simulation ####
# function to simulate #
simulate <- function(beta,s2,phi,dt,nsim){
     s <- matrix(ncol=2,nrow=nsim)
     times <- numeric(nsim)
     s[1,] <- runif(2, c(0.7,0.9), c(0.85,0.9))
     times[1] <- dt
     
     for(t in 2:nsim){
          times[t] <- dt*t
          s[t,] <- rnorm(2, s[t-1,]-c(sin(phi),cos(phi))*beta*dt, sqrt(s2*dt)) 
     }
     sim=cbind(times, x=s[,1], y=s[,2])
}

## Simulate tracks ##
# parameters #
set.seed(100)
beta.true <- rnorm(20, 1.8, sqrt(0.1))  #draw betas and log(sigma)'s from 'population' distributions
log.sig.true <- rnorm(20,-2.2,sqrt(0.01)) 
s2.true <- exp(log.sig.true)^2

phi <- runif(20, pi/6, pi/3)
rot <- c(sin(phi), cos(phi))

# get simulations #
sims <- list(20)
set.seed(100)
for(j in 1:20){
     sims[[j]] <- simulate(beta.true[j],s2.true[j],phi[j],.005,100)
}

## Plot sims ##
plot(sims[[1]][,2:3], xlim=c(0,1), ylim=c(0,1),lwd=1.2,col=alpha(1, 0.5), type="l",
     xlab="X", ylab="Y",bty="n")
for(j in 2:20){
     lines(sims[[j]][,2:3], lwd=1.2,col=alpha(j, 0.7))
}

#############################################################
## Fit two stage model ##
source("mcmc.recursive.R")
J=20
n.mcmc=20000 
n.burn=4000
out.arr <- array(dim=c(J,n.mcmc-n.burn,3))

# Stage 1 #
fit.s1.sims <- function(sims){
        fit.j <- s1.mcmc(times=sims[,1],S=sims[,2:3],n.mcmc=n.mcmc,q0=0.01,r0=100,mu.0=0,s2.0=10, phi.tune=0.2)
        beta.vec <- fit.j$beta.save[-(1:n.burn)]
        s2.vec <- fit.j$s2.save[-(1:n.burn)]
        phi.vec <- fit.j$phi.save[-(1:n.burn)]
        list(pars=rbind(beta.vec, s2.vec, phi.vec), accept=fit.j$phi.accept)
}

stage1 <- mclapply(sims, fit.s1.sims, mc.cores=numCores)

# organize output #
for(j in 1:J){
        out.arr[j,,1] <- stage1[[j]]$pars[1,]
        out.arr[j,,2] <- stage1[[j]]$pars[2,]
        out.arr[j,,3] <- stage1[[j]]$pars[3,]
}

# Stage 2 #
stage2 <- s2.mcmc(out.arr, q0=0.01,r0=100, mu.0=0,s2.0=10, mu.1=0,s2.1=10, q1=0.001,r1=1000, q2=0.001,r2=1000, mu.2=0,s2.2=100, n.mcmc=n.mcmc)

## Fit hierarchical algorithm ##
source("mcmc.hierarchical.R")
fit.hier <- hier.mcmc(y.list=sims, mu.1=0,s2.1=10, q1=0.001,r1=1000, q2=0.001,r2=1000, mu.2=0,s2.2=100, n.mcmc=n.mcmc,logsig.tune = rep(0.12,J), phi.tune = rep(0.2,J)) 

###########################
#### Visualize results ####
###########################
three.cols <- c("#000066", "#009966", "#99CC00")
par(mfrow=c(3,1), mar=c(3.5, 5, 1.5, 2.1),mgp=c(2.2, 1, 0),xpd=NA)

## beta ##
# hierarchical #
matplot(1:J,rowMeans(fit.hier$beta.save[,-(1:n.burn)]), type="n", xlim=c(0.75,23), ylim=c(0.5,3.5), xlab="Individual", ylab=expression(beta[j]),xaxt='n') 
axis(side=1, at=c(1:20,22.3), labels=c(1:20, expression(mu[beta])), tck=-0.03)
segments(1:J,apply(fit.hier$beta.save[,-(1:n.burn)],1,quantile,0.025),1:J,apply(fit.hier$beta.save[,-(1:n.burn)],1,quantile,0.975),lwd=2, col=three.cols[1])
points(1:J, rowMeans(fit.hier$beta.save[,-(1:n.burn)]), pch=1, col=three.cols[1])
# stage 1 #
segments(1:J+0.15,apply(out.arr[,,1],1,quantile,0.025),1:J+0.15,apply(out.arr[,,1],1,quantile,0.975),lwd=2, col=three.cols[2])
points(1:J+0.15, rowMeans(out.arr[,,1]), pch=1, col=three.cols[2])
# stage 2 #
segments(1:J+0.3,apply(stage2$beta.save[,-(1:n.burn)],1,quantile,0.025),1:J+0.3,apply(stage2$beta.save[,-(1:n.burn)],1,quantile,0.975),lwd=2, col=three.cols[3])
points(1:J+0.3, rowMeans(stage2$beta.save[,-(1:n.burn)]), pch=1, col=three.cols[3])
# population-level #
segments(21,0.5,21,3.5, lty="dashed",col="gray")
points(22, mean(fit.hier$mu.beta.save[-(1:n.burn)]), pch=1, col=three.cols[1])
segments(22, quantile(fit.hier$mu.beta.save[-(1:n.burn)], 0.025), 22, quantile(fit.hier$mu.beta.save[-(1:n.burn)], 0.975), lwd=2, col=three.cols[1])
points(22.3, mean(stage2$mu.beta.save[-(1:n.burn)]), pch=1, col=three.cols[3])
segments(22.3, quantile(stage2$mu.beta.save[-(1:n.burn)], 0.025), 22.3, quantile(stage2$mu.beta.save[-(1:n.burn)], 0.975), lwd=2, col=three.cols[3])
legend("topright", c("Full", "Stage 1", "Stage 2"),col=three.cols,lty=1,lwd=2, bty="n", seg.len=1.2,y.intersp=.5) 

## log(sigma) ##
# hierarc hcal #
matplot(1:J,rowMeans(fit.hier$log.sig.save[,-(1:n.burn)]), type="n",xlim=c(0.75,23),  ylim=c(-2.6,-1.6),xlab="Individual", ylab=expression(paste("log(",sigma[j],")")),xaxt='n')
axis(side=1, at=c(1:20,22.3), labels=c(1:20, expression(mu[sigma])), tck=-0.03)
segments(1:J,apply(fit.hier$log.sig.save[,-(1:n.burn)],1,quantile,0.025),1:J,apply(fit.hier$log.sig.save[,-(1:n.burn)],1,quantile,0.975),lwd=2, col=three.cols[1])
points(1:J, rowMeans(fit.hier$log.sig.save[,-(1:n.burn)]), pch=1,  col=three.cols[1])
# stage 1 #
segments(1:J+0.15,apply(log(sqrt(out.arr[,,2])),1,quantile,0.025),1:J+0.15,apply(log(sqrt(out.arr[,,2])),1,quantile,0.975),lwd=2,  col=three.cols[2])
points(1:J+0.15, rowMeans(log(sqrt(out.arr[,,2]))), pch=1, col=three.cols[2])
# stage 2 #
segments(1:J+0.3,apply(log(sqrt(stage2$s2.save))[,-(1:n.burn)],1,quantile,0.025),1:J+0.3,apply(log(sqrt(stage2$s2.save))[,-(1:n.burn)],1,quantile,0.975),lwd=2,  col=three.cols[3])
points(1:J+0.3, rowMeans(log(sqrt(stage2$s2.save))[,-(1:n.burn)]), pch=1,  col=three.cols[3])
# population-level #
segments(21,-2.6,21,-1.6, lty="dashed",col="gray")
points(22, mean(fit.hier$mu.sigma.save[-(1:n.burn)]), pch=1, col=three.cols[1])
segments(22, quantile(fit.hier$mu.sigma.save[-(1:n.burn)], 0.025), 22, quantile(fit.hier$mu.sigma.save[-(1:n.burn)], 0.975), lwd=2, col=three.cols[1])
points(22.3, mean(stage2$mu.sigma.save[-(1:n.burn)]), pch=1, col=three.cols[3])
segments(22.3, quantile(stage2$mu.sigma.save[-(1:n.burn)], 0.025), 22.3, quantile(stage2$mu.sigma.save[-(1:n.burn)], 0.975), lwd=2, col=three.cols[3])
legend("topright", c("Full", "Stage 1", "Stage 2"),col=three.cols,lty=1,lwd=2, bty="n", seg.len=1.2,y.intersp=.5)

## phi ##
# hierarchical #
matplot(1:J,rowMeans(fit.hier$phi.save[,-(1:n.burn)]), type="n", xlim=c(0.75,23),ylim=c(0.25,1.75), xlab="Individual", ylab=expression(phi[j]),xaxt='n')
axis(side=1, at=c(1:20), labels=c(1:20), tck=-0.03)
segments(1:J,apply(fit.hier$phi.save[,-(1:n.burn)],1,quantile,0.025),1:J,apply(fit.hier$phi.save[,-(1:n.burn)],1,quantile,0.975),lwd=2, col=three.cols[1])
points(1:J, rowMeans(fit.hier$phi.save[,-(1:n.burn)]), pch=1, col=three.cols[1])
# stage 1 #
segments(1:J+0.15,apply(out.arr[,,3],1,quantile,0.025),1:J+0.15,apply(out.arr[,,3],1,quantile,0.975),lwd=2, col=three.cols[2])
points(1:J+0.15, rowMeans(out.arr[,,3]), pch=1, col=three.cols[2])
# stage 2 #
segments(1:J+0.3,apply(stage2$phi.save[,-(1:n.burn)],1,quantile,0.025),1:J+0.3,apply(stage2$phi.save[,-(1:n.burn)],1,quantile,0.975),lwd=2, col=three.cols[3])
points(1:J+0.3, rowMeans(stage2$phi.save[,-(1:n.burn)]), pch=1, col=three.cols[3])
legend("topright", c("Full", "Stage 1", "Stage 2"),col=three.cols,lty=1,lwd=2, bty="n", seg.len=1.2,y.intersp=.5)

