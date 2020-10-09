# Libraries and mcmc scripts #
library(tictoc) 
library(parallel)
numCores <- detectCores()

source("mcmc.recursive.R")
source("mcmc.hierarchical.R")

## Download and thin data from moveVis ##
library(moveVis)
data("whitestork_data")
rm(m) #detach large movestack data

df <- dplyr::select(df, time=timestamp, lon=`location-long`, lat=`location-lat`, name)
stork <- split(df, df$name)

# function to thin to ~12 hour fixes#
thin_fixes <- function(x){
        df <- x
        df.new <- data.frame()
        df.new <- rbind(df.new,df[1,])
        
        for(i in 1:length(df$time)){
                time.step=df$time[i]-df.new$time[nrow(df.new)]
                units(time.step) <- 'hours'
                
                if(time.step>10) df.new <- rbind(df.new, df[i,])
        }
        df.new
}
stork.list <- lapply(stork, thin_fixes)

## Alternatively, load thinned data as csv ##
     # csv file can be downloaded from  https://github.com/hannamcc/TARB
storks <- read.csv('stork_subset.csv') 
storks <- dplyr::select(storks, -X)
storks$time <- as.POSIXct(storks$time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
stork.list <- split(storks, storks$name)

###############################################3
## Setup ##
set.seed(420)
J=15
n.mcmc=100000 
n.burn=10000

##### Stage 1 #####
out.arr <- array(dim=c(J,n.mcmc-n.burn,3))

# Function to fit in parallel #
fit.s1 <- function(stork.list){
        fit.j <- s1.mcmc(times=stork.list[,1],S=stork.list[,2:3],n.mcmc=n.mcmc,q0=0.001,r0=1000,mu.0=0,s2.0=10, phi.tune=0.6)
        beta.vec <- fit.j$beta.save[-(1:n.burn)]
        s2.vec <- fit.j$s2.save[-(1:n.burn)]
        phi.vec <- fit.j$phi.save[-(1:n.burn)]
        list(pars=rbind(beta.vec, s2.vec, phi.vec), accept=fit.j$phi.accept)
}

# Fit model #
tic("recursive total")
tic("stage1")
stage1 <- mclapply(stork.list, fit.s1, mc.cores=numCores)
toc()

# organize output #
for(j in 1:15){
     out.arr[j,,1] <- stage1[[j]]$pars[1,]
     out.arr[j,,2] <- stage1[[j]]$pars[2,]
     out.arr[j,,3] <- stage1[[j]]$pars[3,]
}

##### Stage 2 #####
stage2 <- s2.mcmc(out.arr, q0=0.001,r0=1000, mu.0=0,s2.0=10, mu.1=0,s2.1=100, q1=0.001,r1=1000, q2=0.001,r2=1000, mu.2=0,s2.2=100, n.mcmc=n.mcmc)
toc() 

##### Hierarchical algortihm #####
tic("hierarch fit")
fit.hier <- hier.mcmc(y.list=stork.list, mu.1=0,s2.1=10, q1=0.001,r1=1000, q2=0.001,r2=1000, mu.2=0.2,s2.2=100, n.mcmc=n.mcmc,logsig.tune = rep(0.12,J), phi.tune = rep(0.7,J)) 
toc()


##### Create results figure #####
three.cols <- c("#000066", "#009966", "#99CC00")
par(mfrow=c(3,1), mar=c(3.5, 5, 1.5, 2.1),mgp=c(2.2, 1, 0),xpd=NA)

## beta ##
# hierarchical #
matplot(1:J,rowMeans(fit.hier$beta.save[,-(1:n.burn)]), type="n", xlim=c(0.75,18), ylim=c(0,0.03), xlab="Individual", ylab=expression(beta[j]),xaxt='n') 
axis(side=1, at=c(1:15,17.3), labels=c(1:15, expression(mu[beta])), tck=-0.03)
segments(1:J,apply(fit.hier$beta.save[,-(1:n.burn)],1,quantile,0.025),1:J,apply(fit.hier$beta.save[,-(1:n.burn)],1,quantile,0.975),lwd=2, col=three.cols[1])
points(1:J, rowMeans(fit.hier$beta.save[,-(1:n.burn)]), pch=1, col=three.cols[1])
# stage 1 #
segments(1:J+0.15,apply(out.arr[,,1],1,quantile,0.025),1:J+0.15,apply(out.arr[,,1],1,quantile,0.975),lwd=2, col=three.cols[2])
points(1:J+0.15, rowMeans(out.arr[,,1]), pch=1, col=three.cols[2])
# stage 2 #
segments(1:J+0.3,apply(stage2$beta.save[,-(1:n.burn)],1,quantile,0.025),1:J+0.3,apply(stage2$beta.save[,-(1:n.burn)],1,quantile,0.975),lwd=2, col=three.cols[3])
points(1:J+0.3, rowMeans(stage2$beta.save[,-(1:n.burn)]), pch=1, col=three.cols[3])
# population-level #
segments(16,0,16,0.03, lty="dashed",col="gray")
points(17, mean(fit.hier$mu.beta.save[-(1:n.burn)]), pch=1, col=three.cols[1])
segments(17, quantile(fit.hier$mu.beta.save[-(1:n.burn)], 0.025), 17, quantile(fit.hier$mu.beta.save[-(1:n.burn)], 0.975), lwd=2, col=three.cols[1])
points(17.3, mean(stage2$mu.beta.save[-(1:n.burn)]), pch=1, col=three.cols[3])
segments(17.3, quantile(stage2$mu.beta.save[-(1:n.burn)], 0.025), 17.3, quantile(stage2$mu.beta.save[-(1:n.burn)], 0.975), lwd=2, col=three.cols[3])
legend("topright", c("Full", "Stage 1", "Stage 2"),col=three.cols,lty=1,lwd=2, bty="n", seg.len=1.2)
text("A", x=-1.2,y=0.032,cex=1.5)

## log(sigma) ##
# hierarchical #
matplot(1:J,rowMeans(fit.hier$log.sig.save[,-(1:n.burn)]), type="n",xlim=c(0.75,18), ylim=c(-3.0,-2.0), xlab="Individual", ylab=expression(paste("log(",sigma[j],")")),xaxt='n')
axis(side=1, at=c(1:15,17.3), labels=c(1:15, expression(mu[sigma])), tck=-0.03)
segments(1:J,apply(fit.hier$log.sig.save[,-(1:n.burn)],1,quantile,0.025),1:J,apply(fit.hier$log.sig.save[,-(1:n.burn)],1,quantile,0.975),lwd=2, col=three.cols[1])
points(1:J, rowMeans(fit.hier$log.sig.save[,-(1:n.burn)]), pch=1,  col=three.cols[1])
# stage 1 #
segments(1:J+0.15,apply(log(sqrt(out.arr[,,2])),1,quantile,0.025),1:J+0.15,apply(log(sqrt(out.arr[,,2])),1,quantile,0.975),lwd=2,  col=three.cols[2])
points(1:J+0.15, rowMeans(log(sqrt(out.arr[,,2]))), pch=1, col=three.cols[2])
# stage 2 #
segments(1:J+0.3,apply(log(sqrt(stage2$s2.save))[,-(1:n.burn)],1,quantile,0.025),1:J+0.3,apply(log(sqrt(stage2$s2.save))[,-(1:n.burn)],1,quantile,0.975),lwd=2,  col=three.cols[3])
points(1:J+0.3, rowMeans(log(sqrt(stage2$s2.save))[,-(1:n.burn)]), pch=1,  col=three.cols[3])
# population-level #
segments(16,-3.0,16,-2.0, lty="dashed",col="gray")
points(17, mean(fit.hier$mu.sigma.save[-(1:n.burn)]), pch=1, col=three.cols[1])
segments(17, quantile(fit.hier$mu.sigma.save[-(1:n.burn)], 0.025), 17, quantile(fit.hier$mu.sigma.save[-(1:n.burn)], 0.975), lwd=2, col=three.cols[1])
points(17.3, mean(stage2$mu.sigma.save[-(1:n.burn)]), pch=1, col=three.cols[3])
segments(17.3, quantile(stage2$mu.sigma.save[-(1:n.burn)], 0.025), 17.3, quantile(stage2$mu.sigma.save[-(1:n.burn)], 0.975), lwd=2, col=three.cols[3])
legend("topright", c("Full", "Stage 1", "Stage 2"),col=three.cols,lty=1,lwd=2, bty="n", seg.len=1.2)
text("B", x=-1.2,y=-1.9,cex=1.5)

## phi ##
# hierarchical #
matplot(1:J,rowMeans(fit.hier$phi.save[,-(1:n.burn)]), type="n", xlim=c(0.75,18),ylim=c(0.25,1.75), xlab="Individual", ylab=expression(phi[j]),xaxt='n')
axis(side=1, at=c(1:15), labels=c(1:15), tck=-0.03)
segments(1:J,apply(fit.hier$phi.save[,-(1:n.burn)],1,quantile,0.025),1:J,apply(fit.hier$phi.save[,-(1:n.burn)],1,quantile,0.975),lwd=2, col=three.cols[1])
points(1:J, rowMeans(fit.hier$phi.save[,-(1:n.burn)]), pch=1, col=three.cols[1])
# stage 1 #
segments(1:J+0.15,apply(out.arr[,,3],1,quantile,0.025),1:J+0.15,apply(out.arr[,,3],1,quantile,0.975),lwd=2, col=three.cols[2])
points(1:J+0.15, rowMeans(out.arr[,,3]), pch=1, col=three.cols[2])
# stage 2 #
segments(1:J+0.3,apply(stage2$phi.save[,-(1:n.burn)],1,quantile,0.025),1:J+0.3,apply(stage2$phi.save[,-(1:n.burn)],1,quantile,0.975),lwd=2, col=three.cols[3])
points(1:J+0.3, rowMeans(stage2$phi.save[,-(1:n.burn)]), pch=1, col=three.cols[3])
legend("topright", c("Full", "Stage 1", "Stage 2"),col=three.cols,lty=1,lwd=2, bty="n", seg.len=1.2)
text("C", x=-1.2,y=1.9,cex=1.5)
