######################################################
### Example code: fit bernoulli glmm to cheatgrass ###
###    occurrence data with two-stage procedure    ###
######################################################
# Supporting information: Hierarchical models for 
#    hierarchical computing in ecology (2020)
# McCaslin, Feuka, and Hooten

### Import data ######################################################
     # This code chunk creates the data used in this example, collected by:
     # Pearson, Dean E. et al. (2020), Data from: Are exotic plants more abundant in the introduced           versus native range?, v3, Dryad, Dataset, https://doi.org/10.5061/dryad.r7v91

# This vector allows user to replicate the response data 
y.vec <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,0,0,1,1,1,0,1,1,1,0,1,0,0,1,0,1,1,1,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,0,0,1,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,0,0,0,1,1,1,1,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1)

brotec <- cbind.data.frame(SITE_ID=rep(1:20,each=20),PLOT=rep(1:20,times=20), y=y.vec)

### Fit model using TARB ################################################

# Split dataframe into list of dataframes #
brotec.list <- split(brotec, brotec$SITE_ID)

# Stage 1 #
     # sample directly from posterior distribution for each site #
     # Note stage 1 posteriors are only to use as proposals in stage 2,
     #    not for inference
set.seed(130)
J=20
alpha=.5
beta=.5 

brotec.s1 <- list()
for(j in 1:J){
     brotec.s1[[j]] <- rbeta(100000, sum(brotec.list[[j]]$y) + alpha, sum(1-brotec.list[[j]]$y) + beta)
}

# Unlist samples into matrix #
s1.matrix <- matrix(nrow=J,ncol=100000)
for(j in 1:J) s1.matrix[j,] <- brotec.s1[[j]]

# Stage 2 #
source("bernoulli_glmm/mcmc.bern.R")
brotec.s2 <- bern.s2(mat=s1.matrix, alpha=.5, beta=.5, n.mcmc=100000)

### Model results ####################################################
library(scales)
nburn=20000 #remove first 20% of samples as burn in

plot(density(brotec.s2$mu.save[-(1:nburn)]), xlim=c(-4,10), ylim=c(0,0.9), 
     lwd=2,main="", xlab=expression(logit(p[j])))

for(i in 20:100){
curve(dnorm(x, brotec.s2$mu.save[1000*i], sqrt(brotec.s2$s2.save[1000*i])), col=alpha("gray", 0.2), add=T)
}

for(j in 1:J) {
        lines(density(brotec.s2$logit.p.save[j,-(1:nburn)]),
              col=alpha("#99CC00", 0.35), lwd=1.5)
}
legend("topright", c(expression(logit(p[j])), expression(mu), expression(N(mu, sigma^2))),col=c(alpha("#99CC00", 0.35), "black", alpha("gray", 0.5)),lty=1,lwd=c(2,1.5), bty="n", seg.len=1.2,y.intersp=.8)

######################################################################
### For comparison, fit model in single algorithm ####################
set.seed(130)
tune.vec <- c(rep(3,4), 2,3,2,3,3,2,3,2,2,3,3,3,1,1,1,2) #adjust tuning parameters for each site
brotec.full <- bern.full(y.list=brotec.list, tune=tune.vec, n.mcmc=100000)

brotec.full$accept #this algorithm must be tuned (optimal acceptance rate ~0.4-0.5)


### Create plot to compare results of 2 approaches ###################
three.cols <- c("#000066", "#009966", "#99CC00")
matplot(1:J,rowMeans(logit(s1.matrix)), type="n", xlim=c(0.75,23), ylim=c(-4,12), xlab="Individual", ylab="logit(p)",xaxt='n') 
axis(side=1, at=c(1:20,22.3), labels=c(1:20, expression(mu)), tck=-0.03)
# single algorithm #
segments(1:J,apply(brotec.full$logit.p.save[,-(1:nburn)],1,quantile,0.025),1:J,apply(brotec.full$logit.p.save[,-(1:nburn)],1,quantile,0.975),lwd=2, col=three.cols[1])
points(1:J, rowMeans(brotec.full$logit.p.save[,-(1:nburn)]), pch=1, col=three.cols[1])
# TARB stage 1 #
segments(1:J+0.15,apply(logit(s1.matrix),1,quantile,0.025),1:J+0.15,apply(logit(s1.matrix),1,quantile,0.975),lwd=2, col=three.cols[2])
points(1:J+0.15, rowMeans(logit(s1.matrix)), pch=1, col=three.cols[2])
# TARB stage 2 #
segments(1:J+0.3,apply(brotec.s2$logit.p.save[,-(1:nburn)],1,quantile,0.025),1:J+0.3,apply(brotec.s2$logit.p.save[,-(1:nburn)],1,quantile,0.975),lwd=2, col=three.cols[3])
points(1:J+0.3, rowMeans(brotec.s2$logit.p.save[,-(1:nburn)]), pch=1, col=three.cols[3])
# group-level parameters #
segments(21,-4,21,12, lty="dashed",col="gray")
points(22, mean(brotec.full$mu.save[-(1:nburn)]), pch=1, col=three.cols[1])
segments(22, quantile(brotec.full$mu.save[-(1:nburn)], 0.025), 22, quantile(brotec.full$mu.save[-(1:nburn)], 0.975), lwd=2, col=three.cols[1])
points(22.3, mean(brotec.s2$mu.save[-(1:nburn)]), pch=1, col=three.cols[3])
segments(22.3, quantile(brotec.s2$mu.save[-(1:nburn)], 0.025), 22.3, quantile(brotec.s2$mu.save[-(1:nburn)], 0.975), lwd=2, col=three.cols[3])

legend("topright", c("Full", "Stage 1", "Stage 2"),col=three.cols,lty=1,lwd=2, bty="n", seg.len=1.2,y.intersp=.5)

