######################################################
### Example code: fit bernoulli glmm to simulated  ###
###    occurrence data with two-stage procedure    ###
######################################################
# Supporting information: Hierarchical models for 
#    hierarchical computing in ecology (2020)
# McCaslin, Feuka, and Hooten

source("bernoulli_glmm/mcmc.bern.R") # R script containing MCMC algorithm functions

### Simulate data ####################################
set.seed(101)

logit.p <- rnorm(10,1.5,1) # simulate logit(p_j) for J=10
p <- logit.inv(logit.p)

pres<- list()
for(j in 1:10){
     pres[[j]] <- rbinom(20,1,p[j]) # simulate N=20 occurrences for each site j
}

sim.dat <- cbind.data.frame(site=rep(1:10,each=20),plot=rep(1:20,times=10), pres=unlist(pres))
sim.list <- split(sim.dat, sim.dat$site)

### Fit model using TARB #############################
# Stage 1 #
alpha <- 5
beta <- 2.5
sim.s1 <- lapply(sim.list, function(x) rbeta(10000,sum(x$pres)+alpha, sum(1-x$pres)+beta))

sim.matrix <- matrix(nrow=10,ncol=10000) 
for(j in 1:10) sim.matrix[j,] <- sim.s1[[j]]

# Stage 2 #
sim.s2 <- bern.s2(mat=sim.matrix, alpha=5, beta=2.5, n.mcmc=10000)

### For comparison, fit model in single algorithm ####################
sim.full <- bern.full(y.list=sim.list, tune=rep(.8,J), n.mcmc=10000)
sim.full$accept

### Create plot to compare results of 2 approaches ###################
three.cols <- c("#000066", "#009966", "#99CC00")
J=10
matplot(1:J,rowMeans(logit(sim.matrix)), type="n", xlim=c(0.75,13), ylim=c(-1,6), xlab="Individual", ylab="logit(p)",xaxt='n') 
axis(side=1, at=c(1:10,12.3), labels=c(1:10, expression(mu)), tck=-0.03)
# single algorithm #
segments(1:J,apply(sim.full$logit.p.save[,-(1:1000)],1,quantile,0.025),1:J,apply(sim.full$logit.p.save[,-(1:1000)],1,quantile,0.975),lwd=2, col=three.cols[1])
points(1:J, rowMeans(sim.full$logit.p.save[,-(1:1000)]), pch=1, col=three.cols[1])
# TARB stage 1 #
segments(1:J+0.15,apply(logit(sim.matrix),1,quantile,0.025),1:J+0.15,apply(logit(sim.matrix),1,quantile,0.975),lwd=2, col=three.cols[2])
points(1:J+0.15, rowMeans(logit(sim.matrix)), pch=1, col=three.cols[2])
# TARB stage 2 #
segments(1:J+0.3,apply(sim.s2$logit.p.save[,-(1:1000)],1,quantile,0.025),1:J+0.3,apply(sim.s2$logit.p.save[,-(1:1000)],1,quantile,0.975),lwd=2, col=three.cols[3])
points(1:J+0.3, rowMeans(sim.s2$logit.p.save[,-(1:1000)]), pch=1, col=three.cols[3])
# group-level parameters #
segments(11,-3,11,5, lty="dashed",col="gray")
points(12.3, mean(sim.s2$mu.save[-(1:1000)]), pch=1, col=three.cols[3])
segments(12.3, quantile(sim.s2$mu.save[-(1:1000)], 0.025), 12.3, quantile(sim.s2$mu.save[-(1:1000)], 0.975), lwd=2, col=three.cols[3])
points(12, mean(sim.full$mu.save[-(1:1000)]), pch=1, col=three.cols[1])
segments(12, quantile(sim.full$mu.save[-(1:1000)], 0.025), 12, quantile(sim.full$mu.save[-(1:1000)], 0.975), lwd=2, col=three.cols[1])
# "truth" (simulated values) #
points(1:J, logit.p, col="red", cex=0.5, pch=2)

legend("topright", c("Full", "Stage 1", "Stage 2"),col=three.cols,lty=1,lwd=2, bty="n", seg.len=1.2,y.intersp=.75)
