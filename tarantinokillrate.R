## BAYESIAN PROJECT about TARANTINO'S Movie KillRate

install.packages("PearsonDS")
install.packages("geoR")
install.pacages("nleqslv")
library(geoR)
library(PearsonDS)
library(HDInterval)
library(nleqslv)


movie <- read.csv('~/SJSU_ClASSES/Math264_Bayesian/project/tarantinodata.csv')
movie_history <- read.csv('~/SJSU_ClASSES/Math264_Bayesian/project/tarantino_historydata.csv')
length(movie_history $KillRate)  ## 9
max(movie_history $KillRate)
min(movie_history $KillRate)
 m <- mean(movie_history $KillRate)
ssq <- var(movie_history $KillRate)

length(movie$KillRate)  ## 
max(movie$KillRate)
min(movie$KillRate)
m_2 <- mean(movie$KillRate)
ssq_2 <- var(movie$KillRate)

##
#pseudo- likelihood of dispersion
(ssq -m)/m^2
(ssq_2 -(m_2))/(m_2)^2


## mle
temp <- 0
v <- movie_history$KillRate[-8]

phi <- function(x){
	for (i in 1:length(v)){
		for(j in 0:(v[i]-1)){
			temp <- temp + j/(1+x*j)
		}
	}
	1/length(v)*temp + x^(-2)*log(1+x*mean(v)) - mean(v)*(mean(v)+x^(-1))/(1+x*mean(v))
}

xstart <- c(0.1,0.99)
nleqslv(xstart, phi, control=list(trace=1,btol=.01,delta="cauchy"))


## Plot of Prior vs Posterior

curve(dpearsonVI(x,247, 14,1, 1/0.67),xlim = c(0,100), add=TRUE,col="violetred",ann=FALSE)
curve(dpearsonVI(x,16, 2.05,1, 1/0.67), xlim = c(0,100),add=TRUE, col="royalblue", ann=FALSE)
title(xlab="lambda", ylab="Density")
legend(65, 0.05, legend = c("Prior of lambda", "Posterior of lambda"), lty = 1:1, lwd = c(1.5, 1.2), col = c("royalblue","violetred"), y.intersp = 1.2, bty = "n", cex = 1)

curve(dpearsonVI(x,247, 14,1, 1/0.67),xlim = c(0,100), add=TRUE,col="violetred",ann=FALSE)       
lines(c(15.6,15.6),c(0,1), lty=2,col="gray")
lines(c(46.07,46.07),c(0,1), lty=2,col="gray")
lines(c(25.5,25.5),c(0,1), lty=2,col="gray")
axis(1, at=c(0,15.6,20,25.5,40,46.1,60,80,100),labels=c(0,15.6,20,25.5,40,46.1,60,80,100))

curve(dpearsonVI(x,253, 15.5,1, 1/0.67),xlim = c(0,100), add=TRUE,col="violetred",ann=FALSE)
curve(dpearsonVI(x,16, 2.05,1, 1/0.67), xlim = c(0,100),add=TRUE, col="royalblue", ann=FALSE)
title(xlab="lambda", ylab="Density")
legend(80, 0.05, legend = c("Prior", "Posterior"), lty = 1:1, lwd = c(1.5, 1.2), col = c("royalblue","violetred"), y.intersp = 1.2, bty = "n", cex = 1)
       
       
       
## HDI interval

hdi(qpearsonVI, a= 247, b=14, location= 1, scale = 1/0.67, credMass= 0.95)

#   lower    upper 
# 15.58805 46.07395 


## posterior predictive

set.seed(210)
rate <- rpearsonVI(1000, 247, 14,1, 1/0.67)
y<- rnbinom(1000, size=1/0.67, mu =rate)

hist(y, breaks=seq(0,170,by=1), right=FALSE, xlim=c(0,170),freq=FALSE,col = "cornflowerblue", border = "wheat",  ann=FALSE)
axis(1, at=seq(0,130,by=5), labels=seq(0,130,by=5), cex=0.5)
title(xlab="Body Count", ylab="Densiy", main="Predictive Body Count In An Hour")




set.seed(210)
y <- rep(0,1000)
rate <- rpearsonVI(1000, 247, 14,1, 1/0.67)
for (i in 1:1000) {
y[i] <- rnbinom(1, size=1/0.67 mu =rate[i])
}
hist(y, breaks=100, xlim=c(0,100),freq=FALSE,col = "cornflowerblue", border = "wheat",,  ann=FALSE)
title(xlab="Body Count", ylab="Densiy", main="Next Movie Body Count Prediction")



plot(factor(y), space=1, col="cornflowerblue")

#lambdas = rgamma(100, shape=2, scale=3)
#samples = rep(0, 100)
#for (i in 1:100)
#samples[i] = rpois(1, lambdas[i])
  
  

##

## posterior predictive checking

n.obs <- length(movie$KillRate[-9])         # number of observations
y.bar <- mean(movie$KillRate[-9])          # sample mean
s.sq <- var(movie$KillRate[-9])            # sample variance

set.seed(2109) 
m <- 10000
rate <- rpearsonVI(m, 247, 14,1, 1/0.67)
yrep <- mapply(rnbinom, n= 8, size=1/0.67, mu=rate)


# minimum
obs.min <- min(movie$KillRate[-9])      # observed minimum
sim.min <- apply(yrep, 2, min)   # simulated minimum
pval_min <- length(sim.min[sim.min >= obs.min]) / m
hist(sim.min, freq=FALSE,right=FALSE, col = "cornflowerblue", border = "wheat", xlim=c(0,50), ann=FALSE)
lines(rep(obs.min, 2), c(0, 1), col = "red", lwd = 1.5)

# maximum
obs.max <- max(movie$KillRate[-9])      # observed 
sim.max <- apply(yrep, 2, max)   # simulated
pval_max <- length(sim.max[sim.max <= obs.max]) / m
hist(sim.max, freq=FALSE,right=FALSE,col = "cornflowerblue", border = "wheat",, xlim=c(0,300), ann=FALSE)
lines(rep(obs.max, 2), c(0, 1), col = "red", lwd = 1.5)

# mean
obs.mean <- mean(movie$KillRate[-9])      # observed 
sim.mean <- apply(yrep, 2, mean)   # simulated 
pval_mean <- length(sim.max[sim.mean >= obs.mean]) / m
hist(sim.mean, freq=FALSE,right=FALSE,col = "cornflowerblue", border = "wheat",, xlim=c(0,150), ann=FALSE)
lines(rep(obs.mean, 2), c(0, 1), col = "red", lwd = 1.5)

# sd
obs.sd <- sqrt(var(movie$KillRate[-9]))      # observed 
sim.sd <- sqrt(apply(yrep, 2, var))   # simulated
pval_sd <- length(sim.sd[sim.sd >= obs.sd]) / m
hist(sim.sd, freq=FALSE,right=FALSE, col = "cornflowerblue", border = "wheat",, xlim=c(0,100), ann=FALSE)
lines(rep(obs.sd, 2), c(0, 1), col = "red", lwd = 1.5)

# var - mean
obs.vm <- var(movie$KillRate[-9])- mean(movie$KillRate[-9])      # observed 
sim.vm <- apply(yrep, 2, var)- apply(yrep, 2, mean)   # simulated
pval_vm <- length(sim.vm[sim.vm >= obs.vm]) / m
hist(sim.vm, freq=FALSE,right=FALSE, col = "cornflowerblue", border = "wheat",, xlim=c(0,5000), ann=FALSE)
lines(rep(obs.vm, 2), c(0, 1), col = "red", lwd = 1.5)


mu <- apply(yrep,2,mean)
obs.stat <- obs.sd^2 - sim.mean
sim.stat <- (sim.sd)^2 - sim.mean
pval <- length(sim.stat[sim.stat >= obs.stat]) / m
plot(obs.stat, sim.stat, cex=0.2,col='royalblue', ann=FALSE)
abline(a = 0, b = 1, lty = 2, col = "gray30")


# Longtail
indx <- c(1,8)
mu <- apply(yrep,2,mean)
obs.ord <- sort(movie$KillRate[-9])[indx]
obs.stat <- abs(obs.ord[2] - mu) - abs(obs.ord[1] - mu)
sim.ord <- apply(yrep, 2, function(y) sort(y)[indx])
sim.stat <- abs(sim.ord[2, ] - mu) - abs(sim.ord[1, ] - mu)
pval_lt <- length(sim.stat[sim.stat >= obs.stat]) / m
hist(sim.stat, freq=FALSE,col = "cornflowerblue", border = "wheat",, xlim=c(-100,100), ann=FALSE)
lines(rep(obs.stat, 2), c(0, 1), col = "red", lwd = 1.5)

plot(obs.stat, sim.stat, cex=0.2,col='royalblue', ann=FALSE)
abline(a = 0, b = 1, lty = 2, col = "gray30")

## Movie length

length <- rep(0,1000)
sgm2 <- rinvchisq(1000,17,0.1391)
for (i in 1:1000) {
length[i] <- rnorm(1, 2, sqrt(sgm2[i]))
}
hist(length, breaks=30, xlim=c(0,4),freq=FALSE,col = "cornflowerblue", border = "wheat",,  ann=FALSE)
title(xlab="Length(Hour)", ylab="Densiy", main="Next Movie Length Prediction")
f <- density(length, from = 0, to = 4, width=0.8)
lines(f$x, f$y, col = "violetred", lty = 2, lwd = 1.2)

legend(3, 0.8, legend = c("Kernel"), lty = 1, lwd = c( 1.2),
       col = c("violetred"), bty = "n", cex = 1)
       




## Next Movie prediction
set.seed(210)
sample <- rep(0,1000)
y <- rep(0,1000)
length <- rep(0,1000)
rate <- rpearsonVI(1000, 261.6, 15,1, 1/0.67)
sgm2 <- rinvchisq(1000,17,0.1391)
for (i in 1:1000) {
y[i] <- rnbinom(1, size=1/0.67, mu =rate[i])
length[i] <- rnorm(1, 2, sqrt(sgm2[i]))
sample[i] <- y[i]*length[i]
}
hist(sample, breaks=300, xlim=c(0,200),freq=FALSE,col = "cornflowerblue", border = "wheat",,  ann=FALSE)
title(xlab="Body Count", ylab="Densiy", main="Next Movie Body Count Prediction")


set.seed(210)
sample <- rep(0,10000)
y <- rep(0,10000)
length <- rep(0,10000)
rate <- rpearsonVI(10000, 261.6, 15,1, 1/0.67)
sgm2 <- rinvchisq(10000,17,0.1391)
for (i in 1:10000) {
y[i] <- rnbinom(1, size=1/0.67, mu =rate[i])
length[i] <- rnorm(1, 2, sqrt(sgm2[i]))
sample[i] <- y[i]*length[i]
}
hist(round(sample), breaks=1000, xlim=c(0,200),freq=FALSE,col = "cornflowerblue", border = "wheat",right=FALSE,  ann=FALSE)
title(xlab="Body Count", ylab="Densiy", main="Next Movie Body Count Prediction")

## 95% HPD
hdi(round(sample),credMass= 0.95)
#lower upper 
#    0   169 
hdi(round(sample), 0.8)
#lower upper 
#    0    93 



marginal <- function(x) {
	nc <- gamma(n.obs / 2) / (gamma((n.obs - 1) / 2) * sqrt((n.obs - 1) * pi * s.sq / n.obs))
	nc * (1 + n.obs * (x - y.bar)^2 / ((n.obs - 1) * s.sq))^(-n.obs / 2)
}


