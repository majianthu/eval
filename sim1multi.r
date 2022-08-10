library(copula)
library(mnormt)
library(corrplot)

library(copent) # Copula Entropy
library(dHSIC) # Hilbert-Schmidt Independence Criterion
library(BET) # Binary Expansion Testing
library(MixedIndTests) # Cramer-von Mises statistics
library(NNS) # Nonlinear Nonparametric Statistics
library(subcopem2D) # supremum dependence

n = 10
ce = rep(0,n)
dhsic = rep(0,n)  # Hilbert-Schmidt Independence Criterion
bet = rep(0,n) # Binary Expansion Testing
mixed = rep(0,n) # Cramer-von Mises statistics
nns = rep(0,n) # Nonlinear Nonparametric Statistics
subcop = rep(0,n) #supremum dependence

for (i in 1:n){
  # normal distribution
  # rho <- (i-1)/n
  # sigma <- matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3)
  # x <- rmnorm(800,c(0,0,0),sigma)
  # copula 
  # mv.NE <- mvdc(claytonCopula(i,dim = 3), c("norm", "exp","exp"), list(list(mean = 0, sd =2), list(rate = 2), list(rate = .5)))
  # mv.NE <- mvdc(frankCopula(i,dim = 3), c("norm", "exp","exp"), list(list(mean = 0, sd =2), list(rate = 2), list(rate = .5)))
  mv.NE <- mvdc(gumbelCopula(i,dim = 3), c("norm", "exp","exp"), list(list(mean = 0, sd =2), list(rate = 2), list(rate = .5)))
  x <- rMvdc(800, mv.NE)
  
  ce[i] = copent(x)
  dhsic[i] = dhsic(list(x1 = x[,1], x2 = x[,2], x3 = x[,3]))$dHSIC
  bet[i] = MaxBET(x, 3, index = list(1,2,3))$z.statistic
  mixed[i] = TestIndCopula(x)$stat$cvm
  nns[i] = NNS.dep(x)$Dependence[1,2]
  subcop[i] = dependence(x)[1,2,2]
}

#joint
joint1 = cbind(ce,dhsic,bet,mixed,subcop)
x11()
corrplot(cor(joint1), method = "shade", order = "hclust", col = COL2(n=200))

#### plot
labels1 = seq(0,0.9, by = 0.1); xlab1 = "rho" # normal and normal copula
w1 = 5; h1 = 5
# ce
x11(width = w1, height = h1)
plot(ce, xlab = xlab1, ylab = "Copula Entropy", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(ce)

# dhsic
x11(width = w1, height = h1)
plot(dhsic, xlab = xlab1, ylab = "dHSIC", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(dhsic)


# BET
x11(width = w1, height = h1)
plot(bet, xlab = xlab1, ylab = "BET", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(bet)

# MixedIndCopula
x11(width = w1, height = h1)
plot(mixed, xlab = xlab1, ylab = "Mixed", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(mixed)

# NNS
x11(width = w1, height = h1)
plot(nns, xlab = xlab1, ylab = "NNS", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(nns)

# subcop
x11(width = w1, height = h1)
plot(subcop, xlab = xlab1, ylab = "subcop", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(subcop)

