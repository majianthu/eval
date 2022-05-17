library(copula)
library(mnormt)
library(corrplot)

library(copent) # Copula Entropy
library(energy) # Distance Correlation
library(dHSIC) # Hilbert-Schmidt Independence Criterion
## for additional tests
library(HHG) # Heller-Heller-Gorfine Tests of Independence
library(independence) # Hoeffding's D test or Bergsma-Dassios T* sign covariance
library(Ball) # Ball correlation
library(qad) # Quantification of Asymmetric Dependence
library(BET) # Binary Expansion Testing
library(MixedIndTests) # Cramer-von Mises statistics
library(NNS) # Nonlinear Nonparametric Statistics
library(subcopem2D) # supremum dependence
library(EDMeasure) # Mutual Independence Measure
library(FOCI) # Dependence Coefficient 

n = 10
ce = rep(0,n)
ktau = rep(0,n) # Kendall's tau
dcor = rep(0,n) # Distance Correlation
dhsic = rep(0,n)  # Hilbert-Schmidt Independence Criterion
hhg1 = rep(0,n)  # Heller-Heller-Gorfine Tests, Pearson chi-squared statistic
hhg2 = rep(0,n)  # Heller-Heller-Gorfine Tests, likelihood ratio statistic
hoeff = rep(0,n)  # Hoeffding's D test 
bdtau = rep(0,n) # Bergsma-Dassios T* sign covariance
ball = rep(0,n) # Ball correlation
qad = rep(0,n) # Quantification of Asymmetric Dependence
bet = rep(0,n) # Binary Expansion Testing
mixed = rep(0,n) # Cramer-von Mises statistics
nns = rep(0,n) # Nonlinear Nonparametric Statistics
subcop = rep(0,n) #supremum dependence
mdm = rep(0,n) # mutual independence measure
codec = rep(0,n) # dependence coefficient

for (i in 1:n){
  # normal distribution
  # rho <- (i-1)/n
  # sigma <- matrix(c(1,rho,rho,1),2,2)
  # x <- rmnorm(800,c(0,0),sigma)
  
  # copula 
  # mv.NE <- mvdc(normalCopula( (i-1)/n ), c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = 2)))
  # mv.NE <- mvdc(claytonCopula(i), c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = 2)))
  # mv.NE <- mvdc(frankCopula(i), c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = 2)))
  mv.NE <- mvdc(gumbelCopula(i), c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = 2)))
  x <- rMvdc(800, mv.NE)

  ce[i] = copent(x)
  ktau[i] = cor(x, method = "kendall")[1,2]
  dcor[i] = dcor(x[,1],x[,2])
  dhsic[i] = dhsic(x[,1],x[,2])$dHSIC
  Dx = as.matrix(dist((x[,1]), diag = TRUE, upper = TRUE))
  Dy = as.matrix(dist((x[,2]), diag = TRUE, upper = TRUE))
  hhg1[i] = hhg.test(Dx,Dy, nr.perm = 500)$sum.chisq
  hhg2[i] = hhg.test(Dx,Dy, nr.perm = 500)$sum.lr
  hoeff[i] = hoeffding.D.test(x[,1],x[,2])$Dn
  bdtau[i] = tau.star.test(x[,1],x[,2])$Tn
  ball[i] = bcor(x[,1],x[,2])
  qad[i] = qad(x[,1],x[,2])$`q(X,Y)`
  bet[i] = MaxBET(x, 3, index = list(1,2))$z.statistic
  mixed[i] = TestIndCopula(x)$stat$cvm
  nns[i] = NNS.dep(x[,1],x[,2])$Dependence
  subcop[i] = dependence(x)[1,2,2]
  mdm[i] = mdc(x[,1],x[,2])
  codec[i] = codec(x[,1],x[,2])
}

#joint
joint1 = cbind(ce,ktau,dcor,dhsic,hhg1,hhg2,hoeff,bdtau,ball,qad,bet,mixed,subcop,mdm,codec,nns)
# x11()
corrplot(cor(joint1), method = "shade", order = "hclust", col = COL2(n=200))

#### plot
# labels1 = seq(0,0.9, by = 0.1); xlab1 = "rho" # normal and normal copula
labels1 = 1:10; xlab1 = "alpha" # Archimedean copula
w1 = 5; h1 = 5
# ce
x11(width = w1, height = h1)
plot(ce, xlab = xlab1, ylab = "Copula Entropy", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(ce)

# kendall's tau
x11(width = w1, height = h1)
plot(ktau, xlab = xlab1, ylab = "Kendall's tau", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(ktau)

# dcor
x11(width = w1, height = h1)
plot(dcor, xlab = xlab1, ylab = "dCor", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(dcor)

# dhsic
x11(width = w1, height = h1)
plot(dhsic, xlab = xlab1, ylab = "dHSIC", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(dhsic)

# hhg
x11(width = w1, height = h1)
plot(hhg1, xlab = xlab1, ylab = "HHG.chisquared", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(hhg1)

x11(width = w1, height = h1)
plot(hhg2, xlab = xlab1, ylab = "HHG.lr", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(hhg2)

# independence
x11(width = w1, height = h1)
plot(hoeff, xlab = xlab1, ylab = "Hoeffding", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(hoeff)

# BD tau
x11(width = w1, height = h1)
plot(bdtau, xlab = xlab1, ylab = "BD tau", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(bdtau)

# ball
x11(width = w1, height = h1)
plot(ball, xlab = xlab1, ylab = "Ball", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(ball)

# qad
x11(width = w1, height = h1)
plot(qad, xlab = xlab1, ylab = "qad", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(qad)

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

# mdm
x11(width = w1, height = h1)
plot(mdm, xlab = xlab1, ylab = "MDC", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(mdm)

# CODEC
x11(width = w1, height = h1)
plot(codec, xlab = xlab1, ylab = "CODEC", xaxt = 'n')
axis(side = 1, at = c(1:10), labels = labels1)
lines(codec)

