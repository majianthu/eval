library(copula)
library(mnormt)
library(corrplot)

library(copent) # Copula Entropy
library(energy) # Distance Correlation
library(dHSIC) # Hilbert-Schmidt Independence Criterion
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


n = 9
ce = rep(0,n)
ktau = rep(0,n) # Kendall's tau
dcor = rep(0,n) # Distance Correlation
dhsic = rep(0,n)  # Hilbert-Schmidt Independence Criterion
hhg1 = rep(0,n)  # Heller-Heller-Gorfine Tests, Pearson chi-squared statistic
hhg2 = rep(0,n)  # Heller-Heller-Gorfine Tests, likelihood ratio statistic
ball = rep(0,n) # Ball correlation
qad = rep(0,n) # Quantification of Asymmetric Dependence
bet = rep(0,n) # Binary Expansion Testing
nns = rep(0,n) # Nonlinear Nonparametric Statistics
subcop = rep(0,n) #supremum dependence
mdm = rep(0,n) # mutual independence measure

for (i in 1:n){
  # normal distribution
  rho <- (i-1)/10
  r12 = 0.8; r34 = 0.75
  sigma <- matrix(c(1,r12,rho,rho,
                    r12,1,rho,rho,
                    rho,rho,1,r34,
                    rho,rho,r34,1),4,4)
  x <- rmnorm(800,c(0,0,0,0),sigma)
  x1 = x[,c(1,2)]; x2 = x[,c(3,4)]

  ce[i] = copent(x) - copent(x1) - copent(x2)
  dcor[i] = dcor(x1,x2)
  dhsic[i] = dhsic(x1,x2)$dHSIC
  Dx = as.matrix(dist((x1), diag = TRUE, upper = TRUE))
  Dy = as.matrix(dist((x2), diag = TRUE, upper = TRUE))
  hhg = hhg.test(Dx,Dy, nr.perm = 500)
  hhg1[i] = hhg$sum.chisq
  hhg2[i] = hhg$sum.lr
  ball[i] = bcor(x1,x2)
  qad[i] = qad(x1,x2)$`q(X,Y)`
  bet[i] = MaxBET(x, 3, index = list(c(1,2),c(3,4)))$z.statistic
  nns[i] = NNS.dep(x1,x2)$Dependence
  mdm[i] = mdc(x1,x2)
}

#joint
joint1 = cbind(ce,dcor,dhsic,hhg1,hhg2,ball,qad,bet,mdm,nns)
x11()
corrplot(cor(joint1), method = "shade", order = "hclust", col = COL2(n=200))

#### plot
labels1 = seq(0,0.8, by = 0.1); xlab1 = "rho" # normal and normal copula
w1 = 5; h1 = 5
# ce
x11(width = w1, height = h1)
plot(ce, xlab = xlab1, ylab = "Copula Entropy", xaxt = 'n')
axis(side = 1, at = c(1:n), labels = labels1)
lines(ce)

# dcor
x11(width = w1, height = h1)
plot(dcor, xlab = xlab1, ylab = "dCor", xaxt = 'n')
axis(side = 1, at = c(1:n), labels = labels1)
lines(dcor)

# dhsic
x11(width = w1, height = h1)
plot(dhsic, xlab = xlab1, ylab = "dHSIC", xaxt = 'n')
axis(side = 1, at = c(1:n), labels = labels1)
lines(dhsic)

# hhg
x11(width = w1, height = h1)
plot(hhg1, xlab = xlab1, ylab = "HHG.chisquared", xaxt = 'n')
axis(side = 1, at = c(1:n), labels = labels1)
lines(hhg1)

x11(width = w1, height = h1)
plot(hhg2, xlab = xlab1, ylab = "HHG.lr", xaxt = 'n')
axis(side = 1, at = c(1:n), labels = labels1)
lines(hhg2)

# ball
x11(width = w1, height = h1)
plot(ball, xlab = xlab1, ylab = "Ball", xaxt = 'n')
axis(side = 1, at = c(1:n), labels = labels1)
lines(ball)

# qad
x11(width = w1, height = h1)
plot(qad, xlab = xlab1, ylab = "qad", xaxt = 'n')
axis(side = 1, at = c(1:n), labels = labels1)
lines(qad)

# BET
x11(width = w1, height = h1)
plot(bet, xlab = xlab1, ylab = "BET", xaxt = 'n')
axis(side = 1, at = c(1:n), labels = labels1)
lines(bet)

# NNS
x11(width = w1, height = h1)
plot(nns, xlab = xlab1, ylab = "NNS", xaxt = 'n')
axis(side = 1, at = c(1:n), labels = labels1)
lines(nns)

# mdm
x11(width = w1, height = h1)
plot(mdm, xlab = xlab1, ylab = "MDC", xaxt = 'n')
axis(side = 1, at = c(1:n), labels = labels1)
lines(mdm)

