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
library(corrplot)
library(mnormt)

# red1 = read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv",sep = ";")
wine1 = as.matrix(read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv",sep = ";"))
n = dim(wine1)[1]
d = dim(wine1)[2] -1

ce = rep(0,d)
ktau = rep(0,d) # Kendall's tau
dcor = rep(0,d) # Distance Correlation
dhsic = rep(0,d)  # Hilbert-Schmidt Independence Criterion
hhg1 = rep(0,d)  # Heller-Heller-Gorfine Tests, Pearson chi-squared statistic
hhg2 = rep(0,d)  # Heller-Heller-Gorfine Tests, likelihood ratio statistic
ind = rep(0,d)  # Hoeffding's D test 
bdtau = rep(0,d) # Bergsma-Dassios T* sign covariance
ball = rep(0,d) # Ball correlation
qad = rep(0,d) # Quantification of Asymmetric Dependence
bet = rep(0,d) # Binary Expansion Testing
mixed = rep(0,d) # Cramer-von Mises statistics
nns = rep(0,d) # Nonlinear Nonparametric Statistics
subcop = rep(0,d) #supremum dependence
mdm = rep(0,d) # mutual independence measure
codec = rep(0,d) # dependence coefficient

l = 500
for (i in 1:d){
  x = wine1[1:l,c(i,12)]

  for(k in 1:20){
    x1 = x
    x1[,1] = x1[,1] + runif(l) * 0.0000001 * max(abs(x[,1]))
    x1[,2] = x1[,2] + runif(l) * 0.0000001 * max(abs(x[,2]))
    ce[i] = ce[i] + copent(x1)/20
  }
  
  ktau[i] = cor(x, method = "kendall")[1,2]
  dcor[i] = dcor(x[,1],x[,2])
  dhsic[i] = dhsic(x[,1],x[,2])$dHSIC
  Dx = as.matrix(dist((x[,1]), diag = TRUE, upper = TRUE))
  Dy = as.matrix(dist((x[,2]), diag = TRUE, upper = TRUE))
  hhg = hhg.test(Dx,Dy, nr.perm = 500)
  hhg1[i] = hhg$sum.chisq
  hhg2[i] = hhg$sum.lr
  ind[i] = hoeffding.D.test(x[,1],x[,2])$Dn
  bdtau[i] = tau.star.test(x[,1],x[,2])$Tn
  ball[i] = bcor(x[,1],x[,2])
  qad[i] = qad(x[,1],x[,2])$`q(X,Y)`
  bet[i] = MaxBET(as.matrix(x), 3, index = list(1,2))$z.statistic
  mixed[i] = TestIndCopula(x)$stat$cvm
  nns[i] = NNS.dep(x[,1],x[,2])$Dependence
  subcop[i] = dependence(x)[1,2,2]
  mdm[i] = mdc(x[,1],x[,2])
  codec[i] = codec(x[,1],x[,2])
}

#joint
joint1 = cbind(ce,ktau,dcor,dhsic,hhg1,hhg2,ind,bdtau,ball,qad,bet,mixed,subcop,mdm,codec,nns)
x11()
corrplot(cor(joint1), method = "shade", order = "hclust", col = COL2(n=200))

# normalization
normalize <- function(x){  xn = (x-x[1])/(x[11]-x[1]) }
ce1 = normalize(ce)
ktau1 = normalize(ktau)
dcor1 = normalize(dcor)
dhsic1 = normalize(dhsic)
hhg1a = normalize(hhg1)
hhg2a = normalize(hhg2)
ind1 = normalize(ind)
bdtau1 = normalize(bdtau)
ball1 = normalize(ball)
qad1 = normalize(qad)
bet1 = normalize(bet)
mixed1 = normalize(mixed)
subcop1 = normalize(subcop)
mdm1 = normalize(mdm)
codec1 = normalize(codec)
nns1 = normalize(nns)

all1 = rbind(ce1,ktau1,dcor1,dhsic1,hhg1a,hhg2a,ind1,bdtau1,ball1,qad1,bet1,mixed1,subcop1,mdm1,codec1,nns1)
x11(width = 10, height = 6)
col1 = rainbow(16)
plot(ce1, xlab = "", ylab = "", ylim = c(-2,4), col = col1[1], xaxt = "n", pch = 1); lines(ce1, col = col1[1])
points(ktau1, col = col1[2], pch = 2); lines(ktau1, col = col1[2])
points(dcor1, col = col1[3], pch = 3); lines(dcor1, col = col1[3])
points(dhsic1, col = col1[4], pch = 4); lines(dhsic1, col = col1[4])
points(hhg1a, col = col1[5], pch = 5); lines(hhg1a, col = col1[5])
points(hhg2a, col = col1[6], pch = 6); lines(hhg2a, col = col1[6])
points(ind1, col = col1[7], pch = 7); lines(ind1, col = col1[7])
points(bdtau1, col = col1[8], pch = 8); lines(bdtau1, col = col1[8])
points(ball1, col = col1[9], pch = 9); lines(ball1, col = col1[9])
points(qad1, col = col1[10], pch = 10); lines(qad1, col = col1[10])
points(bet1, col = col1[11], pch = 11); lines(bet1, col = col1[11])
points(mixed1, col = col1[12], pch = 12); lines(mixed1, col = col1[12])
points(subcop1, col = col1[13], pch = 13); lines(subcop1, col = col1[13])
points(mdm1, col = col1[14], pch = 14); lines(mdm1, col = col1[14])
points(codec1, col = col1[15], pch = 15); lines(codec1, col = col1[15])
points(nns1, col = col1[16], pch = 16); lines(nns1, col = col1[16])
legend(0.65,4,legend = c("CE","Ktau","dCor","dHSIC","HHG.chisq","HHG.lr","Hoeff","BDtau","Ball","QAD","BET","mixed","subcop","MDM","CODEC","NNS"), 
       col = col1, pch = 1:16, horiz = T, cex = 0.46, box.col = "white")
axis(1, at = 1:11, labels = FALSE) # colnames(wine1)[1:11], cex = 0.46)
text(1:11+0.15, par("usr")[3] - 0.18, labels = colnames(wine1)[1:11], srt = 45, pos = 2, xpd = TRUE, cex = 0.75)

#### plot
# ce
x11(width = 5, height = 5)
plot(ce, xlab = "rho", ylab = "Copula Entropy")
lines(ce)

# kendall's tau
x11(width = 5, height = 5)
plot(ktau, xlab = "rho", ylab = "Kendall's tau", xaxt = 'n')
lines(ktau)

# dcor
x11(width = 5, height = 5)
plot(dcor, xlab = "rho", ylab = "dCor", xaxt = 'n')
lines(dcor)

# dhsic
x11(width = 5, height = 5)
plot(dhsic, xlab = "rho", ylab = "dHSIC", xaxt = 'n')
lines(dhsic)

# hhg
x11(width = 5, height = 5)
plot(hhg1, xlab = "rho", ylab = "HHG.chisquared", xaxt = 'n')
lines(hhg1)

x11(width = 5, height = 5)
plot(hhg2, xlab = "rho", ylab = "HHG.lr", xaxt = 'n')
lines(hhg2)

# independence
x11(width = 5, height = 5)
plot(ind, xlab = "rho", ylab = "Hoeffding", xaxt = 'n')
lines(ind)

# BD tau
x11(width = 5, height = 5)
plot(bdtau, xlab = "rho", ylab = "BD tau", xaxt = 'n')
lines(bdtau)

# ball
x11(width = 5, height = 5)
plot(ball, xlab = "rho", ylab = "Ball", xaxt = 'n')
lines(ball)

# qad
x11(width = 5, height = 5)
plot(qad, xlab = "rho", ylab = "qad", xaxt = 'n')
lines(qad)

# BET
x11(width = 5, height = 5)
plot(bet, xlab = "rho", ylab = "BET", xaxt = 'n')
lines(bet)

# MixedIndCopula
x11(width = 5, height = 5)
plot(mixed, xlab = "rho", ylab = "Mixed", xaxt = 'n')
lines(mixed)

# NNS
x11(width = 5, height = 5)
plot(nns, xlab = "rho", ylab = "NNS", xaxt = 'n')
lines(nns)

# subcop
x11(width = 5, height = 5)
plot(subcop, xlab = "rho", ylab = "subcop", xaxt = 'n')
lines(subcop)

# mdm
x11(width = 5, height = 5)
plot(mdm, xlab = "rho", ylab = "MDC", xaxt = 'n')
lines(mdm)

# CODEC
x11(width = 5, height = 5)
plot(codec, xlab = "rho", ylab = "CODEC", xaxt = 'n')
lines(codec)
