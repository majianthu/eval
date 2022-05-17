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

scan_heart_data <-function(filename1, nl = 0){
  data1 = scan(filename1, nlines = nl, what = c(as.list(rep(0,75)),list("")))
  l = length(data1[[1]])
  data1m = matrix(unlist(data1), l, 76)
  matrix(as.numeric(data1m[,1:75]), l, 75)
}
#### load heart disease data (899 samples)
h1 = scan_heart_data("http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/cleveland.data", 282*10)
h2 = scan_heart_data("http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/hungarian.data")
h3 = scan_heart_data("http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/switzerland.data")
h4 = scan_heart_data("http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/long-beach-va.data")

heart1 = as.matrix( rbind(h1,h2,h3,h4) )
m = dim(heart1)[1]
n = dim(heart1)[2]

## statistical dependence with attr #58
# Copula Entropy
l = 30
ce = rep(0,n)
for (i in 1:n){
  for (j in 1:l){
    data2 = heart1[,c(i,58)]
    data2[,1] = data2[,1] + max(abs(data2[,1])) * 0.000001 * rnorm(m)
    data2[,2] = data2[,2] + max(abs(data2[,2])) * 0.000001 * rnorm(m)
    ce[i] = ce[i] + copent(data2)/l
  }
}
ce[c(1,2,58)] = min(ce)

# Other Independence Measures
ktau = rep(0,n) # Kendall's tau
dcor = rep(0,n) # Distance Correlation
dhsic = rep(0,n)  # Hilbert-Schmidt Independence Criterion
hhg1 = rep(0,n)  # Heller-Heller-Gorfine Tests
hhg2 = rep(0,n)  # Heller-Heller-Gorfine Tests
hoeff = rep(0,n)  # Hoeffding's D test or 
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
  ktau[i] = cor(heart1[,c(i,58)], method = "kendall")[1,2]
  dcor[i] = dcor(heart1[,i],heart1[,58])
  dhsic[i] = dhsic(heart1[,i],heart1[,58])$dHSIC
  Dx = as.matrix(dist((heart1[,i]), diag = TRUE, upper = TRUE))
  Dy = as.matrix(dist((heart1[,58]), diag = TRUE, upper = TRUE))
  hhg = hhg.test(Dx,Dy, nr.perm = 500)
  hhg1[i] = hhg$sum.chisq
  hhg2[i] = hhg$sum.lr
  hoeff[i] = hoeffding.D.test(heart1[,i],heart1[,58])$Dn
  bdtau[i] = tau.star.test(heart1[,i],heart1[,58])$Tn
  ball[i] = bcor(heart1[,i],heart1[,58])
  qad[i] = qad(heart1[,i],heart1[,58])$`q(X,Y)`
  bet[i] = MaxBET(heart1[,c(i,58)], 3, index = list(1,2))$z.statistic
  mixed[i] = TestIndCopula(heart1[,c(i,58)])$stat$cvm
  nns[i] = NNS.dep(heart1[,i],heart1[,58])$Dependence
  subcop[i] = dependence(heart1[,c(i,58)])[1,2,2]
  mdm[i] = mdc(heart1[,i],heart1[,58])
  codec[i] = codec(heart1[,i],heart1[,58])
}
ktau[c(1,2,58)] = 0
dcor[c(1,2,58)] = 0
dhsic[c(1,2,58)] = 0
hhg1[c(1,2,58)] = 0
hhg2[c(1,2,58)] = 0
hoeff[c(1,2,58)] = 0
bdtau[c(1,2,58)] = 0
ball[c(1,2,58)] = 0
qad[c(1,2,58)] = 0
bet[c(1,2,58)] = 0
mixed[c(1,2,58)] = 0
nns[c(1,2,58)] = 0
subcop[c(1,2,58)] = 0
mdm[c(1,2,58)] = 0
codec[c(1,2,58)] = 0

# joint
joint1 = cbind(ce,ktau,dcor,dhsic,hhg1,hhg2,hoeff,bdtau,ball,qad,bet,mixed,subcop,mdm,codec,nns)
x11()
corrplot(cor(joint1), method = "shade", order = "hclust", col = COL2(n=200))

#### plot
w1 = 8; h1 = 5
# ce
x11(width = w1, height = h1)
plot(ce, xlab = "Variable", ylab = "CE", xaxt = 'n')
lines(ce)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(ce[16],75)
lines(th16, col = "red")


# doc
x11(width = w1, height = h1)
plot(dcor, xlab = "Variable", ylab = "dCor", xaxt = 'n')
lines(dcor)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(dcor[16],75)
lines(th16, col = "red")

# dhsic
x11(width = w1, height = h1)
plot(dhsic, xlab = "Variable", ylab = "dHSIC", xaxt = 'n')
lines(dhsic)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(dhsic[16],75)
lines(th16, col = "red")

# ball
x11(width = w1, height = h1)
plot(ball, xlab = "Variable", ylab = "Ball", xaxt = 'n')
lines(ball)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(ball[16],75)
lines(th16, col = "red")

# bdtau
x11(width = w1, height = h1)
plot(bdtau, xlab = "Variable", ylab = "BDtau", xaxt = 'n')
lines(bdtau)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(bdtau[16],75)
lines(th16, col = "red")

# bet
x11(width = w1, height = h1)
plot(bet, xlab = "Variable", ylab = "BET", xaxt = 'n')
lines(bet)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(bet[16],75)
lines(th16, col = "red")

# codec
x11(width = w1, height = h1)
plot(codec, xlab = "Variable", ylab = "CODEC", xaxt = 'n')
lines(codec)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(codec[16],75)
lines(th16, col = "red")

# hhg1
x11(width = w1, height = h1)
plot(hhg1, xlab = "Variable", ylab = "HHG.chisq", xaxt = 'n')
lines(hhg1)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(hhg1[16],75)
lines(th16, col = "red")

# hhg2
x11(width = w1, height = h1)
plot(hhg2, xlab = "Variable", ylab = "HHG.lr", xaxt = 'n')
lines(hhg2)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(hhg2[16],75)
lines(th16, col = "red")

# hoeff
x11(width = w1, height = h1)
plot(hoeff, xlab = "Variable", ylab = "Hoeff", xaxt = 'n')
lines(hoeff)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(hoeff[16],75)
lines(th16, col = "red")

# ktau
x11(width = w1, height = h1)
plot(ktau, xlab = "Variable", ylab = "Ktau", xaxt = 'n')
lines(ktau)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(ktau[16],75)
lines(th16, col = "red")

# mdm
x11(width = w1, height = h1)
plot(mdm, xlab = "Variable", ylab = "mdm", xaxt = 'n')
lines(mdm)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(mdm[16],75)
lines(th16, col = "red")

# mixed
x11(width = w1, height = h1)
plot(mixed, xlab = "Variable", ylab = "mixed", xaxt = 'n')
lines(mixed)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(mixed[16],75)
lines(th16, col = "red")

# nns
x11(width = w1, height = h1)
plot(nns, xlab = "Variable", ylab = "NNS", xaxt = 'n')
lines(nns)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(nns[16],75)
lines(th16, col = "red")

# qad
x11(width = w1, height = h1)
plot(qad, xlab = "Variable", ylab = "qad", xaxt = 'n')
lines(qad)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(qad[16],75)
lines(th16, col = "red")

# subcop
x11(width = w1, height = h1)
plot(subcop, xlab = "Variable", ylab = "subcop", xaxt = 'n')
lines(subcop)
axis(side = 1, at = c(1,seq(15,75, by = 15)), labels = c(1,seq(15,75, by = 15)))
th16 = rep(subcop[16],75)
lines(th16, col = "red")
