library(corrplot)
library(mnormt)
library(copent) # transfer entropy (TE)
library(CondIndTests) # kernel-based CI test (KCI)
library(RCIT) # Randomized conditional Correlation Test (RCoT)
library(cdcsis) # conditional distance correlation (CDC)
library(FOCI) # conditional dependence coefficient (CODEC)
library(GeneralisedCovarianceMeasure) # Generalised Covariance Measure (GCM)
library(weightedGCM) # weighted GCM
library(KPC) # Kernel Partial Correlation (KPC)
library(ppcor) # Partial Correlation (pcor)
library(CondCopulas) # Conditional Kendall's Tau (CKT)
library(EDMeasure) # Conditional Mean Dependence (CMD)
source("https://raw.githubusercontent.com/lassepetersen/partial-copula-CI-test/main/parCopCITest.R") # partial copula based CI test

ce = 0
kci = 0
rcot = 0
cdc = 0
codec = 0
gcm = 0
wgcm = 0
kpc = 0
pcor = 0
ckt = 0
cmd = 0
pcop = 0

for (i in 1:10){
  rxy = 0.7
  ryz = 0.6
  rxz = (i-1)/10
  # normal
  sigma = matrix(c(1,rxy,rxz,rxy,1,ryz,rxz,ryz,1),3,3)
  xyz = rmnorm(800,c(0,0,0),sigma)
  # copula
  # norm.cop <- normalCopula(param = c(rxy, ryz, rxz), dim = 3, dispstr = "un")
  # mv.NE <- mvdc(norm.cop, c("norm", "exp","exp"), list(list(mean = 0, sd =2), list(rate = 2), list(rate = .5)))
  # xyz <- rMvdc(800, mv.NE)
  
  x = xyz[,1]; y = xyz[,2]; z = xyz[,3]
  
  ce[i] = ci(x,y,z)
  kci[i] = KCI(x,y,z)$testStatistic
  rcot[i] = RCIT(x,y,z)$Sta
  cdc[i] = cdcor(x,y,z)$statistic
  codec[i] = codec(x,y,z)
  gcm[i] = gcm.test(x,y,z)$test.statistic
  wgcm[i] = wgcm.est(x,y,z, regr.meth = "xgboost", beta = 0.7)
  kpc[i] = KPCRKHS(x,z,y)
  pcor[i] = pcor.test(x,y,z)$statistic
  cmd[i] = cmdm_test(x,y,z, compute = "R")$stat
  pcop[i] = test_CI(x,y,z)$statistic
}

# python results
py1 = read.csv("~/Rworks/bench/py1.csv")
kci = py1$kci
cmi1 = py1$knn
cmi2 = py1$cmi
fcit = py1$fcit
ccit = py1$ccit
pcit = py1$pcit

# joint
joint1 = cbind(ce,kci,rcot,cdc,gcm,wgcm,codec,kpc,pcor,cmd,pcop,cmi1,cmi2,fcit,ccit,pcit)
x11()
corrplot(cor(joint1), method = "shade", order = "hclust", col = COL2(n=200))

# plotting
w1 = 5; h1 = 5

# TE via CE
x11(width = w1, height = h1)
plot(ce, xlab = "rho", ylab = "CI")
lines(ce)

# KCI
x11(width = w1, height = h1)
plot(kci, xlab = "rho", ylab = "KCI")
lines(kci)

# RCoT
x11(width = w1, height = h1)
plot(rcot, xlab = "rho", ylab = "RCoT")
lines(rcot)

# CDC
x11(width = w1, height = h1)
plot(cdc, xlab = "rho", ylab = "CDC")
lines(cdc)

# CODEC
x11(width = w1, height = h1)
plot(codec, xlab = "rho", ylab = "CODEC")
lines(codec)

# GCM
x11(width = w1, height = h1)
plot(gcm, xlab = "rho", ylab = "GCM")
lines(gcm)

# wGCM
x11(width = w1, height = h1)
plot(wgcm, xlab = "rho", ylab = "wGCM")
lines(wgcm)

# KPC
x11(width = w1, height = h1)
plot(kpc, xlab = "rho", ylab = "KPC")
lines(kpc)

# Partial Correlation
x11(width = w1, height = h1)
plot(pcor, xlab = "rho", ylab = "Partial Correlation")
lines(pcor)

# Conditional Mean Dependence
x11(width = w1, height = h1)
plot(cmd, xlab = "rho", ylab = "CMD")
lines(cmd)

# partial copula based CI test
x11(width = w1, height = h1)
plot(pcop, xlab = "rho", ylab = "PartialCopula")
lines(pcop)

# cmi1
# x11(width = w1, height = h1)
plot(cmi1, xlab = "rho", ylab = "CMI1")
lines(cmi1)

# cmi2
x11(width = w1, height = h1)
plot(cmi2, xlab = "rho", ylab = "CMI2")
lines(cmi2)

# ccit
x11(width = w1, height = h1)
plot(ccit, xlab = "rho", ylab = "CCIT")
lines(ccit)

# fcit
x11(width = w1, height = h1)
plot(fcit, xlab = "rho", ylab = "FCIT")
lines(fcit)

# pcit
x11(width = w1, height = h1)
plot(pcit, xlab = "rho", ylab = "PCIT")
lines(pcit)

