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
# py1 = read.csv("~/Rworks/bench/py1a.csv") # simulation 1 - normal
py1 = read.csv("~/Rworks/bench/py1b.csv") # simulation 2 - copula
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
x11(width = 10, height = 12)
par(mfrow = c(4,4))
plot(ce, xlab = "rho", ylab = "stats", main = "CE");lines(ce)
plot(kci, xlab = "rho", ylab = "stats", main = "KCI");lines(kci)
plot(rcot, xlab = "rho", ylab = "stats", main = "RCoT");lines(rcot)
plot(cdc, xlab = "rho", ylab = "stats", main = "CDC");lines(cdc)
plot(codec, xlab = "rho", ylab = "stats", main = "CODEC");lines(codec)
plot(gcm, xlab = "rho", ylab = "stats", main = "GCM");lines(gcm)
plot(wgcm, xlab = "rho", ylab = "stats", main = "wGCM");lines(wgcm)
plot(kpc, xlab = "rho", ylab = "stats", main = "KPC");lines(kpc)
plot(pcor, xlab = "rho", ylab = "stats", main = "Partial Correlation");lines(pcor)
plot(cmd, xlab = "rho", ylab = "stats", main = "CMD");lines(cmd)
plot(pcop, xlab = "rho", ylab = "stats", main = "PartialCopula");lines(pcop)
plot(cmi1, xlab = "rho", ylab = "stats", main = "CMI1");lines(cmi1)
plot(cmi2, xlab = "rho", ylab = "stats", main = "CMI2");lines(cmi2)
plot(ccit, xlab = "rho", ylab = "stats", main = "CCIT");lines(ccit)
plot(fcit, xlab = "rho", ylab = "stats", main = "FCIT");lines(fcit)
plot(pcit, xlab = "rho", ylab = "stats", main = "PCIT");lines(pcit)
