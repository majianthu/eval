library(copula)
library(corrplot)
library(latex2exp)
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

ce = kci = rcot = cdc = codec = gcm = wgcm = kpc = pcor = ckt = cmd = pcop = 0

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
py1 = read.csv("~/Rworks/bench/py1a.csv") # simulation 1 - normal distribution
#py1 = read.csv("~/Rworks/bench/py1b.csv") # simulation 2 - copula function
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
rho1 = seq(0,0.9,0.1); xlab1 = TeX(r'($\rho_{xz}$)')
plot(rho1,ce, xlab = xlab1, ylab = "stats", main = "CE");lines(rho1,ce)
plot(rho1,kci, xlab = xlab1, ylab = "stats", main = "KCI");lines(rho1,kci)
plot(rho1,rcot, xlab = xlab1, ylab = "stats", main = "RCoT");lines(rho1,rcot)
plot(rho1,cdc, xlab = xlab1, ylab = "stats", main = "CDC");lines(rho1,cdc)
plot(rho1,codec, xlab = xlab1, ylab = "stats", main = "CODEC");lines(rho1,codec)
plot(rho1,gcm, xlab = xlab1, ylab = "stats", main = "GCM");lines(rho1,gcm)
plot(rho1,wgcm, xlab = xlab1, ylab = "stats", main = "wGCM");lines(rho1,wgcm)
plot(rho1,kpc, xlab = xlab1, ylab = "stats", main = "KPC");lines(rho1,kpc)
plot(rho1,pcor, xlab = xlab1, ylab = "stats", main = "Partial Correlation");lines(rho1,pcor)
plot(rho1,cmd, xlab = xlab1, ylab = "stats", main = "CMD");lines(rho1,cmd)
plot(rho1,pcop, xlab = xlab1, ylab = "stats", main = "PartialCopula");lines(rho1,pcop)
plot(rho1,cmi1, xlab = xlab1, ylab = "stats", main = "CMI1");lines(rho1,cmi1)
plot(rho1,cmi2, xlab = xlab1, ylab = "stats", main = "CMI2");lines(rho1,cmi2)
plot(rho1,ccit, xlab = xlab1, ylab = "stats", main = "CCIT");lines(rho1,ccit)
plot(rho1,fcit, xlab = xlab1, ylab = "stats", main = "FCIT");lines(rho1,fcit)
plot(rho1,pcit, xlab = xlab1, ylab = "stats", main = "PCIT");lines(rho1,pcit)


