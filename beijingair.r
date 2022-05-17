library(corrplot)
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

prsa2010data = read.csv("~/Rworks/beijingair/PRSA_data_2010.1.1-2014.12.31.csv")
# id: 6(PM2.5), 7(Dew Point), 8(Temperature), 9(Pressure), 11(Cumulative Wind Speed)
idx = c(6,9)
data = prsa2010data[2200:2700, idx]

ce = 0
# kci = 0
rcot = 0
cdc = 0
codec = 0
gcm = 0
wgcm = 0
kpc = 0
pcor = 0
cmd = 0
pcop = 0

for (lag in 1:24){
  pm25a = data[1:(501-lag),1]
  pm25b = data[(lag+1):501,1]
  v1 = data[1:(501-lag),2]
  
  ce[lag] = ci(pm25b,v1,pm25a)
  # kci[lag] = KCI(pm25b,v1,pm25a)$testStatistic
  rcot[lag] = RCIT(pm25b,v1,pm25a)$Sta
  cdc[lag] = cdcor(pm25b,v1,pm25a)$statistic
  codec[lag] = codec(pm25b,v1,pm25a)
  gcm[lag] = gcm.test(pm25b,v1,pm25a)$test.statistic
  wgcm[lag] = wgcm.est(pm25b,v1,pm25a, regr.meth = "xgboost")
  kpc[lag] = KPCRKHS(pm25b,pm25a,v1)
  pcor[lag] = pcor.test(pm25b,v1,pm25a)$statistic
  cmd[lag] = cmdm_test(pm25b,v1,pm25a)$stat
  pcop[lag] = test_CI(pm25b,v1,pm25a)$statistic
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
plot(ce, xlab = "lag (hours)", ylab = "Transfer Entropy", main = "Pressure")
lines(ce)

# KCI
x11(width = w1, height = h1)
plot(kci, xlab = "lag (hours)", ylab = "KCI", main = "Pressure")
lines(kci)

# RCoT
x11(width = w1, height = h1)
plot(rcot, xlab = "lag (hours)", ylab = "RCoT", main = "Pressure")
lines(rcot)

# CDC
x11(width = w1, height = h1)
plot(cdc, xlab = "lag (hours)", ylab = "CDC", main = "Pressure")
lines(cdc)

# CODEC
x11(width = w1, height = h1)
plot(codec, xlab = "lag (hours)", ylab = "CODEC", main = "Pressure")
lines(codec)

# GCM
x11(width = w1, height = h1)
plot(gcm, xlab = "lag (hours)", ylab = "GCM", main = "Pressure")
lines(gcm)

# wGCM
x11(width = w1, height = h1)
plot(wgcm, xlab = "lag (hours)", ylab = "wGCM", main = "Pressure")
lines(wgcm)

# KPC
x11(width = w1, height = h1)
plot(kpc, xlab = "lag (hours)", ylab = "KPC", main = "Pressure")
lines(kpc)

# Partial Correlation
x11(width = w1, height = h1)
plot(pcor, xlab = "lag (hours)", ylab = "Partial Correlation", main = "Pressure")
lines(pcor)

# Conditional Mean Dependence
x11(width = w1, height = h1)
plot(cmd, xlab = "lag (hours)", ylab = "CMD", main = "Pressure")
lines(cmd)

# partial copula based CI test
x11(width = w1, height = h1)
plot(pcop, xlab = "lag (hours)", ylab = "PartialCopula", main = "Pressure")
lines(pcop)

# cmi1
x11(width = w1, height = h1)
plot(cmi1, xlab = "lag (hours)", ylab = "CMI1", main = "Pressure")
lines(cmi1)

# cmi2
x11(width = w1, height = h1)
plot(cmi2, xlab = "lag (hours)", ylab = "CMI2", main = "Pressure")
lines(cmi2)

# ccit
x11(width = w1, height = h1)
plot(ccit, xlab = "lag (hours)", ylab = "CCIT", main = "Pressure")
lines(ccit)

# fcit
x11(width = w1, height = h1)
plot(fcit, xlab = "lag (hours)", ylab = "FCIT", main = "Pressure")
lines(fcit)

# pcit
x11(width = w1, height = h1)
plot(pcit, xlab = "lag (hours)", ylab = "PCIT", main = "Pressure")
lines(pcit)

