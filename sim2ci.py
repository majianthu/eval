from copent import ci
from knncmi import cmi
from causallearn.utils.KCI.KCI import KCI_CInd
from pycit.estimators import ksg_cmi
from fcit import fcit
from CCIT import CCIT
from pcit.IndependenceTest import PCIT  ##  the 'mlxtend' package required
from numpy.random import multivariate_normal as mnorm
from pycop import simulation
from scipy.stats import norm, expon

import pandas as pd
import numpy as np

te1 = np.zeros(10)
cmi1 = np.zeros(10)
kci1 = np.zeros(10)
ci1 = np.zeros(10)
fcit1 = np.zeros(10)
ccit1 = np.zeros(10)
pcit1 = np.zeros(10)
for i in range(0,2):
	rxy = 0.7
	ryz = 0.6
	rxz = i/10
	m1 = [0,0,0]
	sigma1 = [ [1,rxy,rxz],[rxy,1,ryz],[rxz,ryz,1] ]

	# normal
	xyz = mnorm(m1, sigma1, 800) # tri-variate gaussian 
	x = xyz[:,0]
	y = xyz[:,1]
	z = xyz[:,2]
	
	# copula
	# ncop1 = simulation.simu_gaussian(3,800,sigma1).T
	# x = norm.ppf(ncop1[:,0], loc = 0, scale = 2)
	# y = expon.ppf(ncop1[:,1], scale = 0.5)
	# z = expon.ppf(ncop1[:,2], scale = 2)
	
	## copent
	te1[i] = ci(x,y,z)
	## knncmi
	df = pd.DataFrame(np.vstack((x,y,z)).T, columns = ['x','y','z'])
	cmi1[i] = cmi(['x'], ['y'], ['z'], 3, df)	
	## pycit
	ci1[i] = ksg_cmi(x, y, z)
	## fcit
	len1 = len(x)
	xa = np.reshape(x,[len1,1])
	ya = np.reshape(y,[len1,1])
	za = np.reshape(z,[len1,1])
	fcit1[i] = fcit.test(xa,ya,za)
	## causal-learn
	kci_cind = KCI_CInd()
	_,kci1[i] = kci_cind.compute_pvalue(xa,ya,za)
	## CCIT
	ccit1[i] = CCIT.CCIT(xa,ya,za)
	## pcit
	p,tmp,tmp = PCIT(xa,ya,z = za)
	pcit1[i] = p[0]
	
	str = "i: %d; TE:%.3f; CMI:%.3f; KCI:%.3f; CI:%.3f; fcit:%.3f; ccit:%.3f; pcit:%.3f" %(i,te1[i],cmi1[i],kci1[i],ci1[i],fcit1[i],ccit1[i],pcit1[i])
	print(str)

data1 = np.vstack((cmi1,kci1,ci1,fcit1,ccit1,pcit1))
df1 = pd.DataFrame(data1.T)
df1.columns = ["cmi","kci","knn","fcit","ccit","pcit"]
df1.to_csv("~/Rworks/bench/py1.csv")

