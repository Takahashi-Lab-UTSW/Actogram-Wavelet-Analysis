#!/usr/bin/env ipython
'''
usage:
	PRCKernelEstimationImputation.py\
		PRC_table_wo_Binning_Npas4++.csv\
		PRC_table_wo_Binning_Npas4--.csv\
		PRC_table_Binning_Npas4++.csv\
		PRC_table_Binning_Npas4--.csv\
		1 0.01 3 1.5
		1: bandwidth (h)
		0.01: smoothed curve step size (xstep)
		3: binwidth, grouping CT levels every 3 hours
		1.5: binoffset, 24 - 1.5 = 22.5 h
'''

import sys
import numpy as np
import scipy.stats
from operator import itemgetter
from pyqt_fit import kde
import scipy.integrate as integrate


# Kernel function (Gaussian)
def funcK(x):
	return 1/(np.sqrt(2*np.pi))*np.exp(-0.5*(x**2))
	#return 0.75*(1-x**2)*I

# Rosenblatt-Parzen kernel density estimator
def funch(x,xi,h,N):
	return (1/(N*h))*np.sum(funcK((x-xi)/h))

# Pointwise approximate confidence interval
def funcCI(xi, yi, h, L2K):
	N = len(xi)

	# addtional term to prevent division by zero
	dn=0.00000000001
	alpha = 0.05
	ppf = scipy.stats.norm.ppf(1-alpha/2)
	#ppf(q, loc=0, scale=1) Percent point function (inverse of cdf — percentiles).

	cnti=0

	for i in np.arange(-24,48,xstep):
		YKE = np.sum(funcK((i-xi)/h)*yi)/(np.sum(funcK((i-xi)/h))+dn)

		Wh = (funcK((i-xi)/h)/h)/(funch(i,xi,h,N)+dn)
		varest = (1/N)*np.sum(Wh*(yi-YKE)**2)

		CI_lb = YKE - ppf*np.sqrt(varest*L2K/(N*h*funch(i,xi,h,N)+dn))
		CI_ub = YKE + ppf*np.sqrt(varest*L2K/(N*h*funch(i,xi,h,N)+dn))

		CI = np.hstack([CI_lb, CI_ub])
		if cnti==0:
			CIAll = CI
		else:
			CIAll = np.vstack([CIAll, CI])
		cnti+=1
	return CIAll


WT_PRC = np.loadtxt(sys.argv[1])
KO_PRC = np.loadtxt(sys.argv[2])

WT_PRC_Binned = np.loadtxt(sys.argv[3])
KO_PRC_Binned = np.loadtxt(sys.argv[4])

h = float(sys.argv[5])
xstep = float(sys.argv[6])
Binwidth = float(sys.argv[7])
Binoffset = float(sys.argv[8])

WT_PRC = np.asarray(sorted(WT_PRC, key=itemgetter(0)))
KO_PRC = np.asarray(sorted(KO_PRC, key=itemgetter(0)))

XWT = WT_PRC[:,0]
YWT = WT_PRC[:,1]
XKO = KO_PRC[:,0]
YKO = KO_PRC[:,1]

XWT2xb = XWT - 24
XWT2xa = XWT + 24
XWT3x = np.hstack([XWT2xb, XWT, XWT2xa])
YWT3x = np.hstack([YWT, YWT, YWT])

XKO2xb = XKO - 24
XKO2xa = XKO + 24
XKO3x = np.hstack([XKO2xb, XKO, XKO2xa])
YKO3x = np.hstack([YKO, YKO, YKO])

# addtional term to prevent division by zero
dn=0.00000000001

x = np.arange(-24,48,xstep)
cnti=0

for i in np.arange(-24,48,xstep):
	YKEWT = np.sum(funcK((i-XWT3x)/h)*YWT3x)/(np.sum(funcK((i-XWT3x)/h))+dn)
	if cnti==0:
		YKEWTAll = YKEWT
	else:
		YKEWTAll = np.hstack([YKEWTAll, YKEWT])
	cnti+=1

cnti=0
for i in np.arange(-24,48,xstep):
	YKEKO = np.sum(funcK((i-XKO3x)/h)*YKO3x)/(np.sum(funcK((i-XKO3x)/h))+dn)
	if cnti==0:
		YKEKOAll = YKEKO
	else:
		YKEKOAll = np.hstack([YKEKOAll, YKEKO])
	cnti+=1

L2K = integrate.quad(lambda x: funcK(x)**2, -np.inf, np.inf)
CIWT = funcCI(XWT3x,YWT3x,h,L2K[0])
CIKO = funcCI(XKO3x,YKO3x,h,L2K[0])


#================================================================================================================
# Missing values imputation
Tcycle = 24.
NumGroups = int(Tcycle/Binwidth)

DiffNumRowsWTKO = WT_PRC_Binned.shape[0] - KO_PRC_Binned.shape[0]
print ('DiffNumRowsWTKO',DiffNumRowsWTKO)

if (DiffNumRowsWTKO > 0):
	WT_PRC_Binned_Complete = np.copy(WT_PRC_Binned)
else:
	WT_PRC_Binned_Complete = np.vstack([np.copy(WT_PRC_Binned), np.empty([DiffNumRowsWTKO, NumGroups])])
	WT_PRC_Binned_Complete[WT_PRC_Binned.shape[0]:,:] = np.NaN

for i in range(NumGroups):
	xgsWT = i*Binwidth - Binoffset	
	xgeWT = (i+1)*Binwidth - Binoffset

	CIWTg = CIWT[((x >= xgsWT) & (x < xgeWT))]
	YgWT = YKEWTAll[((x >= xgsWT) & (x < xgeWT))]

	YgWTMax = np.max(YgWT)
	YgWTMin = np.min(YgWT)

	CIWTgDiffMax = np.max(CIWTg[:,1] - YgWT)

	YgImpMax = YgWTMax + CIWTgDiffMax
	YgImpMin = YgWTMin - CIWTgDiffMax

	WT_PRC_Binned_group = WT_PRC_Binned[:,i]
	nanlist = np.where(np.isnan(WT_PRC_Binned_group))[0]

	if (len(nanlist) > 0):
		for j in range(len(nanlist)):
			WT_PRC_Binned_Complete[j + nanlist[0],i] = np.random.uniform(YgImpMin, YgImpMax)

if (DiffNumRowsWTKO > 0):
	KO_PRC_Binned_Complete = np.vstack([np.copy(KO_PRC_Binned), np.empty([DiffNumRowsWTKO, NumGroups])])
	KO_PRC_Binned_Complete[KO_PRC_Binned.shape[0]:,:] = np.NaN
else:
	KO_PRC_Binned_Complete = np.copy(KO_PRC_Binned)

for i in range(NumGroups):
	xgsKO = i*Binwidth - Binoffset	
	xgeKO = (i+1)*Binwidth - Binoffset

	CIKOg = CIKO[((x >= xgsKO) & (x < xgeKO))]
	YgKO = YKEKOAll[((x >= xgsKO) & (x < xgeKO))]

	YgKOMax = np.max(YgKO)
	YgKOMin = np.min(YgKO)

	CIKOgDiffMax = np.max(CIKOg[:,1] - YgKO)

	YgImpMax = YgKOMax + CIKOgDiffMax
	YgImpMin = YgKOMin - CIKOgDiffMax

	KO_PRC_Binned_group = KO_PRC_Binned_Complete[:,i]
	nanlist = np.where(np.isnan(KO_PRC_Binned_group))[0]

	if (len(nanlist) > 0):
		for j in range(len(nanlist)):
			KO_PRC_Binned_Complete[j + nanlist[0]] = np.random.uniform(YgImpMin, YgImpMax)

np.savetxt('PRC_table_Binning_Npas4++_Imputed.csv', WT_PRC_Binned_Complete, fmt='%10.5f', delimiter="\t")
np.savetxt('PRC_table_Binning_Npas4--_Imputed.csv', KO_PRC_Binned_Complete, fmt='%10.5f', delimiter="\t")
	

