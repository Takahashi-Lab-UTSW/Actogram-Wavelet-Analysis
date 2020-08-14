#!/usr/bin/env ipython
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import pandas as pd

A2KO_PhaseInfo = np.loadtxt(sys.argv[1])
A2WT_PhaseInfo = np.loadtxt(sys.argv[2])
A3KO_PhaseInfo = np.loadtxt(sys.argv[3])
A4WT_PhaseInfo = np.loadtxt(sys.argv[4])
C2KO_PhaseInfo = np.loadtxt(sys.argv[5])
C2WT_PhaseInfo = np.loadtxt(sys.argv[6])

Rsqthreshold = float(sys.argv[7])
RMSEthreshold = float(sys.argv[8])

nf = int(sys.argv[9])
# nf = 14 (number of columns (phase information) in each file)

BinWidth = float(sys.argv[10])
BinOffset = float(sys.argv[11])

NumFiles_A2KO = int(A2KO_PhaseInfo.shape[1]/nf)
NumFiles_A2WT = int(A2WT_PhaseInfo.shape[1]/nf)
NumFiles_A3KO = int(A3KO_PhaseInfo.shape[1]/nf)
NumFiles_A4WT = int(A4WT_PhaseInfo.shape[1]/nf)
NumFiles_C2KO = int(C2KO_PhaseInfo.shape[1]/nf)
NumFiles_C2WT = int(C2WT_PhaseInfo.shape[1]/nf)

NumFilesAll = ([NumFiles_A2KO, NumFiles_A3KO, NumFiles_C2KO],
			[NumFiles_A2WT, NumFiles_A4WT, NumFiles_C2WT])

NumLightPulses_A2KO = int(A2KO_PhaseInfo.shape[0])
NumLightPulses_A2WT = int(A2WT_PhaseInfo.shape[0])
NumLightPulses_A3KO = int(A3KO_PhaseInfo.shape[0])
NumLightPulses_A4WT = int(A4WT_PhaseInfo.shape[0])
NumLightPulses_C2KO = int(C2KO_PhaseInfo.shape[0])
NumLightPulses_C2WT = int(C2WT_PhaseInfo.shape[0])

NumLPAll = ([NumLightPulses_A2KO, NumLightPulses_A3KO, NumLightPulses_C2KO],
			[NumLightPulses_A2WT, NumLightPulses_A4WT, NumLightPulses_C2WT])

WTtuple = (A2WT_PhaseInfo, A4WT_PhaseInfo, C2WT_PhaseInfo)
KOtuple = (A2KO_PhaseInfo, A3KO_PhaseInfo, C2KO_PhaseInfo)

Tcycle = 24.

cntKO=0
cntWT=0

for i in range(2):

	cntidx=0
	if (i==0):
		WTKO = KOtuple 
	else:
		WTKO = WTtuple 
	cnti=0
	for j in range(3):
		#     j=0    j=1       j=2
		#i=0: A2_KO, A3_KO, or C2_KO
		#i=1: A2_WT, A4_WT, or C2_WT

		numLP = NumLPAll[i][j]
	
		for l in range(NumFilesAll[i][j]):

			ck=0
			for k in range(numLP):
				X = WTKO[j][k,nf*l]		# phase (CT)
				Y = WTKO[j][k,nf*l+1]	# phase shift
					
				WTKO_R2b = WTKO[j][k,nf*l+2]
				WTKO_R2a = WTKO[j][k,nf*l+3]
				WTKO_R2_ba = np.hstack([WTKO_R2b, WTKO_R2a, k+1])

				WTKO_SSEb = WTKO[j][k,nf*l+8]
				WTKO_SSEa = WTKO[j][k,nf*l+9]
				WTKO_SSE_ba = np.hstack([WTKO_SSEb, WTKO_SSEa, k+1])

				WTKO_RMSEb = WTKO[j][k,nf*l+12]
				WTKO_RMSEa = WTKO[j][k,nf*l+13]
				WTKO_RMSE_ba = np.hstack([WTKO_RMSEb, WTKO_RMSEa, k+1])

				WTKO_taub = WTKO[j][k,nf*l+5]
				WTKO_taua = WTKO[j][k,nf*l+7]
				WTKO_tau_ba = np.hstack([WTKO_taub, WTKO_taua, k+1])


				if ck==0:
					XY_all = np.asarray([X,Y,k+1])
					WTKO_R2_all = WTKO_R2_ba
					WTKO_SSE_all = WTKO_SSE_ba
					WTKO_RMSE_all = WTKO_RMSE_ba
					WTKO_tau_all = WTKO_tau_ba
				else:
					XY_all = np.vstack([XY_all, np.asarray([X,Y,k+1])])
					WTKO_R2_all = np.vstack([WTKO_R2_all, WTKO_R2_ba])
					WTKO_SSE_all = np.vstack([WTKO_SSE_all, WTKO_SSE_ba])
					WTKO_RMSE_all = np.vstack([WTKO_RMSE_all, WTKO_RMSE_ba])
					WTKO_tau_all = np.vstack([WTKO_tau_all, WTKO_tau_ba])

				ck+=1

			if ck>0:			
				if ck==1:
					XY_all = XY_all.reshape((1,3))
					WTKO_R2_all = WTKO_R2_all.reshape((1,2))
					WTKO_SSE_all = WTKO_SSE_all.reshape((1,2))
					WTKO_RMSE_all = WTKO_RMSE_all.reshape((1,2))
					WTKO_tau_all = WTKO_tau_all.reshape((1,2))

				X_all = XY_all[:,0]
				Y_all = XY_all[:,1]

				idx = np.where(	(WTKO_tau_all[:,0]*WTKO_tau_all[:,1]>0) &\
								(WTKO_R2_all[:,0] >= Rsqthreshold) & (WTKO_R2_all[:,1] >= Rsqthreshold) &\
								(WTKO_RMSE_all[:,0] <= RMSEthreshold) & (WTKO_RMSE_all[:,1] <= RMSEthreshold))[0]

				XX = X_all[idx]
				YY = Y_all[idx]

				XYR = np.hstack([np.reshape(XX,[XX.shape[0],1]), np.reshape(YY,[YY.shape[0],1])])

				if i==0:	
					if cntidx==0:
						XYALLKOR = XYR
					else:
						XYALLKOR = np.vstack([XYALLKOR, XYR])

					cntKO+=idx.shape[0]
					cntidx+=1	
					
				else:
					if cntidx==0:
						XYALLWTR = XYR
					else:
						XYALLWTR = np.vstack([XYALLWTR, XYR])

					cntWT+=idx.shape[0]
					cntidx+=1	

XYALLKOWTR = np.vstack([XYALLKOR, XYALLWTR])
XYALLKOWTRLabelonly = np.hstack([np.repeat('KO', XYALLKOR.shape[0]), np.repeat('WT', XYALLWTR.shape[0])])
XYALLKOWTRwLabel = np.hstack([XYALLKOWTR, XYALLKOWTRLabelonly.reshape(XYALLKOWTRLabelonly.shape[0],1)])

np.savetxt('PRC_table_wo_Binning_Npas4++.csv', XYALLWTR, fmt='%10.5f', delimiter="\t")
np.savetxt('PRC_table_wo_Binning_Npas4--.csv', XYALLKOR, fmt='%10.5f', delimiter="\t")


data_tempR = pd.DataFrame(data=XYALLKOWTRwLabel, columns=['p', 'dp', 'WTKO'])
dataR = data_tempR[data_tempR['dp'].astype(float).notnull()]

XWTR = dataR[dataR['WTKO']=='WT']['p'].astype(float)
YWTR = dataR[dataR['WTKO']=='WT']['dp'].astype(float)
XKOR = dataR[dataR['WTKO']=='KO']['p'].astype(float)
YKOR = dataR[dataR['WTKO']=='KO']['dp'].astype(float)

#===================================================================================================
#===================================================================================================
# PRC table
GKOR = dataR[dataR['WTKO']=='KO']['WTKO']
GWTR = dataR[dataR['WTKO']=='WT']['WTKO']

XKORnew = XKOR.values.reshape(XKOR.shape[0],1)
YKORnew = YKOR.values.reshape(YKOR.shape[0],1)
GKORnew = GKOR.values.reshape(GKOR.shape[0],1)

XWTRnew = XWTR.values.reshape(XWTR.shape[0],1)
YWTRnew = YWTR.values.reshape(YWTR.shape[0],1)
GWTRnew = GWTR.values.reshape(GWTR.shape[0],1)

XYLKOR = np.hstack([XKORnew, YKORnew, GKORnew])
XYLWTR = np.hstack([XWTRnew, YWTRnew, GWTRnew])

pstart = Tcycle - BinOffset
interval = BinWidth

maxrows = 100 # number of observations in each BinWidth (arbitrary large number)
maxcols = int(Tcycle/BinWidth)

WT_group = np.empty((maxrows,maxcols))
WT_group[:] = np.nan

KO_group = np.empty((maxrows,maxcols))
KO_group[:] = np.nan

XYLWTR_group = np.empty((maxrows,maxcols))
XYLWTR_group[:] = np.nan

XYLKOR_group = np.empty((maxrows,maxcols))
XYLKOR_group[:] = np.nan

for i in range(maxcols):

	ps = pstart + i*interval
	pe = pstart + (i+1)*interval

	if (ps>Tcycle):
		ps-=Tcycle
	if (pe>Tcycle):
		pe-=Tcycle

	if i==0:
		XYLWTR_group = XYLWTR[np.where(((XYLWTR[:,0] >= ps) & (XYLWTR[:,0] < Tcycle)) | ((XYLWTR[:,0] >=0) & (XYLWTR[:,0] < pe)))[0],1]
		XYLKOR_group = XYLKOR[np.where(((XYLKOR[:,0] >= ps) & (XYLKOR[:,0] < Tcycle)) | ((XYLKOR[:,0] >=0) & (XYLKOR[:,0] < pe)))[0],1]
	else:
		XYLWTR_group = XYLWTR[np.where((XYLWTR[:,0] >= ps) & (XYLWTR[:,0] < pe))[0], 1]
		XYLKOR_group = XYLKOR[np.where((XYLKOR[:,0] >= ps) & (XYLKOR[:,0] < pe))[0], 1]

	WT_group[:XYLWTR_group.shape[0],i] = np.sort(XYLWTR_group[:])
	KO_group[:XYLKOR_group.shape[0],i] = np.sort(XYLKOR_group[:])

WT_group_woNaN = WT_group[:np.sum(np.sum(np.isnan(WT_group),axis=1)<8),:]
KO_group_woNaN = KO_group[:np.sum(np.sum(np.isnan(KO_group),axis=1)<8),:]

np.savetxt('PRC_table_Binning_Npas4++.csv', WT_group_woNaN, fmt='%10.5f', delimiter="\t")
np.savetxt('PRC_table_Binning_Npas4--.csv', KO_group_woNaN, fmt='%10.5f', delimiter="\t")



#================================================================================================================
#================================================================================================================
# plot PRC & PTC
fig = plt.figure(num=None, figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')

xlR = 0; xrR = Tcycle; ybR = -Tcycle/2.; ytR = Tcycle/2.
xlT = 0; xrT = Tcycle*2; ybT = 0; ytT = Tcycle*2

ticklength = 6
dotsize1 = 30
dotsize2 = 12

#================================================================================================================
# PRC (Npas4++)
gs= gridspec.GridSpec(2, 2)
ax1 = fig.add_subplot(gs[0,0])

ax1.scatter(XWTR, YWTR, s=dotsize1, color='black', label='WT', facecolors='black' )

ax1.axhline(0, color='k', lw=2)
ax1.axis([xlR, xrR, ybR, ytR])
ax1.set_aspect('equal')

ax1.set_xlabel('Circadian Time (h)', fontsize=25, labelpad=5)
ax1.set_ylabel('Phase shift (h)', fontsize=25, labelpad=-2)

ax1.set_xticks(np.arange(0, 24.1, 6))
ax1.set_xticklabels(['0', '6', '12', '18', '24'], fontsize=20)
ax1.set_yticks(np.arange(-12, 12.1, 6))
ax1.set_yticklabels(['-12', '-6', '0', '6', '12'], fontsize=20)

ax1.tick_params(top=True, right=True)
ax1.tick_params(axis='x', labelsize=20, direction='in', length=ticklength)
ax1.tick_params(axis='y', labelsize=20, direction='in', length=ticklength)

for axis in ['top','bottom','left','right']:
	ax1.spines[axis].set_linewidth(3)

#================================================================================================================
# PRC (Npas4--)
ax2 = fig.add_subplot(gs[0,1])

ax2.scatter(XKOR, YKOR, s=dotsize1, color='black', label='KO', facecolors='black')

ax2.axhline(0, color='k')
ax2.axis([xlR, xrR, ybR, ytR])
ax2.set_aspect('equal')

ax2.set_xlabel('Circadian Time (h)', fontsize=25, labelpad=5)
ax2.set_ylabel('Phase shift (h)', fontsize=25, labelpad=-12)

ax2.set_xticks(np.arange(0, 24.1, 6))
ax2.set_xticklabels(['0', '6', '12', '18', '24'], fontsize=20)
ax2.set_yticks(np.arange(-12, 12.1, 6))
ax2.set_yticklabels(['-12', '-6', '0', '6', '12'], fontsize=20)

ax2.tick_params(top=True, right=True)
ax2.tick_params(axis='x', labelsize=20, direction='in', length=ticklength)
ax2.tick_params(axis='y', labelsize=20, direction='in', length=ticklength)

for axis in ['top','bottom','left','right']:
	ax2.spines[axis].set_linewidth(3)

#================================================================================================================
# PTC (Npas4++)
oldphaseWT = XWTR
newphaseWT = XWTR + YWTR

dotcolor = 'black'
dotalpha = 1

ax3 = fig.add_subplot(gs[1,0])

ax3.scatter(oldphaseWT, newphaseWT, s=dotsize2, color=dotcolor, label='WT', alpha=dotalpha)
### Double plot
# Top left
ax3.scatter(oldphaseWT, newphaseWT + Tcycle, s=dotsize2, color=dotcolor, label='WT', alpha=dotalpha)
# Top right
ax3.scatter(oldphaseWT + Tcycle, newphaseWT + Tcycle, s=dotsize2, color=dotcolor, label='WT', alpha=dotalpha)
# Bottom right
ax3.scatter(oldphaseWT + Tcycle, newphaseWT, s=dotsize2, color=dotcolor, label='WT', alpha=dotalpha)

ax3.axis([xlT, xrT, ybT, ytT])
ax3.set_aspect('equal')

ax3.axhline(Tcycle, color='k')
ax3.axvline(Tcycle, color='k')

plw =3
pcolor='red'

l = Line2D([0, Tcycle*2],[0, Tcycle*2], color=pcolor, lw=plw, zorder=-1)
ax3.add_line(l)

l = Line2D([0, Tcycle],[Tcycle, Tcycle*2], color=pcolor, lw=plw, zorder=-1)
ax3.add_line(l)

l = Line2D([Tcycle, Tcycle*2],[0, Tcycle], color=pcolor, lw=plw, zorder=-1)
ax3.add_line(l)

ax3.set_xlabel('Old Phase (h)', fontsize=25, labelpad=5)
ax3.set_ylabel('New Phase (h)', fontsize=25, labelpad=5)

ax3.set_xticks(np.arange(0, 48.1, 12))
ax3.set_xticklabels(['0', '12', '24', '36', '48'], fontsize=20)
ax3.set_yticks(np.arange(0, 48.1, 12))
ax3.set_yticklabels(['0', '12', '24', '36', '48'], fontsize=20)

ax3.tick_params(axis='x', labelsize=20, direction='in', length=ticklength)
ax3.tick_params(axis='y', labelsize=20, direction='in', length=ticklength)

ax3.tick_params(top=True, right=True)
for axis in ['top','bottom','left','right']:
	ax3.spines[axis].set_linewidth(3)

#================================================================================================================
# PTC (Npas4--)
oldphaseKO = XKOR
newphaseKO = XKOR + YKOR

ax4 = fig.add_subplot(gs[1,1])

ax4.scatter(oldphaseKO, newphaseKO, s=dotsize2, color=dotcolor, label='KO', alpha=dotalpha)
### Double plot
# Top left
ax4.scatter(oldphaseKO, newphaseKO + Tcycle, s=dotsize2, color=dotcolor, label='KO', alpha=dotalpha)
# Top right
ax4.scatter(oldphaseKO + Tcycle, newphaseKO + Tcycle, s=dotsize2, color=dotcolor, label='KO', alpha=dotalpha)
# Bottom right
ax4.scatter(oldphaseKO + Tcycle, newphaseKO, s=dotsize2, color=dotcolor, label='KO', alpha=dotalpha)

ax4.axhline(Tcycle, color='k')
ax4.axvline(Tcycle, color='k')
ax4.axis([xlT, xrT, ybT, ytT])
ax4.set_aspect('equal')

l = Line2D([0, Tcycle*2],[0, Tcycle*2], color=pcolor, lw=plw, zorder=-1)
ax4.add_line(l)
l = Line2D([0, Tcycle],[Tcycle, Tcycle*2], color=pcolor, lw=plw, zorder=-1)
ax4.add_line(l)
l = Line2D([Tcycle, Tcycle*2],[0, Tcycle], color=pcolor, lw=plw, zorder=-1)
ax4.add_line(l)

ax4.set_xlabel('Old Phase (h)', fontsize=25, labelpad=5)
ax4.set_ylabel('New Phase (h)', fontsize=25, labelpad=5)

ax4.set_xticks(np.arange(0, 48.1, 12))
ax4.set_xticklabels(['0', '12', '24', '36', '48'], fontsize=20)
ax4.set_yticks(np.arange(0, 48.1, 12))
ax4.set_yticklabels(['0', '12', '24', '36', '48'], fontsize=20)

ax4.tick_params(axis='x', labelsize=20, direction='in', length=ticklength)
ax4.tick_params(axis='y', labelsize=20, direction='in', length=ticklength)
ax4.tick_params(top=True, right=True)

for axis in ['top','bottom','left','right']:
	ax4.spines[axis].set_linewidth(3)

ax1.text(-.18, 1,'F', ha='center', va='center', transform=ax1.transAxes, fontsize=30)
ax2.text(-.18, 1,'G', ha='center', va='center', transform=ax2.transAxes, fontsize=30)
ax3.text(-.18, 1,'H', ha='center', va='center', transform=ax3.transAxes, fontsize=30)
ax4.text(-.18, 1,'I', ha='center', va='center', transform=ax4.transAxes, fontsize=30)

gs.tight_layout(fig, rect=[0, 0, 1, 1])

plt.savefig('Pin_Figure6_PRC_PTC.pdf', dpi=300);

