#!/usr/bin/env ipython
# usage:
#	ConvertPRCInfotoTableforARTAnalysisImputed.py\
#		PRC_table_Binning_Npas4++_Imputed.csv\
#		PRC_table_Binning_Npas4--_Imputed.csv\
#		3
#	3: binwidth, grouping CT levels every 3 hours

import numpy as np
import sys

WT = np.loadtxt(sys.argv[1])
KO = np.loadtxt(sys.argv[2])
BinWidth = int(sys.argv[3])

NumRowsWT = WT.shape[0]
NumRowsKO = KO.shape[0]
NumCols = WT.shape[1]

CT = list(range(0, NumCols*BinWidth, BinWidth))
Labels = ['Npas4++', 'Npas4--']

NumRowsDiff = NumRowsWT - NumRowsKO 

NumRowsWT = WT.shape[0]
NumRowsKO = KO.shape[0]
NumRows = NumRowsWT
NumRowsWTKO = [NumRowsWT, NumRowsKO]

WTKO = (WT, KO)
for k in range(2):
	for j in range(NumCols):
		for i in range(NumRowsWTKO[k]):
			line = '{0}, CT{1}, {2}'.format(WTKO[k][i,j], CT[j], Labels[k])
			if i==0:
				lineAllRows = line
			else:
				lineAllRows = np.vstack([lineAllRows, line])
		
		if j==0:
			lineAllRowsCols = lineAllRows
		else:
			lineAllRowsCols = np.vstack([lineAllRowsCols, lineAllRows])
		
	if k==0:
		lineAllRowsColsGenotype = lineAllRowsCols
	else:
		lineAllRowsColsGenotype = np.vstack([lineAllRowsColsGenotype, lineAllRowsCols])

PRCTableforSAS = np.vstack(['Phaseshift,CT,Genotype', lineAllRowsColsGenotype])
np.savetxt("PRC_Table_for_ART-analysis_Imputed.csv", PRCTableforSAS, fmt='%s')
