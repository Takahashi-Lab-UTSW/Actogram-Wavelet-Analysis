## Maximal Overlap Discrete Wavelet Transform (MODWT)
The activity onsets were determined by maximal overlap discrete wavelet transform (MODWT) based on a Matlab (Mathworks) code created by Leise, T.L. [1, 2]. With a 4-tap Daubechies mother wavelet, MODWT decomposes the 15-min sampled actogram data into 7 details (D1 ~ D7) associated with particular period ranges (e.g., D3: 2 ~ 4 h) and an approximation representing all other remaining coarse scales. Then, an onset time on a particular day was defined to be one of several local peaks in a D3 detail time series around the actual onset.

	Reference:
	  [1] Leise, T.L. et al., Wavelet Meets Actogram, J. Biol. Rhythms 28, 62-68 (2013)
	  [2] Discrete wavelet transform to find activity onset
	  	https://tleise.people.amherst.edu/CircadianWaveletAnalysis.html
			WaveletActivityOnsetDetection.m

## Sensitivity
As in the Leise’s paper, the position of the onset detected was dependent on the percentage of the sensitivity which controlled the level of threshold for the local peak detection. Because the amplitude of the actogram was not maintained constant during a long period of recording time (~ 5 months per mouse), the different values of sensitivity were applied during the time interval before and after each light pulse in the same actogram. The optimal value of the sensitivity was determined as one that gave the smallest product of two root mean square error (RMSE) values from the two lines which fit onsets before and after each light pulse. In this procedure, the first three days after the light pulse was skipped to improve the fitting score due to an unstable nature of the activity right after the light pulse. Moreover, outliers which were significantly deviated from the fitting lines were automatically suppressed by assigning lower weights than others using iteratively reweighted least squares (“fitlm” function with robust linear regression option in Matlab).

![Alt text](README_figures/Actogram_onsets_with_different_sensitivities.png?raw=true "Actogram onsets with different sensitivities")
	Red/blue dot: onset before/after the light pulse
	Green line: fitting line
	Cyan rectangle: light

## Goodness of fit (adjusted R-squared and RMSE) vs. Sensitivity

![Alt text](README_figures/R2_RMSE_vs_Sensitivity.png?raw=true "R-squared and RMSE vs. Sensitivity")

	Green dot: sensitivity at which the product of two RMSE values before and after the light pulse is the lowest.


## Actogram with different sensitivity per light pulse (scan the lowest RMSE product)

![Alt text](README_figures/Actogram_Npas4--_C2-155-11.png?raw=true "R-squared and RMSE vs. Sensitivity")
![Alt text](https://latex.codecogs.com/svg.latex?\phi): CT, 
![Alt text](https://latex.codecogs.com/svg.latex?\Delta\phi): phaseshift

![Alt text](https://latex.codecogs.com/svg.latex?R^2_b),
![Alt text](https://latex.codecogs.com/svg.latex?R^2_a): R-squared before and after the light pulse

![Alt text](https://latex.codecogs.com/svg.latex?RMSE_b),
![Alt text](https://latex.codecogs.com/svg.latex?RMSE_a): RMSE before and after the light pulse

![Alt text](https://latex.codecogs.com/svg.latex?\tau_b),
![Alt text](https://latex.codecogs.com/svg.latex?\tau_a): period before and after the light pulse

S: sensitivity

