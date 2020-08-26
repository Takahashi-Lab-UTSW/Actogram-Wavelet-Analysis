## Maximal Overlap Discrete Wavelet Transform (MODWT)
The activity onsets were determined by maximal overlap discrete wavelet transform (MODWT) based on a Matlab (Mathworks) code created by Leise, T.L. [1, 2]. With a 4-tap Daubechies mother wavelet, MODWT decomposes the 15-min sampled actogram data into 7 details (D1 ~ D7) associated with particular period ranges (e.g., D3: 2 ~ 4 h) and an approximation representing all other remaining coarse scales. Then, an onset time on a particular day was defined to be one of several local peaks in a D3 detail time series around the actual onset.

	References:
		[1] Percival, D.B. Wavelet Methods for Time Series Analysis,
				Cambridge University Press (2000) pp. 169 ~ 179
		[2] Leise, T.L. et al., Wavelet Meets Actogram, J. Biol. Rhythms 28, 62-68 (2013)

![Alt text](README_figures/C2-155-11_MODWT.png?raw=true "Actogram MODWT decomposition")

	x(t): actogram count per 15 min	
	TW_i: vector (time series) of circulary shifted MODWT wavelet (W) coefficients
	TV: vector of circulary shifted MODWT approximation (V) coefficients
	X(t): actogram reconstruction from wavelet and approximation coefficients by inverse MODWT

![Alt text](https://latex.codecogs.com/svg.latex?\Large&space;W_{j,t}&space;=&space;\sum_{l=0}^{L_j-1}&space;h_{j,l}&space;~x_{t-l&space;~mod~&space;N}~~~~~~~~~~V_{j,t}&space;=&space;\sum_{l=0}^{L_j-1}&space;g_{j,l}&space;~x_{t-l&space;~mod~&space;N})

	h_j,l : periodized j-th level MODWT wavelet filter ('d4', 4-tap Daubechies filter)
	g_j,l : periodized j-th level MODWT scale filter (quadrature mirror filter that corresponds to h_j,l)
	L_j: filter width
	N: total number of frames (length)

![Alt text](https://latex.codecogs.com/svg.latex?\Large&space;X(t)&space;=&space;S_{J0,t}&space;&plus;&space;\sum_{j=1}^{J0=7}&space;D_{j,t})

![Alt text](https://latex.codecogs.com/svg.latex?\Large&space;S_{J0,t}&space;=&space;\sum_{l=0}^{N-1}&space;g_{j,l}&space;~V_{j,~t&plus;l~&space;mod~&space;N}~~~~~~~~~~D_{j,t}&space;=&space;\sum_{l=0}^{N-1}&space;h_{j,l}&space;~W_{j,~t&plus;l~&space;mod~&space;N})

	S_J0: MODWT smooth
	D_j,t: j-th level MODWT detail


![Alt text](README_figures/C2-155-11_Day45_Sensitivity0.0_0.8.png?raw=true "Actogram MODWT decomposition")

	black thick line: actogram count per 15 min (15 min binning, raw data: count/min)
	blue thick line (dwcoeff (3 h)): difference between a time series of a circulary shifted wavelet coefficient
		(TW3 or wcoeff (3 h), period range: 2 ~ 4h) and it's 5-frame (5*15 min) adjacent vector
		dwcoeff (3 h) = TW_3 (6:end) - TW_3 (1:end-5)

	black dotted rectangle: peak detection interval
		(variable depending on the number of peaks of dwcoeff (3 h) signal for 24 hrs)
	light green horizontal line: global mean of actogram count (mean count)
	rectangle (magenta and cyan): sum of count for 0 ~ 1 hr and 1 hr ~ 3 hrs after the peak time.

	dots (red, green, gray): peaks within the interval of dwcoeff (3 h) signal

	The onset (red dot) is the earliest one of all peaks within the black dotted rectangle, which meets:
		1) peaks > threshold (= maximum peak (the highest green) value *sensitivity (0 ~ 1))
			→  red and green dots (gray dot: below the threshold)
		2) sum of count (0 ~ 1 hr from a peak, magenta rect) > mean count (light green line)
		3) sum of count (0 ~ 1 hr from a peak, magenta rect)
			< sum of count (1 hr ~ 3 hrs from a peak, cyan rect) *0.46

		Reference: Discrete wavelet transform to find activity onset
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

