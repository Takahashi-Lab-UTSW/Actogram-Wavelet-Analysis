<!--Phase response curve
	The phase shift in response to a 6-h light pulse was determined from the steady-state phase of activity onset for 14 days preceding and after the day of the light pulse.!-->

## Maximal Overlap Discrete Wavelet Transform (MODWT)
The activity onsets were determined by maximal overlap discrete wavelet transform (MODWT) based on a Matlab (Mathworks) code created by Leise, T.L. [1, 2]. With a 4-tap Daubechies mother wavelet, MODWT decomposes the 15-min sampled actogram data into 7 details (D1 ~ D7) associated with particular period ranges (e.g., D3: 2 ~ 4 h) and an approximation representing all other remaining coarse scales. Then, an onset time on a particular day was defined to be one of several local peaks in a D3 detail time series around the actual onset.

	Reference:
	  [1] Leise, T.L. et al., Wavelet Meets Actogram, J. Biol. Rhythms 28, 62-68 (2013)
	  [2] Discrete wavelet transform to find activity onset
	  	https://tleise.people.amherst.edu/CircadianWaveletAnalysis.html
			WaveletActivityOnsetDetection.m


## Sensitivity
As in the Leise’s paper, the position of the onset detected was dependent on the percentage of the sensitivity which controlled the level of threshold for the local peak detection. Because the amplitude of the actogram was not maintained constant during a long period of recording time (~ 5 months per mouse), the different values of sensitivity were applied during the time interval before and after each light pulse in the same actogram. The optimal value of the sensitivity was determined as one that gave the smallest product of two root mean square error (RMSE) values from the two lines which fit onsets before and after each light pulse. In this procedure, the first three days after the light pulse was skipped to improve the fitting score due to an unstable nature of the activity right after the light pulse. Moreover, outliers which were significantly deviated from the fitting lines were automatically suppressed by assigning lower weights than others using iteratively reweighted least squares (“fitlm” function with robust linear regression option in Matlab).

![Alt text](ActogramAnalysis/Actogram_onsets_with_different_sensitivities.png?raw=true "Actogram onsets with different sensitivities")
<!--![Alt text](pipeline.png?raw=true "piecewise rigid motion correction pipeline")!-->
	Red/blue dot: onset before/after the light pulse
	Green line: fitting line
	Cyan rectangle: light

## Goodness of fit (adjusted R-squared and RMSE) vs. Sensitivity

![Alt text](ActogramAnalysis/R2_RMSE_vs_Sensitivity.png?raw=true "R-squared and RMSE vs. Sensitivity")

	Green dot: sensitivity at which the product of two RMSE values before and after the light pulse is the lowest.


## Actogram with different sensitivity per light pulse (scan the lowest RMSE product)

![Alt text](ActogramAnalysis/Actogram_Npas4--_C2-155-11.png?raw=true "R-squared and RMSE vs. Sensitivity")


<!--Then, the phase of the light pulse was calculated from the circadian period (slope of the fitting line) and the time difference between CT 12 (phase of activity onset) and the beginning of the light pulse. The magnitude of the phase shift was the time difference between the intersection points of the two lines with the day of the light pulse. The phase response curve (PRC, Figure 6F and 6G) to the light pulse was plotted with the data which met the two conditions: 1) two adjusted R-squared values of the line fitting onsets before and after the light pulses were both larger than 0.5, 2) two RMSE values of the line fitting onsets before and after the light pulses were both smaller than 3.0 h. !-->

<!--
<img src="https://raw.githubusercontent.com/vokkmi/Takahashi-Lab-UTSW/Wavelet-Analysis/.github/images/Sensitivity_Fig1.png" align="center">

	black thick line: actogram count per 15 min (15 min binning, raw data: count/min)
	blue thick line: discrete wavelet-transformed data with a period range 2 ~ 4 h (D3, dwcoeff (3h))
	black vertical lines: peak detection interval (variable depending on the number of peaks of dwcoeff (3h) signal for 24 hrs)
	black horizontal dotted lines: mean of actogram count between 3 days before and after the interval
	rectangle (magenta and cyan): sum of count for 0 ~ 1 hr and 1 hr ~ 3 hrs after the peak time.
	dots (red, green, gray): peaks within the interval of dwcoeff (3h) signal

	onset requirement: 1) peaks > threshold, 2) sum of count > mean count, 3) sum of count (0 ~ 1 hr) < sum of count (1 hr ~ 3 hrs) *0.46
		* threshold: maximum peak value*sensitivity (0 ~ 1)
		requirement 1): red and green dots (gray dot: under the threshold)
		requirements 2) & 3): all dots
	→ red dot (onset): the earliest peak of the peaks fulfilling all the requirements 


In this procedure, the first three days after the light pulse was skipped to improve the fitting score due to an unstable nature of the activity right after the light pulse. Moreover, outliers which are significantly deviated from the fitting lines were automatically removed with a robust linear regression model (“fitlm” in Matlab).
	Then, the phase of the light pulse was calculated from the circadian period (slope of the fitting line) and the time difference between CT12 (phase of activity onset) and the beginning of the light pulse. The magnitude of the phase shift was the time difference between the intersection points of the two lines with the day of the light pulse. The phase response curve (PRC, Figure 6F and 6G) to the light pulse was plotted with the data which meet the two conditions: 1) two R-squared values of the line fitting before and after the light pulses are both larger than 0.5, 2) two RMSE values of the line fitting before and after the light pulses are both smaller than 3.0 h.
	The numbers of mice used for the PRC are 25 for Npas4+/+ and 21 for Npas4-/- mice. Individual mice were given 3 ∼ 7 light pulses every ∼ 21 days in the middle of bedding change. The total number (N) of light pulses plotted were 206 (126 for Npas4+/+ and 80 for Npas4-/- mice). The raw actogram data from all mice and the source code for the time series analysis with MODWT are available in https://github.com/Takahashi-Lab-UTSW/Wavelet-Analysis. Circadian period in constant darkness (DD) or constant light (LL) in Figure S6B was estimated using the Chi-square periodogram in ClockLab for 10-day intervals.

    Statistics
	For statistical comparisons between Npas4+/+ and Npas4-/- PRC, the data on the graph were grouped every CT 3 hours (CT 0, CT 3, CT 6, ..., CT 21; e.g. a group CT0 is the data between CT 22.5 an CT 1.5) separately in Npas4+/+ and Npas4-/- PRC. Because the number of data in each CT group was different, an unbalanced two-way factorial ANOVA (using SAS University Edition) was applied to the data. The analysis with type III sums of squares showed that there was a significant interaction between CT and genotype (p-value < 0.0001). Further analysis showed that the effect of CT on the phaseshift depended on the genotype significantly at CT 0, 9, 12, 15, and 21 (p-values < 0.05, < 0.001, < 0.01, < 0.0001, and < 0.001 respectively). A SAS code for the analysis and the result are available in the above github repository as well.



Wavelet-Analysis


	RunAllActogramData.m
		IterativeWaveletAnalysisOfActogram.m
			WaveletAnalysisOfActogram.m
				FittingOnsetsWrapped.m
			WaveletAnalysisOfActogramPlot.m
				FittingOnsetsWrapped.m
				Actogram2xOuterNorthWesternStyle.m


IterativeWaveletAnalysisOfActogram.m
	Iterative wavelet analysis process to find onsets (peaks) in circidian actogram data
	with different sensitivity (~ threshold) values which are separately applied to
	time intervals (before and after) of individual light pulses

	Usage: IterativeWaveletAnalysisOfActogram('Npas4++', 'A4-96-10_Counts.csv', 'LightPulseInfo.csv', 15, 14, 3, 6, 7, 0.1)
		  Npas4++
			  Genotype
		  A4-96-10_Counts.csv
			  ClockLab actogram data
				  Day, Hr, Min, Cnts/min, Lights
		  LightPulseInfo.csv
			  Light pulse information
				  LightCycleDescription Year Month Day Hour Minute Second
				  Start: initial time to start recording actogram in ClockLab
				  DD: constant darkness
				  LP1: time at the 1st light pulse
		  15: actogram binning time (minute)
		  14: time interval (days) before or after the light pulse for fitting onsets 
			  ** assume that the bedding change doesn't affect the phase shift very much.
		  3:  time offset (days) after the light pulse 
		  6: light pulse duration (hours)
		  7: number of light pulses
		  0.1: binning width of the sensitivity values (0 ~ 1)
			  0, 0.1, 0.2, ..., 1.0






# NoRMCorre: Non-Rigid Motion Correction 
This package provides a Matlab implementation of the NoRMCorre algorithm [[1]](#ref), and can be used for online piecewise rigid motion correction of 2d (planar) or 3d (volumetric) calcium imaging data. 

## Citation

If you find this package useful please cite the companion paper [[1]](#ref):

```
@article{pnevmatikakis2017normcorre,
  title={NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data},
  author={Pnevmatikakis, Eftychios A and Giovannucci, Andrea},
  journal={Journal of neuroscience methods},
  volume={291},
  pages={83--94},
  year={2017},
  publisher={Elsevier}
}
```

## Synopsis

The algorithm operates by splitting the field of view into a set of overlapping patches. For each patch and each frame a rigid translation is estimated by aligning the patch against a template using an efficient, FFT based, algorithm for subpixel registration [[2]](#reg). The estimated set of translations is further upsampled to a finer resolution to create a smooth motion field that is applied to a set of smaller overlapping patches. Extra care is taken to avoid smearing caused by interpolating overlapping patches with drastically different motion vectors. The registered frame is used to update the template in an online fashion by calculating a running/mean of past registered frames. The pipeline is summarized in the figure below.

![Alt text](pipeline.png?raw=true "piecewise rigid motion correction pipeline")

## Code details

See the function [```demo.m```](https://github.com/simonsfoundation/NoRMCorre/blob/master/demo.m) for an example of the code. The algorithm is implemented in the function ```normcorre.m```. If you have access to the parallel computing toolbox, then the function ```normcorre_batch.m``` can offer speed gains by enabling within mini-batch parallel processing. The user gives a dataset (either as 3D or 4D tensor loaded in RAM or memory mapped, or a pointer to a .tiff stack or .hdf5 file), and a parameters struct ```options```. Optionally, an initial template can also be given. The algorithm can also be used for motion correction of 1p micro-endoscopic data, by estimating the shifts on high pass spatially filtered version of the data. See the script [```demo_1p.m```](https://github.com/simonsfoundation/NoRMCorre/blob/master/demo_1p.m) for an example.

The algorithm can also be ran using the MotionCorrection object. See [```demo_mc_class.m```](https://github.com/simonsfoundation/NoRMCorre/blob/master/demo_mc_class.m) for an example on how to use the object for 2p and 1p data.

The ```options``` struct can be set either manually or by using the function ```NoRMCorreSetParms.m```. Some parameters of the ```options``` struct are the following:

| Parameter name | Description |
|----------------|-------------|
| ```d1,d2,d3``` | dimensions of field of view |
| ```grid_size``` | size of non-overlapping portion of each patch the grid in each direction (x-y-z)|
| ```overlap_pre```| size of overlapping region in each direction before upsampling  |
| ```mot_uf```    | upsampling factor for smoothing and refinement of motion field |
| ```overlap_post ``` | size of overlapping region in each direction after upsampling |
| ```max_shift``` | maximum allowed shift for rigid translation | 
| ```max_dev``` | maximum deviation of each patch from estimated rigid translation |
| ```upd_template``` | update the template online after registering some frames |
| ```bin_width``` | length of bin over which the registered frames are averaged to update the template |
| ```init_batch``` | number of frames to be taken for computing initial template |
| ```iter``` | number of times to go over the dataset |
| ```output_type``` | type of output registered file |
| ```phase_flag``` | flag for using phase correlation |
| ```correct_bidir``` | check for offset due to bidirectional scanning (default: true) |

The performance of registration can be evaluated using the function ```motion_metrics.m```. The function simply computes the correlation coefficient of each (registered) frame, with the mean (registered) frame across time, the mean registered frame, and its crispness.

## Developers

[Eftychios A. Pnevmatikakis](https://github.com/epnev), Flatiron Institure, Simons Foundation

## External packages

This package includes functions from the following packages
- [Save and load a multipage tiff file](https://www.mathworks.com/matlabcentral/fileexchange/35684-save-and-load-a-multiframe-tiff-image/content/loadtiff.m)
- [Savefast](https://www.mathworks.com/matlabcentral/fileexchange/39721-save-mat-files-more-quickly) for saving (and then loading) MAT files more quickly without compressing their contents. 
- [Eficient subpixel registration by cross-correlation](https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation) for fast alignment of an image against a template.

## Integrations 
This package will be integrated with the [Matlab code](https://www.github.com/epnev/ca_source_extraction) for source extraction and deconvolution using CNMF.

## References 
% Discrete wavelet transform to find activity onset
% https://tleise.people.amherst.edu/CircadianWaveletAnalysis.html
% WaveletActivityOnsetDetection.m
%A website for   [here](https://tleise.people.amherst.edu/CircadianWaveletAnalysis.html).

<a name="ref"></a>[1] Tanya L Leise, Premananda Indic, Matthew J Paul, William J Schwartz, Wavelet Meets Actogram, *J. Biol. Rhythms* 28, 62-68 (2013); doi: [https://doi.org/10.1177/0748730412468693].


