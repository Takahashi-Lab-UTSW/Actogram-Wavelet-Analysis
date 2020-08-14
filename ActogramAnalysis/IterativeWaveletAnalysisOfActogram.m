% Iterative wavelet analysis process to find onsets (peaks) in circidian actogram data
% with different sensitivity (~ threshold) values which are separately applied to
% time intervals (before and after) of individual light pulses

% Usage: IterativeWaveletAnalysisOfActogram('Npas4++', 'A4-96-10_Counts.csv', 'LightPulseInfo.csv', 15, 14, 3, 6, 7, 0.1)
%		Npas4++
%			Genotype
%		A4-96-10_Counts.csv
%			ClockLab actogram data
%				Day, Hr, Min, Cnts/min, Lights
%		LightPulseInfo.csv
%			Light pulse information
%				LightCycleDescription Year Month Day Hour Minute Second
%				Start: initial time to start recording actogram in ClockLab
%				DD: constant darkness
%				LP1: time at the 1st light pulse
%		15: actogram binning time (minute)
%		14: time interval (days) before or after the light pulse for fitting onsets 
%			** assume that the bedding change doesn't affect the phase shift very much.
%		3:	time offset (days) after the light pulse 
%		6: light pulse duration (hours)
%		7: number of light pulses
%		0.1: binning width of the sensitivity values (0 ~ 1)
%			0, 0.1, 0.2, ..., 1.0

% Reference:
%	Leise, T.L. et al., Wavelet Meets Actogram, J. Biol. Rhythms 28, 62-68 (2013)
%	Discrete wavelet transform to find activity onset
%		https://tleise.people.amherst.edu/CircadianWaveletAnalysis.html
% 		WaveletActivityOnsetDetection.m
	
% byeongha.jeong@utsouthwestern.edu, 08/11/20

function [MinRMSEProdSensitivityAll, PhaseInfoLPMinRMSEProdAll] = IterativeWaveletAnalysisOfActogram(Genotype, ActogramFilename, LightInfoFilename, ActogramTimeBin, OnsetFittingDay, OnsetFittingDayOffset, LightPulseDuration, NumLightPulses, SensitivityBinWidth, R2Cutoff, RMSECutoff)

warning('off','all');

NumSensitivityBin = 1/SensitivityBinWidth;

for i=0:NumSensitivityBin
	Sensitivity = i*SensitivityBinWidth;
	
	[PhaseInfo, CountPerMinAllBin, Onsets] = WaveletAnalysisOfActogram(ActogramFilename, LightInfoFilename, ActogramTimeBin, OnsetFittingDay, OnsetFittingDayOffset, NumLightPulses, Sensitivity);

	if i==0
		PhaseInfoAll = PhaseInfo;	
		OnsetsAllCell = {Onsets};
	else
		PhaseInfoAll = [PhaseInfoAll; PhaseInfo];
		OnsetsAllCell = [OnsetsAllCell, {Onsets}];
	end

end

SensitivityAll = [0:SensitivityBinWidth:1.0];

for i=1:NumLightPulses
	PhaseInfoLP = PhaseInfoAll([i:NumLightPulses:NumLightPulses*length(SensitivityAll)],:);

	RMSEProd = prod(PhaseInfoLP(1:end,end-1:end), 2)
    [MinRMSEProdVal, MinRMSEProdIdx]=min(RMSEProd);
    PhaseInfoLPMinRMSEProd = PhaseInfoLP(MinRMSEProdIdx,:);
	sprintf('The %d-th light pulse: the best sensitivity = %0.1f', i, SensitivityAll(MinRMSEProdIdx))

	if i==1
		MinRMSEProdIdxAll = MinRMSEProdIdx;
		PhaseInfoLPMinRMSEProdAll = PhaseInfoLPMinRMSEProd;
	else
		MinRMSEProdIdxAll = [MinRMSEProdIdxAll; MinRMSEProdIdx];
		PhaseInfoLPMinRMSEProdAll = [PhaseInfoLPMinRMSEProdAll; PhaseInfoLPMinRMSEProd];
	end	
end

MinRMSEProdSensitivityAll = SensitivityAll(MinRMSEProdIdxAll);

WaveletAnalysisOfActogramPlot(Genotype, ActogramFilename, LightInfoFilename, ActogramTimeBin, OnsetFittingDay, OnsetFittingDayOffset, LightPulseDuration, NumLightPulses, MinRMSEProdIdxAll, CountPerMinAllBin, OnsetsAllCell, SensitivityAll, R2Cutoff, RMSECutoff);
