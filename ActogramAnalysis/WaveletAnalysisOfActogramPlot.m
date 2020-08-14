function WaveletAnalysisOfActogramPlot(Genotype, ActogramFilename, LightInfoFilename, ActogramTimeBin, OnsetFittingDay, OnsetFittingDayOffset, LightPulseDuration, NumLightPulses, MaxR2ProdIdxAll, CountPerMinAllBin, OnsetsAllCell, SensitivityAll, R2Cutoff, RMSECutoff)

% Import a file of light pulse information
LightInfo = CalLightInfo(LightInfoFilename);
DD = LightInfo(1);
for i=1:NumLightPulses
    LightPulseDay(i) = LightInfo(2*i);
    LightPulseTime(i) = LightInfo(2*i+1);
end

LastLDDay = DD-1;
ArbitraryMaxDay = 1000;
LightPulseDayAll = [0 LightPulseDay(LightPulseDay>0) ArbitraryMaxDay];

% Actogram data
x = CountPerMinAllBin;
NstepsPerHr=double(60./ActogramTimeBin);
Tcycle=24.;
Nsteps=length(x);

t=(0:Nsteps-1)'/NstepsPerHr;

% Plot actogram
ll=0.12; bb=0.06; wx=-0.45; wh=-0.12;
Actogram2xOuterNorthWesternStyle(x, NstepsPerHr, t(end)-t(1), Tcycle, ll, bb, wx, wh);

hold on;

ActogramFilenamePrefixTemp = strsplit(ActogramFilename, '_');
ActogramFilenamePrefix = cell2mat(ActogramFilenamePrefixTemp(1));
title(sprintf('%s   %s', Genotype, ActogramFilenamePrefix), 'FontSize', 20, 'HorizontalAlignment', 'center');

% light on in LD cycle
rectangle('Position', [6 0 12 LastLDDay], 'FaceColor', [0 1 1 0.5], 'EdgeColor', 'blue');
% light on in LD cycle (double plot)
rectangle('Position', [6+24 0 12 LastLDDay-1], 'FaceColor', [0 1 1 0.5], 'EdgeColor', 'blue');

% light pulses
for i=1:NumLightPulses
	rectangle('Position', [LightPulseTime(i) LightPulseDay(i) LightPulseDuration 1], 'FaceColor', [0 1 1 0.5], 'EdgeColor', 'blue');
end
% light pulses double plot
for i=1:NumLightPulses
	rectangle('Position', [LightPulseTime(i)+24 LightPulseDay(i)-1 LightPulseDuration 1], 'FaceColor', [0 1 1 0.5], 'EdgeColor', 'blue');
end


% Plot onsets
% 	onset color: alternates every other light pulse. For example,
% 	onsets before and after the 1st light pulse: red
% 	onsets before and after the 2nd light pulse: blue
% 	onsets before and after the 3nd light pulse: red
coloro = ['r', 'b'];
coloro = repmat(coloro,1,ceil(NumLightPulses/2));
% 	onset marker size
ms = 50;

for i=1:NumLightPulses
    if LightPulseDay(i) > 0
		onsetall = OnsetsAllCell{MaxR2ProdIdxAll(i)};
		onsetallpart = onsetall((onsetall(:,2) >= LightPulseDayAll(i+1) - OnsetFittingDay) & (onsetall(:,2) < LightPulseDayAll(i+1) + OnsetFittingDayOffset + OnsetFittingDay),:);
		
		XYb = onsetall((onsetall(:,2) >= LightPulseDayAll(i+1) - OnsetFittingDay) & (onsetall(:,2) < LightPulseDayAll(i+1)),:);
		XYa = onsetall((onsetall(:,2) >= LightPulseDayAll(i+1) + OnsetFittingDayOffset) & (onsetall(:,2) < LightPulseDayAll(i+1) + OnsetFittingDayOffset + OnsetFittingDay),:);

		scatter(onsetallpart(:,1), onsetallpart(:,2), ms, 'MarkerFaceColor', coloro(i), 'MarkerEdgeColor', coloro(i));
		%double plot
		scatter(onsetallpart(:,1)+24, onsetallpart(:,2)-1, ms, 'MarkerFaceColor', coloro(i), 'MarkerEdgeColor', coloro(i));

		Xbo = XYb(:,1);
		Ybo = XYb(:,2);
		Xao = XYa(:,1);
		Yao = XYa(:,2);

        [Xbo, Ybo, Xao, Yao, pbo, FitSSEbo, Rsquaredbo, Fstatbo, RMSEbo, pao, FitSSEao, Rsquaredao, Fstatao, RMSEao, PhaseAtLightPulseo, PhaseShifto] = FittingOnsetsWrapped(Xbo, Ybo, Xao, Yao, LightPulseDay(i), LightPulseTime(i), Tcycle, 1);

        PhaseInfoo = horzcat(PhaseAtLightPulseo, PhaseShifto, Rsquaredbo, Rsquaredao, pbo, pao, FitSSEbo, FitSSEao, Fstatbo, Fstatao, RMSEbo, RMSEao);

        if i==1
            PhaseInfoAllo = PhaseInfoo;
			xyo_new = [Xbo Ybo; Xao, Yao];
        else
            PhaseInfoAllo = vertcat(PhaseInfoAllo, PhaseInfoo);
			xyo_new = [xyo_new; [Xbo Ybo; Xao, Yao]];
        end

		taub = pbo(2) + 24;
		taua = pao(2) + 24;

		if (pbo(2)*pao(2)>0)
			if ((RMSEbo <=RMSECutoff) & (RMSEao <=RMSECutoff))
				if ((Rsquaredbo >=R2Cutoff) & (Rsquaredao >=R2Cutoff))
					text(50, LightPulseDayAll(i+1), sprintf('\\phi = %.1f (h), \\Delta\\phi = %.1f (h)\nR^2b = %.2f, {R^2}a = %.2f\n{RMSE_b} = %.2f, {RMSE_a} = %.2f\n{\\tau_b} = %.1f (h), {\\tau_a} = %.1f (h)\nS = %.1f', PhaseAtLightPulseo, PhaseShifto, Rsquaredbo, Rsquaredao, RMSEbo, RMSEao, taub, taua, SensitivityAll(MaxR2ProdIdxAll(i))), 'clipping', 'off', 'FontSize', 7, 'Color', 'red');
				elseif ((Rsquaredbo <R2Cutoff) | (Rsquaredao <R2Cutoff))
					text(50, LightPulseDayAll(i+1), sprintf('\\phi = %.1f (h), \\Delta\\phi = %.1f (h)\nR^2b = %.2f, {R^2}a = %.2f\n{RMSE_b} = %.2f, {RMSE_a} = %.2f\n{\\tau_b} = %.1f (h), {\\tau_a} = %.1f (h)\nS = %.1f', PhaseAtLightPulseo, PhaseShifto, Rsquaredbo, Rsquaredao, RMSEbo, RMSEao, taub, taua, SensitivityAll(MaxR2ProdIdxAll(i))), 'clipping', 'off', 'FontSize', 7, 'Color', 'blue');
				end

			elseif (((Rsquaredbo >=R2Cutoff) & (Rsquaredao >=R2Cutoff)) & ((RMSEbo >RMSECutoff) | (RMSEao >RMSECutoff)))
					text(50, LightPulseDayAll(i+1), sprintf('\\phi = %.1f (h), \\Delta\\phi = %.1f (h)\nR^2b = %.2f, {R^2}a = %.2f\n{RMSE_b} = %.2f, {RMSE_a} = %.2f\n{\\tau_b} = %.1f (h), {\\tau_a} = %.1f (h)\nS = %.1f', PhaseAtLightPulseo, PhaseShifto, Rsquaredbo, Rsquaredao, RMSEbo, RMSEao, taub, taua, SensitivityAll(MaxR2ProdIdxAll(i))), 'clipping', 'off', 'FontSize', 7, 'Color', 'green');

			else	
				text(50, LightPulseDayAll(i+1), sprintf('\\phi = %.1f (h), \\Delta\\phi = %.1f (h)\nR^2b = %.2f, {R^2}a = %.2f\n{RMSE_b} = %.2f, {RMSE_a} = %.2f\n{\\tau_b} = %.1f (h), {\\tau_a} = %.1f (h)\nS = %.1f', PhaseAtLightPulseo, PhaseShifto, Rsquaredbo, Rsquaredao, RMSEbo, RMSEao, taub, taua, SensitivityAll(MaxR2ProdIdxAll(i))), 'clipping', 'off', 'FontSize', 7);
			end
		else
			text(50, LightPulseDayAll(i+1), sprintf('\\phi = %.1f (h), \\Delta\\phi = %.1f (h)\nR^2b = %.2f, {R^2}a = %.2f\n{RMSE_b} = %.2f, {RMSE_a} = %.2f\n{\\tau_b} = %.1f (h), {\\tau_a} = %.1f (h)\nS = %.1f', PhaseAtLightPulseo, PhaseShifto, Rsquaredbo, Rsquaredao, RMSEbo, RMSEao, taub, taua, SensitivityAll(MaxR2ProdIdxAll(i))), 'clipping', 'off', 'FontSize', 7, 'Color', 'magenta');
		end

	end

end

%{
ActogramFilenamePrefixTemp = strsplit(ActogramFilename, '_');
ActogramFilenamePrefix = cell2mat(ActogramFilenamePrefixTemp(1));
dlmwrite(sprintf('%s_OnsetAll_BestSensitivity.csv', ActogramFilenamePrefix), xyo_new, 'delimiter', '\t');
%}

export_fig(sprintf('Actogram_%s_%s.pdf', Genotype, ActogramFilenamePrefix));
close all;
