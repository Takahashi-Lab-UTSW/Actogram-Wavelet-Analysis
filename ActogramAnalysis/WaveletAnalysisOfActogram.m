function [PhaseInfoAllo, CountPerMinAllBin, Onsets] = WaveletAnalysisOfActogram(ActogramFilename, LightInfoFilename, ActogramTimeBin, OnsetFittingDay, OnsetFittingDayOffset, NumLightPulses, Sensitivity) 

% Import an actogram file
ActogramData = importdata(ActogramFilename);
AD = ActogramData.data;
AD(isnan(AD))=0;
DayAll = AD(:,1);
HourAll = AD(:,2);
MinAll = AD(:,3);
CountPerMinAll = AD(:,4);
LightAll = AD(:,5);

Tcycle=24.;
edges = 0:ActogramTimeBin:(length(AD) + ActogramTimeBin);
StartTime = DayAll(1)*Tcycle*60 + HourAll(1)*60 + MinAll(1);
DateMin = DayAll*Tcycle*60 + HourAll*60 + MinAll - StartTime;

% Binning the actogramd data
[~, idx] = histc(DateMin,edges);
BinSums = accumarray(idx,CountPerMinAll);
ActogramBin = [DayAll(1:ActogramTimeBin:end) HourAll(1:ActogramTimeBin:end) MinAll(1:ActogramTimeBin:end) BinSums LightAll(1:ActogramTimeBin:end) edges(1:end-1)'];

DayAllBin = ActogramBin(:,1);
HourAllBin = ActogramBin(:,2);
MinAllBin = ActogramBin(:,3);
CountPerMinAllBin = ActogramBin(:,4);
LightAllBin = ActogramBin(:,5);

% Import a file of light pulse information
LightInfo = CalLightInfo(LightInfoFilename);
DD = LightInfo(1);
for i=1:NumLightPulses
	LightPulseDay(i) = LightInfo(2*i);
	LightPulseTime(i) = LightInfo(2*i+1);
end

%================================================================================================================================
% Discrete wavelet transform to find activity onset
% Leise, T.L. et al., Wavelet Meets Actogram, J. Biol. Rhythms 28, 62-68 (2013)
% https://tleise.people.amherst.edu/CircadianWaveletAnalysis.html
% WaveletActivityOnsetDetection.m

x= CountPerMinAllBin;
NstepsPerHr=double(60./ActogramTimeBin);
Nsteps=length(CountPerMinAllBin);

t=(0:Nsteps-1)'/NstepsPerHr;

filt='d4'; % 4-tap Daubechies filter works well
J0=ceil(log2(Tcycle*NstepsPerHr)); % J0 will be one level beyond circadian level; will be 7 for 15 min bins

[WJt, VJt, att] = modwt(x,filt,J0,'circular'); % calculate modwt coefficients
% WJt, DJt: N x 7 columns (N: length of x, 7: total number of details (J0), D1 - D7)
[DJt, ~]=imodwt_details(WJt, filt);

% Assumes 15-minute bins, so circadian level is 6
[mm,ind]=findpeaks(DJt(:,6));
ind=ind(2:end-1);
mm=mm(2:end-1);

peaktimes=t(ind);  % peaks of circadian component
peak2peak=peaktimes(2:end)-peaktimes(1:end-1);

% Calculates wavelet coefficients for 24h and 3h scales
[TWJt, TVJt] = modwt_cir_shift(WJt, VJt, filt, J0); % adjusted wavelet coefficients (undo shift caused by filtering)
%   TWJt          -  NxJ array of circularly shifted (advanced)
%	 				MODWT wavelet coefficents
% J0-1=6, J0-4=3

wcoeff_24h=TWJt(:,J0-1); % wavelet coefficients for levels of interest
wcoeff_3h=TWJt(:,J0-4);
dwcoeff_3h=wcoeff_3h(6:end)-wcoeff_3h(1:end-5);
dwcoeff_3h=[dwcoeff_3h;zeros(5,1)];

[pks24,ipks24]=findpeaks(wcoeff_24h);  % find activity peaks at 24h scale to delineate circadian days
ipks24=ipks24(2:end-1);
pks24=pks24(2:end-1); % remove edge effected values
pks24=pks24(t(ipks24)>Tcycle);
ipks24=ipks24(t(ipks24)>Tcycle);

% Removes spurious peaks with extremely low amplitude
ind=find(pks24<mean(pks24)/10);     
if ~isempty(ind)
    pks24(ind)=[];ipks24(ind)=[];
end

% Removes spurious peaks that are too closely spaced
spacing=(ipks24(2:end)-ipks24(1:end-1))/NstepsPerHr;

while min(spacing)<Tcycle/2       
    [~,ind]=min(spacing);
    if pks24(ind)<pks24(ind+1)
        pks24(ind)=[];
        ipks24(ind)=[];
    else
        pks24(ind+1)=[];
        ipks24(ind+1)=[];
    end
    spacing=(ipks24(2:end)-ipks24(1:end-1))/NstepsPerHr;
end

Npeaks24h=numel(ipks24);

% Finds troughs at 24h scale
trs24=zeros(Npeaks24h-1,1)+NaN;
itrs24=zeros(Npeaks24h-1,1)+NaN; 

[mm,ind]=min(wcoeff_24h(1:ipks24(1)));

itrs24(1)=ind;
trs24(1)=mm;

for day=2:Npeaks24h
    [mm,ind]=min(wcoeff_24h(ipks24(day-1):ipks24(day)));
    itrs24(day)=ind+ipks24(day-1)-1;trs24(day)=mm;
end

% Pre-allocates arrays
ipks3=zeros(Npeaks24h-1,1)+NaN;
itrs3=zeros(Npeaks24h-1,2)+NaN;

pks3=zeros(Npeaks24h-1,1)+NaN;
trs3=zeros(Npeaks24h-1,2)+NaN;

% Uses 24h scale wavelet coefficients combined with 3h scale to calculate onsets
for day=1:Npeaks24h-1

    quarterday=min(NstepsPerHr*9,ceil(0.55*NstepsPerHr*(t(ipks24(day))-t(itrs24(day)))));
    range3=max(0,ipks24(day)-quarterday)+1:ipks24(day);

    [mm,ind]=findpeaks(dwcoeff_3h(range3),'sortstr','descend');

    sumact=zeros(numel(ind),3);
    for j=1:numel(ind)
        sumact(j,:)=[sum(x(ind(j)+range3(1)-1+(NstepsPerHr:3*NstepsPerHr))) sum(x(ind(j)+range3(1)-1+(0:NstepsPerHr))) sum(x(ind(j)+range3(1)-1+(1*NstepsPerHr:12*NstepsPerHr)))];
    end

    if ~isempty(ind)
        indsig=find(mm>Sensitivity*mm(1) & sumact(:,1)>mean(x) & sumact(:,2)<0.46*sumact(:,1));

        if isempty(indsig)
            range3=max(0,ipks24(day)-round(1.6*quarterday))+1:ipks24(day);

            [mm,ind]=findpeaks(dwcoeff_3h(range3),'sortstr','descend');

            sumact=zeros(numel(ind),2);
            for j=1:numel(ind)
                sumact(j,:)=[sum(x(ind(j)+range3(1)-1+(NstepsPerHr:3*NstepsPerHr))) sum(x(ind(j)+range3(1)-1+(0:NstepsPerHr)))];
            end
            indsig=find(mm>.37*mm(1) & sumact(:,1)>mean(x) & sumact(:,2)<0.4*sumact(:,1));
        end

        if ~isempty(indsig)
            
            [earliest,iearliest]=min(ind(indsig));

            if numel(indsig)>1
                [indsorted,ix]=sort(ind(indsig));
			
				if size(sumact,2)>2
					if indsorted(2)-earliest>5*NstepsPerHr && sumact(iearliest,3)<max(sumact(indsig,3))
						earliest=indsorted(2);
						iearliest=ix(2);
					end
				end
            end
            itr=earliest;
			tr=mm(indsig(iearliest));

            itrs3(day,1)=itr+range3(1)-1;
            trs3(day,1)=tr;

            [tr,itr]=min(wcoeff_3h(ipks24(day)+1:ipks24(day)+quarterday));

            itrs3(day,2)=itr+ipks24(day);
            trs3(day,2)=tr;
        end
    end
    
    [pk,ipk]=findpeaks(wcoeff_3h(max(0,ipks24(day)-quarterday)+1:ipks24(day)+quarterday),'sortstr','descend');
    if ~isempty(ipk)
        ipks3(day)=ipk(1)+max(0,ipks24(day)-quarterday);
        pks3(day)=pk(1);

    end
end

trs3inds=find(~isnan(itrs3(:,1)));

% (x, y) = (time, day) position of onset
% To account for the phase shift caused by the wavelet filter,
% 1.25 h was added to the phase of the D 3 peak to yield the phase of activity onset.
% This offset was determined through testing on both simulated and real activity data.
% Leise, T.L. et al., Wavelet Meets Actogram, J. Biol. Rhythms 28, 62-68 (2013)

xo = mod(t(itrs3(trs3inds,1))-.1,Tcycle) + 0.1 + 1.25;
yo = ceil((t(itrs3(trs3inds,1)))/Tcycle);
%yo = ceil((t(itrs3(trs3inds,1)))/Tcycle) - 0.3;
Onsets = [xo, yo];

%==============================================================================================================================
% Fitting onsets
for i=1:NumLightPulses

	if LightPulseDay(i) > 0
		
		% Onsets before the i-th light pulse
		Xbo = xo((yo >= LightPulseDay(i) - OnsetFittingDay) & (yo < LightPulseDay(i)));
		Ybo = yo((yo >= LightPulseDay(i) - OnsetFittingDay) & (yo < LightPulseDay(i)));

		% Onsets after the i-th light pulse
		Xao = xo((yo >= LightPulseDay(i) + OnsetFittingDayOffset) & (yo < LightPulseDay(i) + OnsetFittingDayOffset + OnsetFittingDay));
		Yao = yo((yo >= LightPulseDay(i) + OnsetFittingDayOffset) & (yo < LightPulseDay(i) + OnsetFittingDayOffset + OnsetFittingDay));
		
		[Xbo, Ybo, Xao, Yao, pbo, FitSSEbo, Rsquaredbo, Fstatbo, RMSEbo, pao, FitSSEao, Rsquaredao, Fstatao, RMSEao, PhaseAtLightPulseo, PhaseShifto] = FittingOnsetsWrapped(Xbo, Ybo, Xao, Yao, LightPulseDay(i), LightPulseTime(i), Tcycle, 0);

		PhaseInfoo = horzcat(PhaseAtLightPulseo, PhaseShifto, Rsquaredbo, Rsquaredao, pbo, pao, FitSSEbo, FitSSEao, Fstatbo, Fstatao, RMSEbo, RMSEao);

		if i==1
			PhaseInfoAllo = PhaseInfoo;
		else
			PhaseInfoAllo = vertcat(PhaseInfoAllo, PhaseInfoo);
		end
	end

end

%{
% Save onsets 
ActogramFilenamePrefixTemp = strsplit(ActogramFilename, '_');
ActogramFilenamePrefix = cell2mat(ActogramFilenamePrefixTemp(1));
dlmwrite(sprintf('%s_OnsetAll_Sensitivity_%0.1f.csv', ActogramFilenamePrefix, Sensitivity), [xo yo], 'delimiter', '\t');
%}

close all;
