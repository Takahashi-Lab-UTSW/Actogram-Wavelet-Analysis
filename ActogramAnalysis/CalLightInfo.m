% To convert the date and time of light pulses to time (days and hours) elapsed from the beginning
function LightInfoOut = CalLightInfo(fileLightPulseInfo)
% LightPulseInfo.csv
%	Light pulse information
%   	LightCycleDescription Year Month Day Hour Minute Second
%    	Start: initial time to start recording actogram in ClockLab
%    	DD: constant darkness
%    	LP1: time at the 1st light pulse

% Output: LightInfoOut
%	[DD-start (day), LP1-start (day), LP1-start (hour), LP2-start (day), LP2-start (hour), ...]
%	size: 1 x (1 + number of light pulses * 2)

fid = fopen(fileLightPulseInfo);
LightPulseInfo = textscan(fid,'%s %d %d %d %d %d %d');
fclose(fid);

LPstr = string(LightPulseInfo{1});
LPstart = find(LPstr =='LP1');

for i=2:size(LightPulseInfo,2)
	LPvec = cell2mat(LightPulseInfo(i));
	if i==2
		LPmat = LPvec;
	else
		LPmat = [LPmat LPvec];
	end
end

datetime_ref = datetime(LPmat(1,:));
datetime_DD = datetime(LPmat(2,:));
daysonlydiffDD = floor(days(datetime(LPmat(2,:)) - datetime_ref));

for i=LPstart:size(LightPulseInfo{1},1)
	
	daysdiff = days(datetime(LPmat(i,:)) - datetime_ref);
	daysonlydiff = floor(daysdiff);
	hoursdiff = (daysdiff-floor(daysdiff))*24;
	dhdiff = [daysonlydiff hoursdiff];
	if i==LPstart
		dhdiff_all = dhdiff;
	else
		dhdiff_all = [dhdiff_all dhdiff];
	end
end
LightInfoOut = [daysonlydiffDD dhdiff_all];
