% Fitting onsets after wrapping them to 24 hours
% Input
%	Xbo, Ybo: onset time and day before the light pulse
%	Xao, Yao: onset time and day after the light pulse
%	LightPulseDay, LightPulseTime: day and time of the light pulse
%	Tcycle: 24 hours
%	BoolShowPlot: 1 for plot, 0 for without plot
% Output
%	Xbo, Ybo, Xao, Yao: Onsets after wrapping
%	pbo, pao: Slope and intercept

%	FitSSEbo, FitSSEao: Sum of squares due to error (SSE)
%	RMSEbo, RMSEao: Root mean square error (RMSE) values
%	Rsquaredbo, Rsquaredao: R-squared values
%	Fstatbo, Fstatao: F-statistic values
% 		of a line fitting onsets before and after the light pulse

%	PhaseAtLightPulse: Phase at the light pulse
%	PhaseShift: Phase shift due to the light pulse

function [Xbo, Ybo, Xao, Yao, pbo, FitSSEbo, Rsquaredbo, Fstatbo, RMSEbo, pao, FitSSEao, Rsquaredao, Fstatao, RMSEao, PhaseAtLightPulse, PhaseShift] = FittingOnsetsWrapped(Xbo, Ybo, Xao, Yao, LightPulseDay, LightPulseTime, Tcycle, BoolShowPlot)

	XYbo = [Xbo Ybo];
	XYao = [Xao Yao];

	XYbo(:,1) = XYbo(:,1) + Tcycle;
	XYbo(:,2) = XYbo(:,2) - 1;
	XYao(:,1) = XYao(:,1) + Tcycle;
	XYao(:,2) = XYao(:,2) - 1;

	[Xbo, Ybo] = AlignXY(XYbo(:,1), XYbo(:,2), Tcycle);
	[Xao, Yao] = AlignXY(XYao(:,1), XYao(:,2), Tcycle);

	[pbo, FitSSEbo, Rsquaredbo, Fstatbo, RMSEbo] = FittingOnsets(Xbo, Ybo, 'ro', 'c', BoolShowPlot);
	[pao, FitSSEao, Rsquaredao, Fstatao, RMSEao] = FittingOnsets(Xao, Yao, 'ro', 'c', BoolShowPlot);

	[PhaseAtLightPulse, PhaseShift] = PhaseShiftOnset (pbo, pao, LightPulseDay, LightPulseTime, Tcycle);

function [p, FitSSE, R_squared, Fstat, FitRMSE] = FittingOnsets(X, Y, DotType, FitLineColor, BoolShowPlot)

	% Fitting onsets to a line: y ~ 1 + x1
	% y: X (time), x1: Y (day)
	% 	** xy axis is reversed.
	mdl = fitlm(Y, X,'RobustOpts', 'on');
	mRR = mdl.Residuals.Raw;

	FitSSE = mdl.SSE;								% Sum of squares due to error (SSE)
	mdl_slope = mdl.Coefficients.Estimate(2);		% slope the fitting line
	mdl_intercept = mdl.Coefficients.Estimate(1);	% y-intercept of the fitting line
	R_squared = mdl.Rsquared.Adjusted;				% R-squared
	%p_value = mdl.Coefficients.pValue;				% p-value
	AnovaTable = anova(mdl);
	FstatAll = AnovaTable.F;
	Fstat = FstatAll(1);
	p_value = AnovaTable.pValue;
	%[pv, Fstat] = mdla.coefTest
	FitRMSE = mdl.RMSE;
		
	p = [mdl_intercept, mdl_slope];

	XFit = mdl_slope*Y + mdl_intercept;
	minXFit = min(XFit);
	maxXFit = max(XFit);

	XFitInterval = minXFit:0.02:maxXFit;
	YFitInterval = (XFitInterval - mdl_intercept)/mdl_slope;

	if (BoolShowPlot==1)
		% plot fitting lines
		% plot(X,Y, 'ko');
		% hold on;
		plot(XFitInterval, YFitInterval, 'g', 'LineWidth', 3);
		% axis ij
	end


function [Xo, Yo] = AlignXY(Xo, Yo, Tcycle)

	mX = median(Xo);

	while (sum(abs(mX - Xo) > Tcycle/2.) > 0)
		idx = find(mX - Xo > Tcycle/2.);
		Xo(idx) = Xo(idx) + Tcycle;
		Yo(idx) = Yo(idx) - 1;

		idx = find(Xo - mX > Tcycle/2.);
		Xo(idx) = Xo(idx) - Tcycle;
		Yo(idx) = Yo(idx) + 1;

		mX = median(Xo);
	end
	
	if (mX > Tcycle + Tcycle/2.)
		Xo = Xo - Tcycle;
		Yo = Yo + 1;
	end	


function [PhaseAtLightPulse, PhaseShift] = PhaseShiftOnset (pb, pa, LightPulseDay, LightPulseTime, Tcycle)

	mb = pb(2);
	cb = pb(1);
	ma = pa(2);
	ca = pa(1);

	yL = LightPulseDay;
	xL = LightPulseTime;

	xb = mb*yL + cb;
	xa = ma*yL + ca;

	if xb > Tcycle
		xb = xb - Tcycle;
	elseif xa > Tcycle
		xa = xa - Tcycle;
	end

	if (xa - xb > Tcycle/2.)
		xa = xa - Tcycle;
	elseif (xa - xb < -Tcycle/2.)
		xa = xa + Tcycle;
	end

	if (xL - xb > Tcycle/2.)
		xL = xL - Tcycle;
	elseif (xL - xb < -Tcycle/2.)
		xL = xL + Tcycle;
	end

	dxbL = xb - xL;
	mb = mb + Tcycle;

	% Phase at the light pulse
	if (dxbL >= 0)
		PhaseAtLightPulse = (mb - dxbL);
	else
		PhaseAtLightPulse = -dxbL;
	end
	PhaseAtLightPulse = PhaseAtLightPulse*(Tcycle/mb);

	if PhaseAtLightPulse > Tcycle/2.
		PhaseAtLightPulse = PhaseAtLightPulse - Tcycle;
	end

	PhaseAtLightPulse = PhaseAtLightPulse + Tcycle/2.;
	PhaseAtLightPulse = mod(PhaseAtLightPulse, Tcycle);

	% Phase shift
	PhaseShift = -(xa - xb);
	PhaseShift = PhaseShift*(Tcycle/mb);

	if (PhaseShift > 12)
		PhaseShift = PhaseShift - Tcycle;
	elseif (PhaseShift < -12)
		PhaseShift = PhaseShift + Tcycle;
	end



