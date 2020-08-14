% Northwestern style for actogram plot
% 4:3 = width:height in 78 days in ClockLab

function Actogram2xOuterNorthWesternStyle(x, NstepsPerHr, T, Tcycle, ll, bb, wx, wh)

graycolor=0;

x=[x;zeros(Tcycle*NstepsPerHr*4,1)];

days0=ceil(T/Tcycle+1);
daylength0=Tcycle*NstepsPerHr;
t0=(1:2*daylength0)'/NstepsPerHr;

seventyeightdays_width = 300;
oneday_height = seventyeightdays_width*(3./4.)/78.;
totalday_height = oneday_height*days0;

xp = 0;
yp = 0;

fig0=figure;set(fig0, 'OuterPosition',[100 50 800 1000]);

for day=1:(days0)

   	today=x((1+(day-1)*daylength0):((day+1)*daylength0));
    today=day-0.9*today/max(today);

    bar_handle=bar(t0(today>0),today(today>0),'FaceColor',graycolor*[1 1 1],'BaseValue',day,'EdgeColor',graycolor*[1 1 1]);
    baseline_handle=get(bar_handle,'BaseLine');
    set(baseline_handle,'Color','w');

    hold on;
	
	if (T>150)
    	plot([0.1 47.9],[day day]+1/2,'w','LineWidth',3);
	else
    	plot([0.1 47.9],[day day]+1/2,'w','LineWidth',6);
	end	

end

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left+ll bottom+bb ax_width+wx ax_height+wh];

xticks([0:Tcycle:2*Tcycle]);

xlh = xlabel('Time (hour)');
ylh = ylabel('Time (day)');

axis([0 2*Tcycle 0 days0-1]);
axis ij;	

xa = 1;
ya = 3/4/78*xa*days0;

set(gca,'Fontsize',20,'PlotBoxAspectRatio',[xa ya 1]);
