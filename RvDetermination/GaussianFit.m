function [out1,out2] = GaussianFit(xdata,ydata,startpoints,plotflag,color,lw)

% temp = xdata';
% a = polyfit(temp([1 end]),ydata([1 end]),1);
% 
% ydata = ydata./polyval(a,temp);

% 
ydata = ydata./max(ydata);

% [~,ind] = max(ydata);
% % centind = xdata==0;
% 
% 
% [~,uind] = min(abs(xdata-5)); %index of closest value
% 
% dind = abs(uind-ind);
% 
% xdata2 = xdata(ind-dind:ind+dind);
% ydata2 = ydata(ind-dind:ind+dind);

% xdata2 =xdata;
% ydata2=ydata;

% ydata2 = detrend(ydata2);
% ydata2 = ydata2+abs(min(ydata2));
% ydata2 = ydata2+(max(ydata)-max(ydata2));
startpoints(1) = min(ydata);
startpoints(4) = max(ydata);
gaussian = 'd1-a1 * exp(-0.5.*((x-b1).^2) / c1.^2) ';


[myfit,gof] = fit(xdata',ydata,gaussian,'Start',startpoints);

out1 = gof.adjrsquare;
params = coeffvalues(myfit);
out2 = params(2);

plotx = xdata(1):(xdata(2)-xdata(1))/100:xdata(end);

if plotflag == 1
    
    figure(1234)
    hold on
    plot(xdata,ydata,'o','color','k')
    plot(plotx,myfit(plotx),'-','color',color,'linewidth',lw)
    box on
    grid on
end