function [fitresult, gof] = fitCpmfvsTemp(T, Cp, wghts, Tc)
%CREATEFIT(T0,CP0,WGTHS)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : T0
%      Y Output: Cp0
%      Weights : wgths
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 15-Mar-2019 17:27:39


%% Fit: 'untitled fit 1'.
[xData, yData, weights] = prepareCurveData( T, Cp, wghts );

% Set up fittype and options.
ft = fittype( sprintf('a*Cpmf_vs_temp(x/%d)',Tc), 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', 29:44 );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.485375648722841;
opts.Weights = weights;
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
Xfit = linspace(0,max(T),1000);
Yfit = fitresult(Xfit);% compute fit over a controlled number of points
figure;
hold on
pfit = plot(Xfit,Yfit,'r-');
pdat = errorbar(xData,yData,1./wghts,'xb','MarkerSize',12,'LineWidth',2);
pexcl = plot(xData(excludedPoints),yData(excludedPoints),'xg',...
    'MarkerSize',12);
legend([pdat,pexcl,pfit],'Cp(H=0) vs T','Excluded','MF fit');
% Label axes
xlabel T(K); ylabel Cp(uJ/K); grid on

