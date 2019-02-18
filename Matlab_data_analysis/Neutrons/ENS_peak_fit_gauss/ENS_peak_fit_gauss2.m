function [fitresult, gof] = ENS_peak_fit_gauss2(hh0, I, hCenter)
%CREATEFIT(HH0,I)
%  Create a fit.
%
%  Data for 'Gaussian 2' fit:
%      X Input (numeric array): hh0, cut of data along hh0 direction
%      Y Input (numeric array): I, neutrons intensity received by detector, in arb. units
%      hCenter (integer): center of peak, in reciprocal space units 
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 14-Feb-2019 15:54:35


%% Fit: '2 Gaussians'.
[xData, yData] = prepareCurveData( hh0, I );
dataExcl = hh0<hCenter-.15 | hh0>hCenter+0.2;
peak1InitCenter = hCenter*1.0025;
peak2InitCenter = hCenter*0.9975;

% Set up fittype and options.
ft = fittype( 'gauss2' );
excludedPoints = excludedata( xData, yData, 'Indices', dataExcl );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
% opts.Upper = [3.5e6 Inf Inf Inf Inf Inf];
opts.StartPoint = [6e6 peak1InitCenter 0.01 3e6 peak2InitCenter 0.015];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Calculate individual gaussians
g1 = fitresult.a1.*exp(-((xData-fitresult.b1)/fitresult.c1).^2);
g2 = fitresult.a2.*exp(-((xData-fitresult.b2)/fitresult.c2).^2);

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
p1 = plot(xData,g1,'-g');
hold on
p2 = plot(xData,g2,'-m');
h = plot( fitresult, xData, yData, excludedPoints, '.k');
legend( h, 'I vs. hh0', 'Excluded I vs. hh0', 'fit 2 Gaussians', 'Location', 'NorthEast' );
% Label axes
xlabel("hh0")
ylabel("I (a.u.)")
grid on


