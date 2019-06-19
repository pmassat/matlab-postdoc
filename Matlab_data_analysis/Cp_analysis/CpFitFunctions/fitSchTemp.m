function [fitresult, gof] = fitSchTemp(T, Cp, wghts, maxtfit)
%CREATEFITS(TFIT,CPFIT)
%  Create fits.
%
%  Data for 'untitled fit 1' fit:
%      X Input : Tfit
%      Y Output: Cpfit
%  Data for 'untitled fit 2' fit:
%      X Input : Tfit
%      Y Output: Cpfit
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 18-Jun-2019 16:51:40

%% Initialization.

% Initialize arrays to store fits and goodness-of-fit.
% fitresult = cell( 2, 1 );
% gof = struct( 'sse', cell( 2, 1 ), ...
%     'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: 'untitled fit 1'.
[xData, yData, weights] = prepareCurveData( T, Cp, wghts );

% Set up fittype and options.
ft = fittype( 'fSchTemp(T,1,D)+c', 'independent', 'T', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Domain', [0 maxtfit] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [5 .008];
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
pdat = errorbar(xData,yData,1./wghts,'.b','MarkerSize',18,'LineWidth',2);
pexcl = plot(xData(excludedPoints),yData(excludedPoints),'xk',...
    'MarkerSize',12);
legend([pdat,pexcl,pfit],'Data','Excluded','Schottky fit','Location','northwest');
% Label axes
xlabel('$T$ (K)')
ylabel('$C_p/R$')
grid on


