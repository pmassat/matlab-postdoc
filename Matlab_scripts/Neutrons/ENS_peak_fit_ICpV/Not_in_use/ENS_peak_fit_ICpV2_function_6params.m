function [fitresult, gof] = ENS_peak_fit_ICpV2_function_6params(hh0, I, hCenter, dataExcl, varargin)
%CREATEFIT(HH067,I67)
%  Create a fit.
%
%  Data for 'ICpV' (Ikeda-Carpenter-pseudo-Voigt) fit:
%      X Input (numeric array): hh0, cut of data along hh0 direction
%      Y Input (numeric array): I, neutrons intensity received by detector, in arb. units
%      hCenter (integer): center of peak, in reciprocal space units 
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 16-Feb-2019 19:01:29


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData(hh0,I);
peak1InitCenter = hCenter-0.03;
peak2InitCenter = hCenter*1.0;

% I1*voigtIkedaCarpenter(x,[140,gamma1,6.6e-3,1,0,0.05,x01])+I2*voigtIkedaCarpenter(x,[140,gamma2,6.6e-3,1,0,0.05,x02])
% Set up fittype and options.
ft = fittype( ['I1*voigtIkedaCarpenter_ord(x,[0,140,0,gamma1,6.6e-3,0.05,x01])+',...
    'I2*voigtIkedaCarpenter_ord(x,[0,140,0,gamma2,6.6e-3,0.05,x02])'],...
    'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', dataExcl );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0 -Inf -Inf];
opts.StartPoint = [2e5 2e5 1e-3 1e-3 peak1InitCenter peak2InitCenter];
% the choice of initial parameters is critical to the convergence of the fit
% *NOTE*: Matlab orders the fit parameters by alphabetical order and puts
% capital letters before small ones
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

c1 = @(x) fitresult.I1*voigtIkedaCarpenter_ord(x,[0,140,0,fitresult.gamma1,6.6e-3,0.05,fitresult.x01]);
c2 = @(x) fitresult.I2.*voigtIkedaCarpenter_ord(x,[0,140,0,fitresult.gamma2,6.6e-3,0.05,fitresult.x02]);

if nargin>4
    % Plot fit with data.
    figure
    p1 = fplot(c1,[hCenter-0.2 hCenter+0.2],'-g','LineWidth',2);
    hold on
    p2 = fplot(c2,[hCenter-0.2 hCenter+0.2],'-m','LineWidth',2);
    h = plot( fitresult, xData, yData, excludedPoints, '.k');
%     h = fplot(@(x)c1(x)+c2(x),[hCenter-1. hCenter+1.],'-r','LineWidth',2);
    legend( h, 'I vs. hh0', 'Excluded I vs. hh0', 'fit ICpV', 'Location', 'NorthEast' );
    % Label axes
    xlabel("hh0")
    ylabel("I (a.u.)")
    grid on
end
