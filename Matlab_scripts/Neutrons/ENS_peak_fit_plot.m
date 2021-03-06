function ENS_peak_fit_plot(hh0,I,hCenter,fitresult,varargin)
%  Plot a fit.
%
%  Data for 'ICpV' (Ikeda-Carpenter-pseudo-Voigt) fit:
%      X Input (numeric array): hh0, cut of data along hh0 direction
%      Y Input (numeric array): I, neutrons intensity received by detector, in arb. units
%      hCenter (integer): center of peak, in reciprocal space units 
%  Output:
%      plot of fit
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 16-Feb-2019 19:01:29


%% Plot fit result
[xData, yData] = prepareCurveData(hh0,I);
if nargin>4
    dataExcl = varargin{1};
%     dataExcl = hh0<hCenter-.15 | hh0>hCenter+varargin{1};
else dataExcl = hh0<hCenter-.15 | hh0>hCenter+0.2;
end
excludedPoints = excludedata( xData, yData, 'Indices', dataExcl );

% Plot fit with data.
figure
h = plot( fitresult, xData, yData, excludedPoints );
legend( h, 'I vs. hh0', 'Excluded I vs. hh0', 'fit ICpV', 'Location', 'NorthEast' );
% Label axes
xlabel("hh0")
ylabel("I (a.u.)")
grid on

end