% Notes on fitICpV class:
% Ideally: 
% * Plot function with subplots if there is a double peak plot, and with
% possibility of input different number of free parameters for each peak
% 
% Steps:
% * Finish adapting compute_fit method for single peak fit
% * Adapt plot_fit method for single peak fit
% * Extend both methods for double peak fit. This might imply creating a
% function for initialization of fit parameters

% for j1=1:nPeaks;peaksParam{j1} = cell(length(eqParamName),2);for k1=1:length(eqParamName);peaksParam{j1}{k1,1} = strcat(eqParamName{k1},sprintf("%i",j1));end;end

% 2019-02-20: start over with debugging of function 'inFuncName'

% Clean up code: remove unused lines