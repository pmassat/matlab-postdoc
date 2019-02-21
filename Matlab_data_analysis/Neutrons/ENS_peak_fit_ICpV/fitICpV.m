classdef fitICpV < handle
    properties 
        X; Y;% X and Y data
        dataExcl;% logical that determines which data points to exclude from fit
        I = 2e5;% default value for the first fit parameter, if constrained (see function fitEqStr)
% Defining these fit parameters as object properties allows to call
% them both for fitting and plotting. The user must be careful however not
% to change them between fitting and plotting, otherwise the plotted fit
% will not correspond to the calculated one
        R = 0;
        alpha = 140;
        beta = 0;
        gamma = 1e-3;
        sigma = 6.6e-3;
        k = 0.05;
        x0;% last fit parameter; 
% Note that x0 is also the estimated peak position and is therefore also
% used to determine the default range of excluded data points (see constructor)
        freeParams; % cell array of cell arrays, each sub-cell array containing
% the names of the free parameters for each fit. By default, it is assumed
% that there is only one fit, i.e. one peak, with 2 free parameters:
% intensity and position
    end
    methods
        function obj = fitICpV(X, Y, xCenter)
%  Create a fit.
%  Data for 'ICpV' (Ikeda-Carpenter-pseudo-Voigt) fit:
%      X Input (numeric array), e.g. hh0, cut of data along hh0 direction
%      Y Input (numeric array), e.g. I, neutrons intensity received by detector, in arb. units
%      hCenter (integer): position (ideally center) of peak, in reciprocal space units 
            obj.X = X;
            obj.Y = Y;
            obj.x0 = xCenter;
            obj.dataExcl = obj.X<obj.x0-.15 | obj.X>obj.x0+0.2; 
            obj.freeParams = {{'I',obj.I;'x0',obj.x0}};
        end
        function [fitresult, gof] = compute_fit(obj,varargin)
% Input: in addition to data defined in object (see constructor)
%       Additional arguments input by the user in varargin should be
%       cell arrays, each containing the names and initial values (StartPoint, see fit execution code)
%       of the free parameters for each peak to be fitted.
%       The user can input as many cell arrays as they want, each
%       corresponding to one peak.
%       Example: varargin = {'gamma',1e-3;'alpha',200;'x0',-8.1},{'gamma',1e-3;'x0',-7.9}
%       This defines free fit parameters for two peaks: 3 parameters for
%       the first one, two for the second one
%       Note: the parameters in each cell array need NOT be in the
%       order defined in the fit parameter function
% Output:
%       fitresult : a fit object representing the fit.
%       gof : structure with goodness-of fit info.
            if nargin>1
                obj.freeParams = varargin;
            end
            nPeaks = length(obj.freeParams);
%             eqParamName = {'I','R','alpha','beta','gamma','sigma','k','x0'};
            eqParamName = {'I';'R';'alpha';'beta';'gamma';'sigma';'k';'x0'};
            % Generic names of free parameters for the Ikeda-Carpenter-pseudo-Voigt fit
% They are listed in the same order as specified in the definition of voigtIkedaCarpenter_ord
%             % Create list of parameter names for each peak, consisting of
%             % the names stored in the above 'eqParamName' cell array, each
%             % of them being concatenated with a peak number
            totalNumFreeParams = 0;% count the total number of free parameters
            freeFitParams = repmat(obj.freeParams,1);
            for j1=1:nPeaks% j1 is the peak number
                if size(obj.freeParams{j1},1)==1
                    warning("Fit is unlikely to converge with only one "+...
                        "free parameter; it is advised to input at least two, "+...
                        "or none, in which case peak intensity 'I' and position 'x0' "+...
                        "will be used.")
                end
                totalNumFreeParams = totalNumFreeParams + size(obj.freeParams{j1},1);
                % the number of free parameters input by the user for
                % peak #j1 equals the number of rows in 'obj.freeParams{j1}'
                if size(obj.freeParams{j1},2)~=2
                    error(strcat("Error: The array of free parameters for each peak ",...
                        "should be formatted as a N x 2 cell array, ",...
                        "each row containing the name and initial value ",...
                        "of a parameter, ",...
                        "e.g. {'gamma',1e-3;'alpha',200;",...
                        "'x0',-8.1},{'gamma',1e-3;'x0',-7.9}"));
                    break
                end
%                 freeFitParams{j1} = cell(length(eqParamName),2);
%                 % for each peak, create a cell array with 2 columns and as
%                 % many rows as there are parameters
                for k1=1:size(obj.freeParams{j1},1)
%                     peaksParam{j1}{k1} = cell(1,2);
                    if ~any(strcmp(obj.freeParams{j1}{k1,1},eqParamName),'all')
                        error(['Error: input parameter names must be one of the following: [' +...
                            char(join(string(eqParamName),', ')) +...
                            '] entered as character vector.'])
                    end
                    if ~ischar(obj.freeParams{j1}{k1,1})
                       error('Error. \nInput parameter names must be char, not %s.',...
                           class(obj.freeParams{j1}{k1,1}))
                    end
                    obj.freeParams{j1}{k1,1} = strcat(obj.freeParams{j1}{k1,1},int2str(j1));
%                     freeFitParams{j1}{k1,1} = strcat(obj.freeParams{j1}{k1,1},sprintf("%i",j1));
                    % rename the free parameter to include the peak number
%                     obj.freeParams{j1}{k1,2} = obj.(eqParamName{k1});
%                     % the second column contains the default value of the parameter,
                end
            end
            
            %% Create full fit equation string 
% depending on the number of peaks and the number of free fit parameters
% for each peak
            function fitstr = fitEqStr(obj,paramName,peakIndex)
% Construct string defining the fit equation for a single peak, depending on how many
% free parameters the object contains, as defined in obj.freeParams
                function inFName = inFuncName(obj,varName,pkIndx)
% This function determines if string varName is contained in cell array obj.freeParams{pkIndx},
% which should be a N x 2 array containing parameter names as strings in the first column, 
% and the associated numerical value in the second column, even though this
% second column is not used here
%                     paramsCellStr = string(obj.freeParams{pkIndx});
                    % conversion to string array guaranties that all
                    % elements of obj.freeParams{pkIndx} will be taken into account
                    % when performing the following 'strcmp' test
                    if find(contains(obj.freeParams{pkIndx}(:,1),varName))
% argument 'all' is needed for function any to test on all elements and
% return a logical scalar, otherwise it returns a logical array, which does
% not work with 'if' statement
                        idx = contains(obj.freeParams{pkIndx}(:,1),varName);
                    % if string varName is contained in cell array obj.freeParams{pkIndx}
                        inFName = obj.freeParams{pkIndx}{idx};
%                         inFName = strcat(varName,sprintf("%i",pkIndx));
%                         idx = strcmp(varName,paramsCellStr);% index of string
%                         inFName = paramsCellStr(idx,1);
                    % use it as is in the string used to define the fit in fittype
                    else; inFName = string(obj.(varName));
                    % otherwise use the associated numerical value defined
                    % as object property
                    end
                end
                infname = inFuncName(obj,paramName{1},peakIndex);
                fitstr = strcat(infname,"*voigtIkedaCarpenter_ord(x,[");
                % beginning of string defining the fit equation
                for j2=2:length(paramName)-1
% loop through each possible additional free parameter, in the order of eqParamName
                    infname = inFuncName(obj,paramName{j2},peakIndex);
                    fitstr = strcat(fitstr,infname,",");
% depending whether a parameter should be free or constrained in the fit,
% add its name or its value, respectively
                end
                infname = inFuncName(obj,paramName{end},peakIndex);
                fitstr = strcat(fitstr,string(infname),"])");
                % end of string defining the fit equation
            end
            
            fullFitStr = fitEqStr(obj,eqParamName,1);% equation string for first peak
            if nPeaks>1
                for i1=2:nPeaks
                    fullFitStr = fullFitStr + " + " + ...% sum equation strings for all peaks
                        fitEqStr(obj,eqParamName,i1);
                end
            end
            fullFitStr
%% Create arrays containing lower bounds and initial values of fit parameters
            lowBounds = zeros(1,totalNumFreeParams);%
            initParams = ones(1,totalNumFreeParams);%
            % initialize arrays containing all the lower bounds and initial
            % values of free fit parameters, so that their length is sum(Ni)
            cntr = 0;% counter that increments at each iteration of the following loop
            for j3=1:nPeaks
                nPrm = size(obj.freeParams{j3},1);% number of free parameters for peak #j3
                lowBnds = zeros(1,nPrm);% Most fit parameters have zero as lower bound
%                     if any(strcmp("x0",string(obj.freeParams{j3})),'all')
                    if find(contains(obj.freeParams{j3}(:,1),'x0'))
                        lowBnds(end) = -Inf;% except x0 which can be negative                        
                    end
                initPrm = ones(1,nPrm);
% Initilize array with length that depends on how many free parameters the user input
                for i=1:nPrm% for each free parameters
                    initPrm(i) = obj.freeParams{j3}{i,2};% initiliaze them using the user input value
                end
                lowBounds(cntr+1:cntr+nPrm) = lowBnds;% store lower bounds values
                initParams(cntr+1:cntr+nPrm) = initPrm;% store initial values
                cntr = cntr + nPrm;% increase cntr for next iteration
            end

%% Perform fit
            [xData, yData] = prepareCurveData(obj.X,obj.Y);
            % Set up fittype and options.
            ft = fittype( fullFitStr,'independent', 'x', 'dependent', 'y' );
            excludedPoints = excludedata( xData, yData, 'Indices', obj.dataExcl );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = lowBounds;% use the above defined arrays for lower bounds
            opts.StartPoint = initParams;% and starting parameter values
% the choice of initial parameters is critical to the convergence of the fit
            opts.Exclude = excludedPoints;
            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
        end
        
        
%% Plot fit
        function plot_fit(obj,fitresult,varargin)
%  Data for 'ICpV' (Ikeda-Carpenter-pseudo-Voigt) fit:
%      X Input (numeric array): hh0, cut of data along hh0 direction
%      Y Input (numeric array): I, neutrons intensity received by detector, in arb. units
%      hCenter (integer): center of peak, in reciprocal space units 
%  Output:
%      plot of fit
            
            [xData, yData] = prepareCurveData(obj.X,obj.Y);
            excludedPoints = excludedata( xData, yData, 'Indices', obj.dataExcl );
            c2 = @(x) fitresult.I2.*voigtIkedaCarpenter_ord(x,[0,140,0,fitresult.gamma2,6.6e-3,0.05,fitresult.x02]);
            % Plot fit with data.
            figure
            h = plot( fitresult, xData, yData, excludedPoints );
            legend( h, 'I vs. hh0', 'Excluded I vs. hh0', 'fit ICpV', 'Location', 'NorthEast' );
            % Label axes
            xlabel("hh0")
            ylabel("I (a.u.)")
            grid on
        end
        
    end
    
end