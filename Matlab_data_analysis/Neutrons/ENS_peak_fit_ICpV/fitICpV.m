classdef fitICpV
    properties 
        X; Y;% X and Y data
        hc;% peak center
        dataExcl;% logical that determines which data points to exclude from fit
        I = 2e5;% first fit parameters startpoint
        R = 0;
        alpha = 140;
        beta = 0;
        gamma = 1e-3;
        sigma = 6.6e-3;
        k = 0.05;% last fit parameter
    end
    methods
        function obj = fitICpV(X, Y, hCenter)
            obj.X = X;
            obj.Y = Y;
            obj.hc = hCenter;
            obj.dataExcl = obj.X<obj.hc-.15 | obj.X>obj.hc+0.2; 
        end
        function [fitresult, gof] = compute_fit(obj,varargin)
%  Create a fit.
%  Data for 'ICpV' (Ikeda-Carpenter-pseudo-Voigt) fit:
%      X Input (numeric array), e.g. hh0, cut of data along hh0 direction
%      Y Input (numeric array), e.g. I, neutrons intensity received by detector, in arb. units
%      hCenter (integer): position (ideally center) of peak, in reciprocal space units 
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
            %% Identify free fit parameters
            function fitstr = fitEqStr(idx,varargin)
% Construct string defining the fit equation, depending on how many
% free parameters the user inputs in addition to intensity and position
% idx is the index of the peak (i.e. of the fit), in case there are several
                varname = {'R','alpha','beta','gamma','sigma'};
% Possible free fit parameters from user input
% They are listed in the same order as specified in the definition of voigtIkedaCarpenter_ord
                infname = cell(1,5);
                fitstr = sprintf("I%i*voigtIkedaCarpenter_ord(x,[",idx);% beginning of string defining the fit equation
                for j=1:length(varname)
% loop through each possible additional free parameter, in the order of varname
                    if any(strcmp(varname{j},varargin))
                        % if the user actually input its name
                        infname{j} = strcat(varname{j},sprintf("%i,",idx));
                        % format it to include it in fitstr, with the fit index
                    else infname{j} = sprintf("%d%i,",obj.(varname{j}),idx);
% otherwise just use the value of the corresponding object property to constrain the fit 
                    end
                    fitstr = strcat(fitstr, string(infname{j}));
% add both the names of the free parameters and the constrained values to the fit equation
                end
                fitstr = strcat(fitstr,sprintf("%d%i,x0])",obj.k,idx));
                % end of string defining the fit equation
            end
%% Create arrays containing lower bounds and initial values of fit parameters
            lowBounds = zeros(length(varargin)+2,1);% Most fit parameters have zero as lower bound
            lowBounds(end) = -Inf;% except hc which can be negative
            initParam = ones(length(varargin)+2,1);
% Initilize array with length that depends on how many free parameters the
% use wants in addition to intensity and position
            initParam(1) = obj.I;% First free parameter is intensity
            initParam(end) = obj.hc;% Last free parameter is position
            for i=1:length(varargin)% If the user inputs other free parameters
                initParam(i+1) = obj.(varargin{i});% initiliaze them
            end
%% Perform fit
            [xData, yData] = prepareCurveData(obj.X,obj.Y);
            % Set up fittype and options.
            ft = fittype( fitstr,'independent', 'x', 'dependent', 'y' );
            excludedPoints = excludedata( xData, yData, 'Indices', obj.dataExcl );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = lowBounds;% use the above defined arrays for lower bounds
            opts.StartPoint = initParam;% and starting parameter values
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