classdef fitICpV
    properties 
        hh0;% X data: hh0 cut in reciprocal space
        I;% Y data
        hc;% peak center
        dataExcl;% logical that determines which data points to exclude from fit
        R = 0;% first fit parameters
        alpha = 140;
        beta = 0;
        gamma = 0;
        sigma = 6.6e-3;
        k = 0.05;% last fit parameter
    end
    methods
        function obj = fitICpV(X, Y, hCenter)
            obj.hh0 = X;
            obj.I = Y;
            obj.hc = hCenter;
            obj.dataExcl = obj.hh0<obj.hc-.15 | obj.hh0>obj.hc+0.2; 
        end
        function [fitresult, gof] = compute_fit(obj,varargin)
            %  Create a fit.
            %  Data for 'ICpV' (Ikeda-Carpenter-pseudo-Voigt) fit:
            %      X Input (numeric array): hh0, cut of data along hh0 direction
            %      Y Input (numeric array): I, neutrons intensity received by detector, in arb. units
            %      hCenter (integer): center of peak, in reciprocal space units 
            %  Output:
            %      fitresult : a fit object representing the fit.
            %      gof : structure with goodness-of fit info.
            %% Initialize fit parameters
            varname = {'R','alpha','beta','gamma','sigma'};
            infname = cell(1,5);
            fitstr = 'I*voigtIkedaCarpenter_ord(x,[';
            for j=1:length(varname)
                if any(strcmp(varname{j},varargin))
                    infname{j} = strcat(varname{j},',');
                else infname{j} = sprintf("%d,",obj.(varname{j}));
                end
                fitstr = strcat(fitstr, string(infname{j}));
            end
            fitstr = strcat(fitstr,sprintf("%d,x0])",obj.k));
            %% Fit: 'untitled fit 1'.
            [xData, yData] = prepareCurveData(obj.hh0,obj.I);
            % Set up fittype and options.
%             'I*voigtIkedaCarpenter_ord(x,[0,140,0,gamma,6.6e-3,0.05,x0])'
%                 '%d,%d,%d,%d,%d,%d,x0])';
            ft = fittype( fitstr,'independent', 'x', 'dependent', 'y' );
            excludedPoints = excludedata( xData, yData, 'Indices', obj.dataExcl );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [0 0 -Inf];
            opts.StartPoint = [5e5 1e-3 obj.hc];
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
            [xData, yData] = prepareCurveData(obj.hh0,obj.I);
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
    end
end