classdef fitICpV < handle
    properties 
        X; Y;% X and Y data
        dataExcl;% logical that determines which data points to exclude from fit
%         I = 2e5;% default value for the first fit parameter, if constrained (see function fitEqStr)
% % Defining these fit parameters as object properties allows to call
% % them both for fitting and plotting. The user must be careful however not
% % to change them between fitting and plotting, otherwise the plotted fit
% % will not correspond to the calculated one
%         R = 0;
%         alpha = 140;
%         beta = 0;
%         gamma = 1e-3;
%         sigma = 6.6e-3;
%         k = 0.05;
%         x0;% last fit parameter; 
% Note that x0 is also the estimated peak position and is therefore also
% used to determine the default range of excluded data points (see constructor)
        allParams;
        freeParams; % cell array of cell arrays, each sub-cell array containing
% the names of the free parameters for each fit. By default, it is assumed
% that there is only one fit, i.e. one peak, with 2 free parameters:
% intensity and position
    end
    methods
%% Class constructor
        function obj = fitICpV(X, Y, xPeaks)
% Create a fit.
% Data for 'ICpV' (Ikeda-Carpenter-pseudo-Voigt) fit:
%      X Input (numeric array), e.g. hh0, cut of data along hh0 direction
%      Y Input (numeric array), e.g. I, neutrons intensity received by detector, in arb. units
%      xPeak (integer): position of peak, in reciprocal space units 
            if ~isa(X,'double') || ~isa(Y,'double') || ~isa(xPeaks,'double')
               error('Input parameters must all be double.')                
            end
            obj.X = X; obj.Y = Y;
            obj.dataExcl = obj.X<min(xPeaks)-.15 | obj.X>max(xPeaks)+0.2; 
            nPeaks = length(xPeaks);
            fixedKeySet = cell(1,nPeaks); fixedValueSet = cell(1,nPeaks); 
            obj.allParams = cell(1,nPeaks); obj.freeParams = cell(1,nPeaks);
            for i0=1:length(xPeaks)
                fixedKeySet{i0} = {'I','R','alpha','beta','gamma','sigma','k','x0'};
                fixedValueSet{i0} = [2e5, 0, 140, 0, 1e-3, 6.6e-3, 0.05, xPeaks(i0)];
                obj.allParams{i0} = containers.Map(fixedKeySet{i0},fixedValueSet{i0});
                obj.freeParams{i0} = {};%containers.Map('KeyType','char','ValueType','double');
            end
        end
        
%% Compute fit method
        function [fitresult, gof] = compute_fit(obj)
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
%             if nargin>1
%                 obj.freeParams = varargin;
%                 obj.freeParams = varargin;
%             end
            nPeaks = length(obj.freeParams);% user should input free parameters manually before calling this method
%             eqParamName = {'I';'R';'alpha';'beta';'gamma';'sigma';'k';'x0'};
            % Generic names of free parameters for the Ikeda-Carpenter-pseudo-Voigt fit
% They are listed in the same order as specified in the definition of voigtIkedaCarpenter_ord
            totalNumFreeParams = 0;% count the total number of free parameters
            for j1=1:nPeaks% j1 is the peak number
                if size(obj.freeParams{j1},1)==0
                    error("User should manually input at least one (preferably two) "+...
                        "free parameter(s) per peak in the 'freeParams' "+...
                        "property of the class prior to calling this method.")
                end
                if size(obj.freeParams{j1},1)==1
                    warning("fit is unlikely to converge with only one "+...
                        "free parameter for any given peak; "+...
                        "it is advised to input at least two free parameters for each peak")%, "+...
%                         "or none *IF there is ONLY ONE peak*, in which case "+...
%                         "peak intensity 'I1' and position 'x01' "+...
%                         "will be used.")
                end
                totalNumFreeParams = totalNumFreeParams + size(obj.freeParams{j1},1);
                % the number of free parameters input by the user for
                % peak #j1 equals the number of rows in 'obj.freeParams{j1}'
                if size(obj.freeParams{j1},2)~=2
                    error(strcat("The 'freeParams' property should be a cell array containing one ",...
                        "cell array with N rows and 2 columns for each peak to fit, ",...
                        "the first column containing the names of the free parameters ",...
                        "and the second column containing the corresponding initial values, ",...
                        "e.g. {{'gamma',1e-3;'alpha',200;",...
                        "'x0',-8.1},{'gamma',1e-3;'x0',-7.9}}"));
                    break
                end
                ksj = keys(obj.allParams{j1});
                for k0=1:size(obj.allParams{j1},1)                
                    flag = 0;
                    if find(contains(obj.freeParams{j1}(:,1),ksj{k0}))
                        flag = 1;
                        break;
                    end
                end
                if ~flag
                    error(['Error: input parameter names must be one of the following: [' +...
                        char(join(string(ksj),', ')) +...
                        '] entered as character vector.'])
                end
                for k1=1:size(obj.freeParams{j1},1)
                    if ~ischar(obj.freeParams{j1}{k1,1})
                       error('Input parameter names must be char, not %s.',...
                           class(obj.freeParams{j1}{k1,1}))
                    end
%                     obj.freeParams{j1}{k1,1} = strcat(obj.freeParams{j1}{k1,1},int2str(j1));
%  No need to if manual user input                   % rename each free parameter to include peak number,
%                     % but keep its type as character vector
                end
            end
            
            %% Create full fit equation string... 
% ... depending on the number of peaks and the number of free fit parameters for each peak
            function fitstr = fitEqStr(referenceMap,paramArray)
% Construct string defining the fit equation for a single peak, depending on how many
% free parameters the object contains, as defined in obj.freeParams
                function inFName = inFuncName(refMap,prmNum,prmArray)
% This function determines if string varName is contained in cell array nameArray,
% which should be a N x 2 array containing parameter names as strings in the first column, 
% and the associated numerical value in the second column, even though this
% second column is not used here
                    ks = keys(refMap);
                    varName = ks{prmNum};
                    if find(contains(prmArray(:,1),varName))
% if string varName is contained in the first column of cell array nameArray
                        idx = contains(prmArray(:,1),varName);
                        % index of string varName in the first column of nameArray
                        inFName = prmArray{idx};
                        % use it as is in the fit equation string
                    else; inFName = string(refMap(ks{prmNum}));
                    % otherwise use the associated default numerical value to constrain the fit
                    end
                end
%                 Ks = keys(obj.allParams{pkIndx});
%                 prmVal = obj.allParams{pkIndx}(Ks{1})
                infname = inFuncName(referenceMap,1,paramArray);
                fitstr = strcat(infname,"*voigtIkedaCarpenter_ord(x,[");
                % beginning of string defining the fit equation
                for j2=2:length(referenceMap)-1
% loop through each possible additional free parameter, in the order of eqParamName
                    infname = inFuncName(referenceMap,j2,paramArray);
                    fitstr = strcat(fitstr,infname,",");
% depending whether a parameter should be free or constrained in the fit,
% add its name or its value, respectively (see function 'inFuncName' above)
                end
                infname = inFuncName(referenceMap,length(referenceMap),paramArray);
                fitstr = strcat(fitstr,string(infname),"])");
                % end of string defining the fit equation
            end
            
            fitEqn = cell(1:nPeaks);
            for i1=1:nPeaks
                fitEqn{i1} = fitEqStr(obj.allParams{i1},obj.freeParams{i1});
            end
            fullFitStr = join(string(fitEqn),'+');% sum equation strings for all peaks
            
%% Create arrays containing lower bounds and initial values of fit parameters
            lowBounds = zeros(1,totalNumFreeParams);%
            initParams = ones(1,totalNumFreeParams);%
            % initialize arrays containing all the lower bounds and initial
            % values of free fit parameters, so that their length is sum(Ni)
            cntr = 0;% counter that increments at each iteration of the following loop
            for j3=1:nPeaks
                nPrm = size(obj.freeParams{j3},1);% number of free parameters for peak #j3
                lowBnds = zeros(1,nPrm);% Most fit parameters have zero as lower bound
                    if find(contains(obj.freeParams{j3}(:,1),'x0'))
                        lowBnds(end) = -Inf;% except x0 which can be negative                        
                    end
                initPrm = ones(1,nPrm);
                for i=1:nPrm% for each free parameter of peak #j3
                    initPrm(i) = obj.freeParams{j3}{i,2};% initiliaze it using the user input value
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
            eqParamName = {'I';'R';'alpha';'beta';'gamma';'sigma';'k';'x0'};
            nPks = length(obj.freeParams);
            fitPrms = ones(1,length(eqParamName));
            funCell = cell(1,nPks);
            for ii=1:nPks
                for jj=1:length(eqParamName)
                    if find(contains(obj.freeParams{ii}(:,1),eqParamName{jj}))
                        idx = contains(obj.freeParams{ii}(:,1),eqParamName{jj});
                        % index of parameter namestring in the first column of obj.freeParams{ii}
                        fitPrms(jj) = fitresult.(obj.freeParams{ii}{idx});
                    else; fitPrms(jj) = obj.(eqParamName{jj});
                    end
                end
                funCell{ii} = @(x) fitPrms(1).*...
                    voigtIkedaCarpenter_ord(x,[fitPrms(2),fitPrms(3),...
                    fitPrms(4),fitPrms(5),fitPrms(6),fitPrms(7),fitPrms(8)]);
            end
            minIndex = find(~obj.dataExcl,1,'first');
            maxIndex = find(~obj.dataExcl,1,'last');
            Xfit = linspace(obj.X(minIndex),obj.X(maxIndex),1000);
            Yfit = fitresult(Xfit);
            
            % Plot fit with data.
            figure
            hold on;
            if nPks>1% if there is more than one peak
                psub = cell(1,nPks);%
                for kk=1:nPks% plot fit curve for each peak
                    psub{kk} = fplot(funCell{kk},...
                        [obj.X(minIndex) obj.X(maxIndex)],'LineWidth',2);
                end
            end
            pdat = plot(xData,yData,'xb','MarkerSize',9);
            pexcl = plot(obj.X(excludedPoints),obj.Y(excludedPoints),'xk',...
                'MarkerSize',9);
            pfit = plot(Xfit,Yfit,'r-');
            legend([pdat,pexcl,pfit],'I vs. hh0','Excluded','fit ICpV');
            % Label axes
            xlabel("hh0"); ylabel("I (a.u.)");
            xlim([obj.X(minIndex)-0.05 obj.X(maxIndex)+0.05]);
            grid on
        end
        
    end
    
end