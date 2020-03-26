classdef fitICpV < handle% handle allows to modify object properties in the class' methods
    properties 
        X; Y; dY% X and Y data, with error bars on Y
        dataExcl;% logical that determines which data points to exclude from fit
        allParams;% cell array of cell arrays, each sub-cell array containing
% the names and fixed values (if parameters are not free) of all parameters for each peak.
        freeParams;% cell array of cell arrays, each sub-cell array containing
% the names and initial value of the free parameters for each fit.
        indepFreeParams;% cell array of independent free parameters; it is computed from the full set of free parameters
    end
    methods
%% Class constructor
        function obj = fitICpV(X, Y, dY, xPeaks)
% Create a fit.
% Data for 'ICpV' (Ikeda-Carpenter-pseudo-Voigt) fit:
%      X Input (numeric array), e.g. hh0, cut of data along hh0 direction
%      Y Input (numeric array), e.g. I, neutrons intensity received by detector, in arb. units
%      xPeak (integer): position of peak, in reciprocal space units 
            if ~isa(X,'double') || ~isa(Y,'double') || ~isa(xPeaks,'double')
               error('Input parameters must all be double.')                
            end
            obj.X = X(dY>0); obj.Y = Y(dY>0); obj.dY = dY(dY>0);
            obj.dataExcl = obj.X<min(xPeaks)-.15 | obj.X>max(xPeaks)+0.2; 
            nPeaks = length(xPeaks);
            fixedKeySet = cell(1,nPeaks); fixedValueSet = cell(1,nPeaks); 
            obj.allParams = cell(1,nPeaks); obj.freeParams = cell(1,nPeaks);
            for i0=1:length(xPeaks)
                fixedKeySet{i0} = {'I','R','alpha','beta','gamma','k','sigma','x0'};
% Important note: Matlab sorts the keys without asking! They are sorted
% alphabetically, with capital letters grouped together before lowercase ones
                fixedValueSet{i0} = [2e5, 0, 140, 0, 0, 0.05, 6.6e-3, xPeaks(i0)];
% Defining fit parameters in object properties allows to call
% them both for fitting and plotting. The user must be careful however not
% to change them between fitting and plotting, otherwise the plotted fit
% will not correspond to the calculated one
                obj.allParams{i0} = containers.Map(fixedKeySet{i0},fixedValueSet{i0});
                obj.freeParams{i0} = {};%containers.Map('KeyType','char','ValueType','double');
            end
        end
        
        function indepFreePrms(obj)
% Compute array of independent free parameters by removing dependent
% free parameters from the array of free paraemters
            obj.indepFreeParams = vertcat(obj.freeParams{:});
            ifp = obj.indepFreeParams;% shorter variable name, easier to handle
            dependentRow = [];% create a list of row indices of the obj.freeParams that contain dependent fit parameters
            for refRow=1:length(ifp)
                if any(dependentRow==refRow)% if the row index of the row to be compared is already contained in list dependentRow
                    continue% continue to the next row
                end%
                for compRow=refRow+1:length(ifp)% loop over cells that are below the cell of comparison
                    if contains(ifp{compRow,1},ifp{refRow,1})% if the parameter name in row #refRow is contained in that of row #compRow
                        dependentRow(end+1) = compRow;% store row #compRow into dependentRow list
                    end
                end
            end
            obj.indepFreeParams(dependentRow,:) = [];% update array of independent free parameters
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
                totalNumFreeParams = totalNumFreeParams + size(obj.freeParams{j1},1);% change to vertcat(obj.freeParams{:})
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
% loop through each possible additional free parameter, in the order of
% keys of referenceMap
                    infname = inFuncName(referenceMap,j2,paramArray);
                    fitstr = strcat(fitstr,infname,",");
% depending whether a parameter should be free or constrained in the fit,
% add its name or its value, respectively (see function 'inFuncName' above)
                end
                infname = inFuncName(referenceMap,length(referenceMap),paramArray);
                fitstr = strcat(fitstr,string(infname),"])");
                % end of string defining the fit equation
            end
            
            fitEqn = cell(1,nPeaks);
            for peakIdx=1:nPeaks
                fitEqn{peakIdx} = fitEqStr(obj.allParams{peakIdx},obj.freeParams{peakIdx});
            end
%             disp(fitEqn)
            fullFitStr = join(string(fitEqn),'+');% sum equation strings for all peaks
%             disp(fullFitStr);% display the fit equation for debugging
            
%% Create arrays containing lower bounds and initial values of fit parameters
            ifp = obj.indepFreeParams;% shorter variable name, easier to handle
            lowBounds = zeros(size(ifp(:,1)));%
            initParams = ones(size(ifp(:,1)));%
            % initialize arrays containing all the lower bounds and initial
            % values of free fit parameters, so that their length is sum(Ni)
            cntr = 0;% counter that increments at each iteration of the following loop
            ksj = keys(obj.allParams{1});
            for jj = 1:length(obj.allParams{1})-1
%                 ifp = 0;
%                 for j3=1:nPeaks
%                     lj3 = length(obj.freeParams{j3});
                if any(contains(ifp(:,1),ksj{jj}))
                    lgc = contains(ifp(:,1),ksj{jj});
                    valArr = ifp(:,2);
                    initParams(cntr+1:cntr+sum(lgc)) = [valArr{lgc}];
                    cntr = cntr+sum(lgc);
                end
%                 end
            end
%             for j3=1:nPeaks
            if any(contains(ifp(:,1),'x0'))
                lgc = contains(ifp(:,1),'x0');
                valArr = ifp(:,2);
                lowBounds(cntr+1:cntr+sum(lgc)) = -Inf;% x0 can be negative 
                initParams(cntr+1:cntr+sum(lgc)) = [valArr{lgc}];
            end
%             end
%             disp(lowBounds);% sprintf('%d\n',initParams)%display for debugging

%% Perform fit
            [xData, yData] = prepareCurveData(obj.X,obj.Y);
%             [xData, yData, wData] = prepareCurveData(obj.X,obj.Y,abs(1./obj.dY));%length(xData)
            % Set up fittype and options.
            ft = fittype( fullFitStr,'independent', 'x', 'dependent', 'y' );
            excludedPoints = excludedata( xData, yData, 'Indices', obj.dataExcl );
% obj.dataExcl should be defined on the x interval where obj.dY>0, otherwise
% the number of excluded points will be higher than the number of data points
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = lowBounds;% use the above defined arrays for lower bounds
            opts.StartPoint = initParams;% and starting parameter values
% the choice of initial parameters is critical to the convergence of the fit
            opts.Exclude = excludedPoints;
%             opts.Weights = wData;% there does not seem to be any good way
% of weighing the data points: 1/dY reduces significantly the importance of
% data points of the peaks, whereas Y/dY reduces the impact of data points
% with low intensity on the sides of the peak. However, both types of data points are important
            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts);
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

% Extract fit parameter values in order to plot fit curves for each peak individually
            nPks = length(obj.freeParams);
            fitPrms = ones(1,length(obj.allParams{1}));
            % no need to re-initialize the array at each iteration since
            % all sub-maps of obj.Params have the same length
            ks2 = keys(obj.allParams{1});% same
            funCell = cell(1,nPks);% initialize cell array containing both 
            % free and fixed parameters of the fit function to plot
            for ii=1:nPks
                for jj=1:length(ks2)
                    if find(contains(obj.freeParams{ii}(:,1),ks2{jj}))
% if the namestring of a parameter is contained in the array of free parameters
                        idxfp = contains(obj.freeParams{ii}(:,1),ks2{jj});
                        % index of parameter namestring in the first column of obj.freeParams{ii}
                        if any(contains(obj.indepFreeParams(:,1),obj.freeParams{ii}{idxfp}))
                            fitPrms(jj) = fitresult.(obj.freeParams{ii}{idxfp});
% if the current free parameter is an independent one, extract its value
% from the fit
                        else
% otherwise there are some additional characters in the free parameter
% name character array, which need be separated from the mere parameter name in
% order to get its fit value
                            idxifp = cell2mat(arrayfun(@(c) contains(...
                                obj.freeParams{ii}{idxfp},c),obj.indepFreeParams(:,1),'Uniform',0));
                            % cycle through array of independent free parameters to find the one that
                            % is contained in the current dependent parameter and store its array index
                            [sdiff,smatch] = strsplit(obj.freeParams{ii}{idxfp},obj.indepFreeParams{idxifp});
                            if ~isempty(sdiff{1})
                                fitPrms(jj) = str2double([sdiff{1} num2str(fitresult.(smatch{1}))]);
                            % Note: str2double does not work here as it 
                            % does not handle mathematical operations,
                            % whereas str2num does
                            elseif ~isempty(sdiff{2})
                                fitPrms(jj) = str2double([num2str(fitresult.(smatch{1})) sdiff{2}]);
                            else; fitPrms(jj) = fitresult.(smatch{1});
                            end
                        end
% use the parameter value resulting from the fit 
                    else; fitPrms(jj) = obj.allParams{ii}(ks2{jj});
% otherwise use the fixed value stored in the allParams property
                    end
                end
                funCell{ii} = @(x) fitPrms(1).*...
                    voigtIkedaCarpenter_ord(x,[fitPrms(2),fitPrms(3),...
                    fitPrms(4),fitPrms(5),fitPrms(6),fitPrms(7),fitPrms(8)]);
% define the function using the above extracted parameter values
            end
            minIndex = find(~obj.dataExcl,1,'first');% index of non-excluded datapoint with lowest x value 
            maxIndex = find(~obj.dataExcl,1,'last');% index of non-excluded datapoint with highest x value 
            Xfit = linspace(obj.X(minIndex),obj.X(maxIndex),1000);
            Yfit = fitresult(Xfit);% compute fit over a controlled number of points
% this allows to have fit plot with better resolution than the calculated one
            
% Plot fit with data.
            figure; hold on;
            if nPks>1% if there is more than one peak
                psub = cell(1,nPks);%
                for kk=1:nPks% plot fit curve for each peak
                    psub{kk} = fplot(funCell{kk},...
                        [obj.X(minIndex) obj.X(maxIndex)],'LineWidth',2);
                end
            end
            pfit = plot(Xfit,Yfit,'r-');
            pdat = errorbar(xData,yData,obj.dY,'.b','MarkerSize',18,'LineWidth',2);
            pexcl = plot(obj.X(excludedPoints),obj.Y(excludedPoints),'xk',...
                'MarkerSize',9);
            legend([pdat,pfit],'I vs. hh0','fit ICpV','Location','northeast');
%             legend([pdat,pexcl,pfit],'I vs. hh0','Excluded','fit ICpV');
            % Label axes
            xlabel("hh0"); ylabel("I (a.u.)");
            xlim([obj.X(minIndex)-0.05 obj.X(maxIndex)+0.05]);
            grid on
        end
        
    end
    
end