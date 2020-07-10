%%     Analyze heat capacity from DR
% This routine is intended at analyzing Cp data acquired on 2020-07-01 on 
% 10Y-TmVO4-LS5435-IF030-A-FE-tip

%% To do: 
% * Error bars:
% *     For dT data, use standard deviation
% *     For dCp values, add the SampHCErrmJmoleK field 

%% Change directory 
cd 'C:\Users\Pierre\Desktop\Postdoc\YTmVO4\YTmVO4_HeatCapacity\Y10-TmVO4_HC\2020-07_10Y-TmVO4-LS5435-IF030-A-FE-tip'

%% Import YTmVO4 data 
Data = ImportCpFisherPPMSTable('2020-07-01_10Y-TmVO4-LS5435-IF030-A-FE-tip.dat');

%% Rename variables
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
% function to test if a column exists in a table
if isTableCol(Data,'SampleTempKelvin')
    Data.Properties.VariableNames{'SampleTempKelvin'} = 'T';% rename the temperature column
end
if isTableCol(Data,'FieldOersted')
    Data.Properties.VariableNames{'FieldOersted'} = 'H';% rename the magnetic field column
end
if isTableCol(Data,'SampHCuJK')% for samples measured in DR
    Data.Properties.VariableNames{'SampHCuJK'} = 'Cp';% rename the heat capacity column
elseif isTableCol(Data,'SampHCuJmoleK')% for samples measured in shared PPMS
    Data.Properties.VariableNames{'SampHCuJmoleK'} = 'Cp';
elseif isTableCol(Data,'SampHCJmoleK')% for samples measured in Fisher He4 PPMS
    Data.Properties.VariableNames{'SampHCJmoleK'} = 'Cp';
end
if isTableCol(Data,'SampHCErrJmoleK')% for samples measured in DR
    Data.Properties.VariableNames{'SampHCErrJmoleK'} = 'CpErr';
elseif isTableCol(Data,'SampHCErruJK')% for samples measured in DR
    Data.Properties.VariableNames{'SampHCErruJK'} = 'CpErr';        
end

%% Remove NaN rows
Data(any(isnan(Data.T), 2), :) = [];% Remove rows where T is NaN

%% Parameters for plotting heat capacity
xlblTemp = '$T$ (K)';
ylblCp = 'C$_p$ (J$\cdot$mol$^{-1}\cdot$K$^{-1}$)';
ttlCp = 'Heat capacity of TmVO$_4$';

%% Compute molar heat capacity
M = 275.87;% Molar mass of 10%Y-TmVO4 sample, in g/mol
m = 1e-3*18.8;% mass in g of the sample on which each dataset was measured
Data.Cpmol = Data.Cp *1e-6*M/m;% molar heat capacity, in J/mol/K
Data.CpmolErr = Data.CpErr *1e-6*M/m;% molar heat capacity, in J/mol/K
% starting from a heat capacity measured in microJoules per Kelvin, as is measured in the DR
% Cpmol is calculated per mole of Tm3+ ions, hence the (1-dpg) factor in the denominator

%% Keep only data under zero magnetic field
uh = unique(round(Data.H,-1));
for i=1:length(uh)%for all datasets
    split(i).unsorted = Data(round(Data.H,-1)==uh(i),:);% split data based on field value
    split(i).sorted = sortrows(split(i).unsorted ,{'T'});
end

%% Compute average of data points taken
clear avgData
R = 8.314;% gas constant in J/mol/K
Tp = 28.61;%  (28.54, 28.68);% Temperature scale of phonons contribution to Cp in TmVO4, 
% in K, obtained from fitting Cp data with R*(T/Tp)^3 between 11.12K and 21.45K
for i=1:length(uh)
    avgData(i) = averageCpwithH2(3e-3,split(i).sorted.T,split(i).sorted.Cpmol,...
        split(i).sorted.CpmolErr, split(i).sorted.H);
end% Not sure why I need to split the two loops here, but Matlab throws an error if I don't
for i=1:length(uh)
    avgData(i).CpFull = avgData(i).Cp;
    avgData(i).Cp = avgData(i).CpFull - R*(avgData(i).T/Tp).^3;% electronic contribution to Cp, after subtracting phonons contribution
    avgData(i).Cpr = avgData(i).Cp/R;
    avgData(i).CprErr = avgData(i).CpFullErr/R;
end

%% Plot averaged data 
i = 1;
figure
errorbar(avgData(i).T,avgData(i).CpFull,avgData(i).CpFullErr,'.','MarkerSize',18,...
    'DisplayName','$C_p^{\mathrm{full}}$')
hold on
plot(avgData(i).T,R*(avgData(i).T/Tp).^3,'.','DisplayName','$C_p^{\mathrm{phonons}}$')
errorbar(avgData(i).T,avgData(i).Cp,avgData(i).CpFullErr,'.','MarkerSize',18,...
    'DisplayName','$C_p^{\mathrm{full}}-C_p^{\mathrm{phonons}}$')
% plot(avgData(i).T,avgData(i).Cp,'.','DisplayName','$C_p^{\mathrm{full}}-C_p^{\mathrm{phonons}}$')
legend('show')
title('TmVO$_4$ heat capacity')
xlabel(xlblTemp); ylabel('$C_p$ (J/K/mol)');

%% Prepare fit of Cp vs Temperature 
% Use these variables in Curve fitting tool
fitT = avgData(i).T;
fitCp = avgData(i).CpFull;
fitCpErr = avgData(i).CpFullErr;
fitwghts = 1./fitCpErr;

%% Fit: 'Cp phonons contribution'.
[xData, yData, weights] = prepareCurveData( avgData(i).T, avgData(i).CpFull, 1./avgData(i).CpFullErr);

% Set up fittype and options.
ft = fittype( 'Schottky(A,D,T)+8.314*(T/Tp)^3', 'independent', 'T', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', [1 35:51] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Lower = [0 0 0];
opts.StartPoint = [1 1 10];
opts.Weights = weights;
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
prm = coeffvalues(fitresult);

% Plot fit with data.
Tfit = linspace(0,30,301);
fitsel = avgData(i).T>1.85 & avgData(i).T<27;
figure( 'Name', 'Cp phonons contribution' ); hold on;
pexcl = plot( avgData(i).T(~fitsel), avgData(i).CpFull(~fitsel), 'xk',...
    'MarkerSize', 6, 'DisplayName', 'Excluded');
h = errorbar( avgData(i).T(fitsel), avgData(i).CpFull(fitsel),...
    avgData(i).CpFullErr(fitsel), '.', 'MarkerSize', 18, 'DisplayName', 'Fitted data');
pfit = plot(Tfit, Schottky(prm(1),prm(2),Tfit)+8.314*(Tfit/prm(3)).^3,'r',...
    'DisplayName', 'Fit Schottky$+T^3$');
lgd = legend( 'Location', 'NorthEast' );
xlabel(xlblTemp); ylabel(ylblCp);
title('Tm$_{0.9}$Y$_{0.1}$VO$_4$ heat capacity')
xlim([0 30]);

%% Print figure to PNG file
% cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_summary_figures\TmVO4_GS-splitting+multipoles+Cp'
% formatFigure;
printPNG([todaystr '_10Y-TmVO4-LS5435-IF030-A-FE-tip_Cp+fit-Schottky+T3'])







