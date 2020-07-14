% Code was derived from 'YTmVO4_primary_analysis.m' on 2019-08-06

%% Analyze heat capacity from DR
% This routine is intended at anaLzing Cp data acquired with our DR used in 
% the Dynacool PPMS of the Lee lab
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2020-07_TmVO4-LS5228-DR-HC180731-2007'

%% To do: 
% * Error bars:
% *     For dT data, use standard deviation

%% Import data YTmVO4
Data = ImportCpDR('2020-07_TmVO4-LS5228-DR-HC180731-2007.dat');

%% Rename variables
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
% function to test if a column exists in a table
if isTableCol(Data,'SampleTempKelvin')
    Data.Properties.VariableNames{'SampleTempKelvin'} = 'T';% rename the temperature column
end
if isTableCol(Data,'FieldOersted')
    Data.Properties.VariableNames{'FieldOersted'} = 'H';% rename the magnetic field column
end
if isTableCol(Data,'SampHCJK')% for samples measured in DR
    Data.Properties.VariableNames{'SampHCJK'} = 'Cp';% rename the heat capacity column
elseif isTableCol(Data,'SampHCmJmoleK')% for samples measured in shared PPMS
    Data.Properties.VariableNames{'SampHCmJmoleK'} = 'Cp';
elseif isTableCol(Data,'SampHCJmoleK')% for samples measured in Fisher He4 PPMS
    Data.Properties.VariableNames{'SampHCJmoleK'} = 'Cp';
end
if isTableCol(Data,'SampHCErrJmoleK')% for samples measured in DR
    Data.Properties.VariableNames{'SampHCErrJmoleK'} = 'CpErr';
elseif isTableCol(Data,'SampHCErrJK')% for samples measured in DR
    Data.Properties.VariableNames{'SampHCErrJK'} = 'CpErr';        
end

%% Remove NaN rows
Data(any(isnan(Data.T), 2), :) = [];% Remove rows where T is NaN

%% Parameters for plotting heat capacity
xlblTemp = '$T$ (K)';
ylblCp = 'C$_p$ (J$\cdot$mol$^{-1}\cdot$K$^{-1}$)';
ttlCp = 'Heat capacity of TmVO$_4$';

%% Compute molar heat capacity
M = 283.87326;% Molar mass of each sample, in g/mol
m = 1e-3*0.38;% mass in g of the sample on which each dataset was measured
Data.Cpmol = Data.Cp *1e-6*M/m;% molar heat capacity, in J/mol/K
Data.CpmolErr = Data.CpErr *1e-6*M/m;% molar heat capacity, in J/mol/K
% starting from a heat capacity measured in microJoules per Kelvin, as is measured in the DR
% Cpmol is calculated per mole of Tm3+ ions, hence the (1-dpg) factor in the denominator

%% Split data by magnetic field
uh = unique(round(Data.H,-1));
for i=1:length(uh)%for all datasets
    split(i).unsorted = Data(round(Data.H,-1)==uh(i),:);% split data based on field value
    split(i).sorted = sortrows(split(i).unsorted ,{'T'});
end

%% Compute average of data points taken
clear avgData
R = 8.314;% gas constant in J/mol/K
Tp = 30;% Temperature scale of phonons contribution to Cp in TmVO4, in K; see 'TmVO4_Cp_phonons.m'
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
i = 2;
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

%% Prepare plot of theoretical heat capacity curve 
% Tc = 2.19;%  Value of transition temperature obtained from fit with Cp_LFIM_NF function in Curve Fitting Tool
Tc = 2.185;%  Value of transition temperature obtained from fit with Cp_LFIM_CW function in Curve Fitting Tool
s =  5e-4;%  Value of longitudinal field obtained from fit with Cp_LFIM_NF and Cp_LFIM_CW and Cp_LFIM_MEC_cst_stress function in Curve Fitting Tool
T_CW = 2.09;% Curie-Weiss temperature obtained from fit with function Cp_LFIM_CW 
A_CW = 0.0136;% amplitude of CW divergence obtained from fit with function Cp_LFIM_CW 
T0 = 1.6;% Divergence temperature obtained from fit with function Cp_LFIM_MEC_cst_stress
A = 5e-3;% Amplitude of "high"-temperature divergence obtained from fit with function Cp_LFIM_MEC_cst_stress
tvec = linspace(1e-3,4,300);% temperature vector for computation of Cp

%% Fit data at H = 10 kOe to determine the temperature offset of the He4 HC puck #666
i = 2;
[xData, yData, weights] = prepareCurveData( avgData(i).T, avgData(i).Cpr, 1./avgData(i).CprErr);
% Set up fittype and options.
ft = fittype( 'Cp_TFIM((T-Toffset)/2.185,0)', 'independent', 'T', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', 3:12 );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint = 0.1;
opts.Weights = weights;
opts.Exclude = excludedPoints;
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
% Plot fit with data.
figure( 'Name', 'Fit TFIM below $T_c$ at $H = 10 kOe$' );
h = plot( fitresult, xData, yData, excludedPoints );
legend( h, 'Fitted data', 'Excluded', 'TFIM fit', 'Location', 'NorthEast' );
xlabel(xlblTemp); ylabel('$C_p/R$');
coeff = coeffvalues(fitresult);
confbounds = confint(fitresult);
printfitres = annotation('textbox',[.4 .4 .4 .1],'FontSize',12,'LineWidth',1,...
    'String',sprintf('Toffset = %.3f (%.3f %.3f)', coeff, confbounds),...
    'FitBoxToText','on');
% Result
puck666lowToffset = 0.04;% Kelvin; result from fit of 10 kOe data: Toffset = 0.03818 with (-0.1037, 0.18) 95% confidence bounds

%% Fit data at H = 70 kOe 
% A Schottky anomaly doesn't fit!
i = 7;
[xData, yData, weights] = prepareCurveData( avgData(i).T, avgData(i).Cpr, 1./avgData(i).CprErr);
% Set up fittype and options.
% Gaussian fit including phonons contribution and the temperature offset of
% the puck calibration
ft1 = fittype( 'A*exp(-((T-b)/c)^2)+((T-0.04)/30)^3', 'independent', 'T', 'dependent', 'y' );
% Schottky fit including phonons contribution and the temperature offset...
ft2 = fittype( 'Schottky(A,D,T-0.04)+((T-0.04)/30)^3', 'independent', 'T', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Weights = weights;
% Fit model to data.
opts.StartPoint = [1 3 1];
[fitresult1, gof1] = fit( xData, yData, ft1, opts );
opts.StartPoint = [1 3];
[fitresult2, gof2] = fit( xData, yData, ft2, opts );
% Plot fit with data.
figure( 'Name', 'Fit data at 70 kOe' ); hold on
dp = plot( xData, yData,'.', 'DisplayName', 'Fitted data');
h1 = plot( fitresult1, 'g');
h2 = plot( fitresult2);
legend([dp, h1, h2], 'Data', 'Gaussian', 'Schottky');
xlabel(xlblTemp); ylabel('$C_p/R$');
title('TmVO$_4$ heat capacity, H=70 kOe // $a$')

%% Prepare fit of Cp vs Temperature 
% Use these variables in Curve fitting tool
i = 9
fitT = avgData(i).T;
fitCp = avgData(i).Cpr;
fitCpErr = avgData(i).CprErr;
fitwghts = 1./fitCpErr;



%% % % Plot & export figure % % % 

%% Figure formatting when combining plots
ax1.YLabel.Position(1) = ax2.YLabel.Position(1);
anngap.FontSize = ax1.XAxis.Label.FontSize;
ax1.XTickLabel={};
ax1.YTick = [];
ax2.XLabel.Position(2) = -.15;

%% Print figure to PNG file
% cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_summary_figures\TmVO4_GS-splitting+multipoles+Cp'
% formatFigure;
printPNG([todaystr '_TmVO4-LS5228-DR-HC180731_He4-2007_70kOe_Cp+fit'])

%% Export averaged Cp data to text file
% savepath = 'C:\Users\Pierre\Desktop\Postdoc\YTmVO4\YTmVO4_HeatCapacity\YTmVO4_Cp_anaLsis\2018-10-17_YTmVO4_averaged_Cp\';
% for i=L-3:L
%     dpg(i);
%     expstr(i) = fullfile(savepath, ['AnaLzeCpDRdoping_2018-11-30_YTmVO4_x=',num2str(dpg(i)),'.txt']);
%     fid = fopen(expstr(i), 'wt');
%     fprintf(fid, '%s\n', ['Exported from AnaLzeCpDRdoping_2018-10-10_YTmVO4_VTmAsO4.mlx on ',date]);  % header
%     fclose(fid);
%     dlmwrite(expstr(i),cat(2,avgData(i).T',avgData(i).Cp',avgData(i).stdCp'),'-append','delimiter','\t')
% end
% 

