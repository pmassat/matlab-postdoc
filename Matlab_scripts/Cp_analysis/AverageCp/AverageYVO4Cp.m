%% Create new workspace for new data analysis
% clear all;% clear workspace, since variables change for each sample
cd 'C:\Users\Pierre\Desktop\Postdoc\YTmVO4\YTmVO4_HeatCapacity\YVO4\2019-05-10_YVO4-LS5259-HC1905-1\';
% M = 283.87;% molar mass of TmVO4, in g/mol
M = 203.85;% molar mass of YVO4, in g/mol
m = 1.79*1e-3;% mass of sample measured, in g
fileImp = '2019-05-10_YVO4-LS5259-HC1905-1'; 
DATA=ImportTmVO4Cp([fileImp '.dat']);% In case of error, first check that the cwd is correct

%%
% H=[DATA.FieldOersted; DATA(2).FieldOersted];
% T=[DATA.SampleTempKelvin; DATA(2).SampleTempKelvin];
% Cp=[DATA.SampHCJmoleK; DATA(2).SampHCJmoleK];
H=[DATA.FieldOersted];
T=[DATA.SampleTempKelvin];
% Cp=[DATA.SampHCJmoleK];
% CpErr=[DATA.SampHCErrJmoleK];
Cp=[DATA.AddendaHCJK];
CpErr=[DATA.AddendaHCErrJK];

whichPoints = isfinite(H) & isfinite(T) & isfinite(Cp);
H=H(whichPoints);
T=T(whichPoints);
Cp=Cp(whichPoints);
CpErr=CpErr(whichPoints);

fields = unique(round(H,-1));%[10,linspace(1000,4000,4),4500,4750,5000];
hmax=max(fields);

%%
figure
plot(T,Cp,'.')

%% Plot 3D scatter (only for measurements at different magnetic fields)
% figure
% plot(T,Cp,'.')
% % scatter3(H,T,Cp,'.')
% % xlim([0 hmax])
% % ylim([0 3])
% xlabel('Field (Oe)')
% ylabel('Temperature (K)')
% zlabel('C$_p$ ($\mu$J/K)')
% view(90,0)

%% (only for measurements at different magnetic fields)
% 
% steps = -50:50;
% x= tstep*steps;
% s=0.05;
% d1Gaussian = -exp(-x.^2/(2*s^2)).*x./sqrt(s^6*2*pi);
% d2Gaussian = exp(-x.^2/(2*s^2)).*(x.^2 - s^2)/sqrt(s^10*2*pi);
% 
% figure
% d1Cpg = conv2(Cpg,d1Gaussian','same');
% d2Cpg = conv2(Cpg,d2Gaussian','same');
% surf(Hg,Tg,-d1Cpg,'EdgeColor','none');
% xlabel('Field (Oe)')
% ylabel('Temperature (K)')
% zlabel('-dCp/dT (uJ/K^2)')

%% Create structure containing data separated according to magnetic field values
clear separatedCpData

for i = 1:length(fields)
    wp = abs(H-fields(i))<50;
    separatedCpData(i).H = H(wp);
    separatedCpData(i).T = T(wp);
    separatedCpData(i).Cp = Cp(wp);
    separatedCpData(i).CpErr = CpErr(wp);
    
    [separatedCpData(i).T,wo] = sort(separatedCpData(i).T);
    separatedCpData(i).H = separatedCpData(i).H(wo);
    separatedCpData(i).Cp = separatedCpData(i).Cp(wo);
    separatedCpData(i).CpErr = separatedCpData(i).CpErr(wo);
end

%% Average datapoints taken at any given temperature and field setpoint
% We want to compute the average of data points that are repetitions of the
% same measurement, i.e. with same temperature and field setpoints
relTsep = 3e-3;% Data points taken within a relative interval of T*(1+-Tsep) are considered to be measured at the same temperature setpoint
% 3e-3 is an empirical estimate of the relative dispersion of data points taken consecutively at a given temperature. 
for i = 1:length(fields)
    averageCpData(i) = averageCp_relT(relTsep,separatedCpData(i).T,...
        separatedCpData(i).Cp,separatedCpData(i).CpErr);
    averageCpData(i).H = separatedCpData(i).H(1:length(averageCpData(i).T));
end

%% Compute molar quantities, if accessible
for i = 1:length(fields)
    averageCpData(i).Cpm = averageCpData(i).Cp*1e-6*M/m;% Cp measured in PPMS He4 is in uJ/K
    averageCpData(i).CpmErr = averageCpData(i).CpFullErr*1e-6*M/m;% Cp measured in PPMS He4 is in uJ/K
end

%% Separate file name into date and sample name
fnamesep = strsplit(fileImp,'_');
date = fnamesep{1};
sample = fnamesep{2};

%% Fit molar heat capacity data
i = 1;
[xData, yData] = prepareCurveData( averageCpData(i).T, averageCpData(i).Cpm/R );

% Set up fittype and options.
ft = fittype( '(T/T1)^p', 'independent', 'T', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Weights', 1./averageCpData(i).CpmErr );
opts.Display = 'Off';
opts.StartPoint = [10 3];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
fitfunction = @(T) (T/fitresult.T1)^fitresult.p;

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = fplot( fitfunction, [0 70], '-r','LineWidth',2);
legend( h, 'avgCp vs. avgT', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'avgT', 'Interpreter', 'none' );
ylabel( 'avgCp', 'Interpreter', 'none' );
grid on

%% Plot averaged data at each field separately
% R = 8.314;% gas constant, in J/K/mol
figure
hold on
for i=1:length(fields)
    if isfield(averageCpData(i),'Cpm')
    errorbar(averageCpData(i).T,averageCpData(i).Cpm/R,averageCpData(i).CpmErr/R,'.','MarkerSize',18)
    ylabel('$C_p/R$')
    else
    errorbar(averageCpData(i).T,averageCpData(i).Cp,averageCpData(i).CpFullErr,'.','MarkerSize',18)
    ylabel('C$_p$ (uJ/K)')
    end
end
% fp = fplot( fitfunction, [0 75], '-r','LineWidth',1.5);
% title(sample)
xlabel('$T$ (K)');
% legend(fp, sprintf('$(T/%0.2g)^{%0.2g}$',fitresult.T1, fitresult.p), 'Location', 'northwest')
% legendCell = cellstr(num2str(fields, '%-d Oe'));
% legend(legendCell,'Location','best')
grid on
hold off

%% Print figure to pdf
% cd '.\2019-05-10_YVO4-LS5259-HC1905-1_analysis';
formatFigure
% printPDF([todaystr '_' sample '_Cpm']);

%% Prepare curve fit tool
avgT = averageCpData.T;
avgTsqr = avgT.^2;
avgCpm = averageCpData.Cpm;
avgCpmErr = averageCpData.CpmErr;
avgCpmT = avgCpm./avgT;
avgCpmTErr = avgCpmErr./avgT;
avgWeights = 1./avgCpmTErr;

%% Fit Cp/T vs T^2
[xTsqr, yCpmT, weightsCpmT] = prepareCurveData( avgTsqr, avgCpmT, avgWeights );

% Set up fittype and options.
ft = fittype( 'poly1' );
excludedCpmT = excludedata( xTsqr, yCpmT, 'Domain', [0 25] );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Weights = weightsCpmT;
opts.Exclude = excludedCpmT;
opts.Lower = [-Inf 0];
opts.Upper = [Inf 0];

% Fit model to data.
fitCpmT = fit( xTsqr, yCpmT, ft, opts );% Adjusted R-square: 0.9695 in Curve Fitting Tool
fCpmT = @(Tsq) fitCpmT.p1*Tsq + fitCpmT.p2;

%% Plot averaged data at each field separately
% R = 8.314;% gas constant, in J/K/mol
figure
hold on
errorbar(avgTsqr(~excludedCpmT), avgCpmT(~excludedCpmT), avgCpmTErr(~excludedCpmT),'.b','MarkerSize',18)
errorbar(avgTsqr(excludedCpmT), avgCpmT(excludedCpmT), avgCpmTErr(excludedCpmT),'.k','MarkerSize',18)
ylabel('$C_p/T$ (J/K$^2$/mol)')
fp = fplot(fCpmT, [0 25], '-r','LineWidth',1.5);
xlabel('$T^2$ (K$^2$)');
xlim([0 25]);
% ylim([0 0.15]);
% legend(fp, sprintf('$%0.2gT^2+%0.2g$',fitCpmT.p1, fitCpmT.p2), 'Location', 'northwest')
grid on
hold off

%% Print figure to pdf
% cd '.\2019-05-10_YVO4-LS5259-HC1905-1_analysis';
formatFigure
% printPDF([todaystr '_' sample '_CpmT']);

%% Concatenate data for exportation 
if isfield(averageCpData(i),'Cpm')
exportArray = [averageCpData(i).H,averageCpData(i).T,...
averageCpData(i).Tsd,averageCpData(i).Cpm,averageCpData(i).CpmErr];
else
exportArray = [averageCpData(i).H,averageCpData(i).T,...
averageCpData(i).Tsd,averageCpData(i).Cp,averageCpData(i).CpFullErr];
end

%% Export data to txt file
fileExp = [fileImp '_avg.txt'];
fileID = fopen(fileExp,'a');
fprintf(fileID,['Magnetic Field\tTemperature\tTemperature_StDev\tCp\tCp_Err\n'...
    'Oe\tK\tK\tJ/K/mol\tJ/K/mol\n']);
dlmwrite(fileExp,exportArray,'-append','delimiter','\t');
fprintf(fileID,'\n');
fclose(fileID);







  