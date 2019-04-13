%% Create new workspace for new data analysis
clear all;% clear workspace, since variables change for each sample
M = 307.85;% molar mass of TmAsO4, in g/mol
fileImp = '2019-04-05_TmAsO4-LS5341-CHESS1-1'; m = 0.84*1e-3;% mass of sample measured, in g
% fileImp = '2019-04-09_TmAsO4-LS5341-CHESS1-2'; m = 2.47*1e-3;% mass of sample measured, in g
% fileImp = '2019-04-10_TmAsO4-LS5341-CHESS1-3'; m = 2.18*1e-3;% mass of sample measured, in g
DATA=ImportTmVO4Cp([fileImp '.dat']);% Use this data to plot Cp vs H, T superimposed on theoretical curves
% Average data points at each temperature + field 
%%
% H=[DATA.FieldOersted; DATA(2).FieldOersted];
% T=[DATA.SampleTempKelvin; DATA(2).SampleTempKelvin];
% Cp=[DATA.SampHCJmoleK; DATA(2).SampHCJmoleK];
H=[DATA.FieldOersted];
T=[DATA.SampleTempKelvin];
Cp=[DATA.SampHCJmoleK];
CpErr=[DATA.SampHCErrJmoleK];

whichPoints = isfinite(H) & isfinite(T) & isfinite(Cp);
H=H(whichPoints);
T=T(whichPoints);
Cp=Cp(whichPoints);
CpErr=CpErr(whichPoints);

fields = unique(round(H,-1));%[10,linspace(1000,4000,4),4500,4750,5000];
hmax=max(fields);

%% Plot 3D scatter (only for measurements at various magnetic fields)
figure
scatter3(H,T,Cp,'.')
% xlim([0 hmax])
ylim([0 3])
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('Cp (uJ/K)')
view(90,0)

%% (only for measurements at various magnetic fields)

steps = -50:50;
x= tstep*steps;
s=0.05;
d1Gaussian = -exp(-x.^2/(2*s^2)).*x./sqrt(s^6*2*pi);
d2Gaussian = exp(-x.^2/(2*s^2)).*(x.^2 - s^2)/sqrt(s^10*2*pi);

figure
d1Cpg = conv2(Cpg,d1Gaussian','same');
d2Cpg = conv2(Cpg,d2Gaussian','same');
surf(Hg,Tg,-d1Cpg,'EdgeColor','none');
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('-dCp/dT (uJ/K^2)')

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
Tsep = 6e-3;% Data points taken within an interval of Tsep are considered to be measured at the same temperature setpoint
% 6mK is an empirical estimate of the dispersion of data points taken consecutively at a given temperature. 
for i = 1:length(fields)
    averageCpData(i) = averageCp(Tsep,separatedCpData(i).T,...
        separatedCpData(i).Cp,separatedCpData(i).CpErr);
    averageCpData(i).H = separatedCpData(i).H(1:length(averageCpData(i).T));
end

%% Compute molar quantities, if accessible
R = 8.314;% gas constant
for i = 1:length(fields)
    averageCpData(i).Cpm = averageCpData(i).Cp*1e-6*M/m;% Cp measured in PPMS He4 is in uJ/K
    averageCpData(i).CpmErr = averageCpData(i).CpFullErr*1e-6*M/m;% Cp measured in PPMS He4 is in uJ/K
    averageCpData(i).Cpr = averageCpData(i).Cpm/R;% Cp measured in PPMS He4 is in uJ/K
    averageCpData(i).CprErr = averageCpData(i).CpmErr/R;% Cp measured in PPMS He4 is in uJ/K
end

%% Separate file name into date and sample name
fnamesep = strsplit(fileImp,'_');
date = fnamesep{1};
sample = fnamesep{2};

%% Plot averaged data at each field separately
figure
hold on
for i=1:length(fields)
    if isfield(averageCpData(i),'Cpr')
    errorbar(averageCpData(i).T,averageCpData(i).Cpr,averageCpData(i).CprErr,'.','MarkerSize',18)
    ylabel('$C_p/R$')
    else
    errorbar(averageCpData(i).T,averageCpData(i).Cp,averageCpData(i).CpFullErr,'.','MarkerSize',18)
    ylabel('C$_p$ (uJ/K)')
    end
end
title(sample)
xlabel('Temperature (K)');
legendCell = cellstr(num2str(fields, '%-d Oe'));
legend(legendCell)
hold off

%% Concatenate data for exportation 
if isfield(averageCpData(i),'Cpr')
exportArray = [averageCpData(i).H,averageCpData(i).T,...
averageCpData(i).Tsd,averageCpData(i).Cpr,averageCpData(i).CprErr];
else
exportArray = [averageCpData(i).H,averageCpData(i).T,...
averageCpData(i).Tsd,averageCpData(i).Cp,averageCpData(i).CpFullErr];
end

%% Export data to txt file
fileExp = [fileImp '_avg.txt'];
fileID = fopen(fileExp,'a');
fprintf(fileID,['Magnetic Field\tTemperature\tTemperature_StDev\tCp/R\tCp_Err/R\n'...
    'Oe\tK\tK\tJ/K/mol\tJ/K/mol\n']);
dlmwrite(fileExp,exportArray,'-append','delimiter','\t');
fprintf(fileID,'\n');
fclose(fileID);







  