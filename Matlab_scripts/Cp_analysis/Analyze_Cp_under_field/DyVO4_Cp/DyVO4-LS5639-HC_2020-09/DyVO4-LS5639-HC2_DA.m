%% Sample properties
m = 0.62e-3;% mass of sample, in g
M = 277.44;% molar mass of TmVO4, in g/mol
% Tc0rf = 2.126;% Value of transition temperature in this sample, in Kelvin units
% Tcnum = 2.15;% in Kelvin units
% Hc = 5000;% in Oersted units; see data taken on needles of TmVO4-LS5200 in July 2017
% e = 1.5e-3;% constant longitudinal field
% Results for Tmaxfit = 2.15;% Tc = 2.126  (2.122, 2.129);% e = 0.001474 (0.001085, 0.001862);

%% Data importation
sample = 'DyVO4-LS5639-HC2';
filename = '2020-09-29_DyVO4-LS5639-HC2.dat';
cd 'C:\Users\Pierre\Desktop\Postdoc\DyVO4\DyVO4_Cp\DyVO4-LS5639-HC_2020-09'
DataDy2=ImportTmVO4Cp(filename);% Use this data to plot color map of phase diagram

%% Assign data to variables
H=[DataDy2.FieldOersted];%*rescaling;
T=[DataDy2.SampleTempKelvin];
Cp=[DataDy2.SampHCJmoleK];
CpErr=[DataDy2.SampHCErrJmoleK];

whichPoints = isfinite(H) & isfinite(T) & isfinite(Cp);
H=H(whichPoints);
T=T(whichPoints);
Cp=Cp(whichPoints).*M/m*1e-6;% factor 1e-6 converts from uJ/K to J/K
CpErr=CpErr(whichPoints).*M/m*1e-6;%
[uhDy2,~,X] = unique(round(H,-1));

%% Gaussian convolution
rescaling = 1;
tstep=0.025; 
hstep=500*rescaling;
steps = -50:50;
x= tstep*steps;
s=0.05;
d1Gaussian = -exp(-x.^2/(2*s^2)).*x./sqrt(s^6*2*pi);
d2Gaussian = exp(-x.^2/(2*s^2)).*(x.^2 - s^2)/sqrt(s^10*2*pi);
d1Cp = conv2(Cp,d1Gaussian','same');

%% Separate data according to value of field
clear separatedDy2CpData
for i = 1:length(uhDy2)
    wp = abs(H-uhDy2(i))<50;
    separatedDy2CpData(i).H = H(wp);
    separatedDy2CpData(i).T = T(wp);
    separatedDy2CpData(i).Cp = Cp(wp);
    separatedDy2CpData(i).d1Cp = d1Cp(wp);
    separatedDy2CpData(i).CpErr = CpErr(wp);
    
    [separatedDy2CpData(i).T,wo] = sort(separatedDy2CpData(i).T);
    separatedDy2CpData(i).H = separatedDy2CpData(i).H(wo);
    separatedDy2CpData(i).Cp = separatedDy2CpData(i).Cp(wo);
    separatedDy2CpData(i).d1Cp = separatedDy2CpData(i).d1Cp(wo);
    separatedDy2CpData(i).CpErr = separatedDy2CpData(i).CpErr(wo);
end

%% Average datapoints taken at any given temperature and field setpoint
% We want to compute the average of data points that are repetitions of the
% same measurement, i.e. with same temperature and field setpoints
R = 8.314;
clear avgCpDy2Data
Tsep = 6e-3;% Data points taken within an interval of Tsep are considered to be measured at the same temperature setpoint
% 6mK is an empirical estimate of the dispersion of data points taken consecutively at a given temperature. 
Tp = 27;% Temperature scale of phonons contribution to Cp in TmVO4, in K; see 'TmVO4_Cp_phonons.m'
for i=length(uhDy2):-1:1
avgCpDy2Data(i) = averageCpwithH2(Tsep,separatedDy2CpData(i).T,separatedDy2CpData(i).Cp,...
    separatedDy2CpData(i).CpErr,separatedDy2CpData(i).H);
end
for i=length(uhDy2):-1:1
avgCpDy2Data(i).Cpel = avgCpDy2Data(i).Cp - R*(avgCpDy2Data(i).T/Tp).^3;% electronic contribution to Cp, after subtracting phonons contribution
avgCpDy2Data(i).Cpelr = avgCpDy2Data(i).Cpel/R;
avgCpDy2Data(i).CpelrErr = avgCpDy2Data(i).CpFullErr/R;
avgCpDy2Data(i).uh = unique(round(avgCpDy2Data(i).H,-2));
end

%% Plot averaged data 
figure; hold on
for i = 1:length(uhDy2)
plot(avgCpDy2Data(i).T,avgCpDy2Data(i).Cp,'.-', 'MarkerSize', 12, 'LineWidth',1,...
    'DisplayName',sprintf('%i kOe',uhDy2(i)/10^3))
% plot(avgHC1Data(i).T,R*(avgHC1Data(i).T/Tp).^3,'.','DisplayName','$C_p^{\mathrm{phonons}}$')
% plot(avgHC1Data(i).T,avgHC1Data(i).Cpel,'.','DisplayName','$C_p^{\mathrm{full}}-C_p^{\mathrm{phonons}}$')
end
legend('show')
title('DyVO$_4$-LS5639-HC2 heat capacity')
xlabel('T (K)'); ylabel('$C_p$ (J/K/mol)');
xlim([0 30])
% zoom

%% Prepare fit of Cp vs Temperature 
i = 9;
% Use these variables in Curve fitting tool
clear fitT fitCp fitCpErr fitwghts 
fitT = avgCpDy2Data(i).T;
fitCp = avgCpDy2Data(i).Cpelr;
fitCpErr = avgCpDy2Data(i).CpFullErr/R;
fitwghts = 1./fitCpErr;
% Tc = 2.15;

