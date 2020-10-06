%% Sample properties
m = 2.0e-3;% mass of sample, in g
M = 283.87;% molar mass of TmVO4, in g/mol
% Tc0rf = 2.126;% Value of transition temperature in this sample, in Kelvin units
% Tcnum = 2.15;% in Kelvin units
Hc = 5000;% in Oersted units; see data taken on needles of TmVO4-LS5200 in July 2017
% e = 1.5e-3;% constant longitudinal field

%% Data importation
sample = 'TmVO4-LS5228-HC20-X-1';
filename = '2020-10-02_TmVO4-LS5228-HC20-X-1.dat';
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2020-10_TmVO4-LS5228-HC20-X-1'
DataTm20X1=ImportTmVO4Cp(filename);% Use this data to plot color map of phase diagram

%% Assign data to variables
H=[DataTm20X1.FieldOersted];%*rescaling;
T=[DataTm20X1.SampleTempKelvin];
Cp=[DataTm20X1.SampHCJmoleK];
CpErr=[DataTm20X1.SampHCErrJmoleK];

whichPoints = isfinite(H) & isfinite(T) & isfinite(Cp);
H=H(whichPoints);
T=T(whichPoints);
Cp=Cp(whichPoints).*M/m*1e-6;% factor 1e-6 converts from uJ/K to J/K
CpErr=CpErr(whichPoints).*M/m*1e-6;%
[uhTm20X1,~,X] = unique(round(H,-1));

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
for i = 1:length(uhTm20X1)
    wp = abs(H-uhTm20X1(i))<50;
    separatedTm20X1CpData(i).H = H(wp);
    separatedTm20X1CpData(i).T = T(wp);
    separatedTm20X1CpData(i).Cp = Cp(wp);
    separatedTm20X1CpData(i).d1Cp = d1Cp(wp);
    separatedTm20X1CpData(i).CpErr = CpErr(wp);
    
    [separatedTm20X1CpData(i).T,wo] = sort(separatedTm20X1CpData(i).T);
    separatedTm20X1CpData(i).H = separatedTm20X1CpData(i).H(wo);
    separatedTm20X1CpData(i).Cp = separatedTm20X1CpData(i).Cp(wo);
    separatedTm20X1CpData(i).d1Cp = separatedTm20X1CpData(i).d1Cp(wo);
    separatedTm20X1CpData(i).CpErr = separatedTm20X1CpData(i).CpErr(wo);
end

%% Average datapoints taken at any given temperature and field setpoint
% We want to compute the average of data points that are repetitions of the
% same measurement, i.e. with same temperature and field setpoints
R = 8.314;
clear avgCpDy2Data
Tsep = 6e-3;% Data points taken within an interval of Tsep are considered to be measured at the same temperature setpoint
% 6mK is an empirical estimate of the dispersion of data points taken consecutively at a given temperature. 
% Tp = 27;% Temperature scale of phonons contribution to Cp in TmVO4, in K; see 'TmVO4_Cp_phonons.m'
for i=length(uhTm20X1):-1:1
avgCpTm20X1Data(i) = averageCpwithH2(Tsep,separatedTm20X1CpData(i).T,separatedTm20X1CpData(i).Cp,...
    separatedTm20X1CpData(i).CpErr,separatedTm20X1CpData(i).H);
end
% for i=length(uhTm20X1):-1:1
% avgCpTm20X1Data(i).Cpel = avgCpTm20X1Data(i).Cp - R*(avgCpTm20X1Data(i).T/Tp).^3;% electronic contribution to Cp, after subtracting phonons contribution
% avgCpTm20X1Data(i).Cpelr = avgCpTm20X1Data(i).Cpel/R;
% avgCpTm20X1Data(i).CpelrErr = avgCpTm20X1Data(i).CpFullErr/R;
% avgCpTm20X1Data(i).uh = unique(round(avgCpTm20X1Data(i).H,-2));
% end

%% Plot averaged data 
figure; hold on
for i = 1:length(uhTm20X1)
plot(avgCpTm20X1Data(i).T,avgCpTm20X1Data(i).Cp,'.-', 'MarkerSize', 12, 'LineWidth',1,...
    'DisplayName',sprintf('%i kOe',uhTm20X1(i)/10^3))
% plot(avgHC1Data(i).T,R*(avgHC1Data(i).T/Tp).^3,'.','DisplayName','$C_p^{\mathrm{phonons}}$')
% plot(avgHC1Data(i).T,avgHC1Data(i).Cpel,'.','DisplayName','$C_p^{\mathrm{full}}-C_p^{\mathrm{phonons}}$')
end
legend('show')
title('DyVO$_4$-LS5639-HC2 heat capacity')
xlabel('T (K)'); ylabel('$C_p$ (J/K/mol)');
xlim([0 30])
% zoom

%% Prepare fit of Cp vs Temperature 
i = 6;
% Use these variables in Curve fitting tool
clear fitT fitCp fitCpErr fitwghts 
fitT = avgCpTm20X1Data(i).T;
fitCp = avgCpTm20X1Data(i).Cp;
fitCpErr = avgCpTm20X1Data(i).CpFullErr/R;
fitwghts = 1./fitCpErr;
Tc = 2.2;

