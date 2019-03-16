% fileImp = '2019-02-27_TmVO4-LS5228-RUS-1';
% fileImp = '2019-02-27_TmVO4-LS5228-RUS-1-S';
fileImp = '2019-03-01_TmVO4-LS5249-RUS-2';
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

%% Plot 3D scatter (only for measurements at different magnetic fields)
figure
scatter3(H,T,Cp,'.')
% xlim([0 hmax])
ylim([0 3])
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('Cp (uJ/K)')
view(90,0)

%% (only for measurements at different magnetic fields)

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
    T2 = repmat(separatedCpData(i).T,1);
    Tm = zeros(length(separatedCpData(i).T),1);
    Tsd = zeros(length(separatedCpData(i).T),1);
    Cpm = zeros(length(separatedCpData(i).Cp),1);
    Cpsd = zeros(length(separatedCpData(i).Cp),1);
    CpmErr = zeros(length(separatedCpData(i).Cp),1);
    for k = 1:length(T2)
        if T2(k)==0
            continue
        elseif length(T2(abs(T2-T2(k))<Tsep))>3
            Tsep2 = Tsep/2;% reduce the temperature interval
%             T2(abs(T2-T2(k))<Tsep2)%print out values of temperature which
%             verify the if statement
            Tm(k) = mean(T2(abs(T2-T2(k))<Tsep2));
            Tsd(k) = std(T2(abs(T2-T2(k))<Tsep2));
            Cpm(k) = mean(separatedCpData(i).Cp(abs(T2-T2(k))<Tsep2));
            Cpsd(k) = std(separatedCpData(i).Cp(abs(T2-T2(k))<Tsep2));
            CpmErr(k) = sum(separatedCpData(i).CpErr(abs(T2-T2(k))<Tsep2))/...
                sqrt(length(separatedCpData(i).CpErr(abs(T2-T2(k))<Tsep2)));
            T2(abs(T2-T2(k))<Tsep2)=0;
        else
            Tm(k) = mean(T2(abs(T2-T2(k))<Tsep));
            Tsd(k) = std(T2(abs(T2-T2(k))<Tsep));
            Cpm(k) = mean(separatedCpData(i).Cp(abs(T2-T2(k))<Tsep));
            Cpsd(k) = std(separatedCpData(i).Cp(abs(T2-T2(k))<Tsep));
            CpmErr(k) = sum(separatedCpData(i).CpErr(abs(T2-T2(k))<Tsep))/...
                sqrt(length(separatedCpData(i).CpErr(abs(T2-T2(k))<Tsep)));
            T2(abs(T2-T2(k))<Tsep)=0;
        end
    end
    separatedCpData(i).Hm = separatedCpData(i).H(Tm>0);
    separatedCpData(i).Tm = Tm(Tm>0);
    separatedCpData(i).Tsd = Tsd(Tm>0);
    separatedCpData(i).Cpm = Cpm(Cpm>0);
    separatedCpData(i).CpmFullErr = Cpsd(Cpm>0) + CpmErr(Cpm>0);
end

%% Plot averaged data at each field separately
figure
hold on
for i=1:length(fields)
    errorbar(separatedCpData(i).Tm,separatedCpData(i).Cpm,separatedCpData(i).CpmFullErr,'.','MarkerSize',18)
end
xlabel('Temperature (K)');
ylabel('Cp (uJ/K)')
legendCell = cellstr(num2str(fields, '%-d Oe'));
legend(legendCell)
hold off

%% Concatenate data for exportation 
exportArray = [separatedCpData(i).Hm,separatedCpData(i).Tm,...
separatedCpData(i).Tsd,separatedCpData(i).Cpm,separatedCpData(i).CpmFullErr];

%% Export data to txt file
fileExp = [fileImp '_avg.txt'];% '_R=' sprintf('%2.e',ap1('R'))
fileID = fopen(fileExp,'a');
fprintf(fileID,['Magnetic Field\tTemperature\tTemperature_StDev\tCp\tCp_Err\n'...
    'Oe\tK\tK\tuJ/K\tuJ/K\n']);
% for i=1:length(fields)
% fprintf(fileID,'%d\n\n',exportArray);
% end
% export(exportArray,'File',fileID,'Delimiter',',')
dlmwrite(fileExp,exportArray,'-append','delimiter','\t');
fprintf(fileID,'\n');
fclose(fileID);







  