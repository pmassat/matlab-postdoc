DATA=ImportTmVO4Cp('TmVO4_Mosaic_2017-07-28.dat');% Use this data to plot Cp vs H, T superimposed on theoretical curves
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

hmax=5500;
scatter3(H,T,Cp,'.')
xlim([0 hmax])
ylim([0 3])
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('Cp (uJ/K)')

%%
tstep=0.05;
hstep=500;
[Hg,Tg] = meshgrid(0:hstep:hmax,0.365:tstep:3);% for griddata
Hgl=0:hstep:hmax;% for gridfit
Tgl=0.365:tstep:3;% for gridfit
figure
Cpg = gridfit(H,T,Cp,Hgl,Tgl);
%Cpg = griddata(H,T,Cp,Hg,Tg);
surf(Hg,Tg,Cpg)
hold on
scatter3(H,T,Cp,100,'.','.k')
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('Cp (uJ/K)')

%%

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


%%
clear separatedCpData

fields = unique(round(H,-1));%[10,linspace(1000,4000,4),4500,4750,5000];

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

%% Plot averaged data
figure; hold on
rng = [1:1:length(fields)];
for i=rng
    errorbar(separatedCpData(i).Tm,separatedCpData(i).Cpm,separatedCpData(i).CpmFullErr,'.','MarkerSize',18)
end
xlabel('Temperature (K)');
ylabel('C_p (J\cdotmol^{-1}\cdotK^{-1})');

%% Prepare MF fit of Cp vs Temperature under field
% Useful to use the curve fitting tool for quick look at the data
index = 1;
H1 = fields(index);
T1 = separatedCpData(index).Tm;
Cp1 = separatedCpData(index).Cpm;
Cp1Err = separatedCpData(index).CpmFullErr;
wghts1 = 1./Cp1Err;

%% Plot averaged data at each field separately
figure; hold on
rng = 1:length(fields)-1;
clr = cell(size(rng));
eb = cell(size(rng));
for i=rng
    fp = fplot(@(t)31.5*Cp_TFIM(t/2.21,fields(i)/(5.15e3)),[0 4],'LineWidth',2);
% 31.5 is the value of amplitude coefficient extracted from fit using curve
% fitting tool
    clr{rng==i} = get(fp,'Color');
end
for i=rng
    eb{rng==i} = errorbar(separatedCpData(i).Tm,separatedCpData(i).Cpm,separatedCpData(i).CpmFullErr,...
        'x','MarkerSize',18,'DisplayName',num2str(fields(i)/1e4,'%.2f T'),...
        'Color',clr{rng==i},'LineWidth',2);
end
xlabel('Temperature (K)'); ylabel('Cp (uJ/K)');
title('Heat capacity of needles of TmVO4 (no demag) at various fields')
legend([eb{:}]);
hold off


%% Gaussian convolution and computation of derivative of raw Cp
tstep=0.025; 
steps = -50:50;
x= tstep*steps;
s=0.05;
d1Gaussian = -exp(-x.^2/(2*s^2)).*x./sqrt(s^6*2*pi);
d2Gaussian = exp(-x.^2/(2*s^2)).*(x.^2 - s^2)/sqrt(s^10*2*pi);
d1Cp = conv2(Cp,d1Gaussian','same');

%% Compute and plot derivative of averaged data
figure; hold on
for i = 1:length(fields)
    separatedCpData(i).d1Cpm = conv2(separatedCpData(i).Cpm,d1Gaussian','same');
    plot(separatedCpData(i).Tm,separatedCpData(i).d1Cpm)
end

%% Identify experimental critical temperature at each field
% Then we can plot Tc vs fields and fit using equation f(h) = Tc/Hc*h/atanh(h/Hc)
% From this, we can see that the experimental value of Hc is ~0.72T 
% and hence correct for demag, knowing the value of critical field Hc0 (see beginning of code)
M = ones(1,length(fields));
Tc = ones(1,length(fields));
for i=1:length(fields)
    [M(i),I] = min(separatedCpData(i).d1Cpm);
    Tc(i) = separatedCpData(i).Tm(I);
end


%% Plot data at each field separately
% Probably useless now that the average is computed. 2019-03-18
figure
for i=1:length(fields)
    plot(separatedCpData(i).T,separatedCpData(i).Cp,'-+')
    hold on
end
xlabel('Temperature (K)')
ylabel('Cp (uJ/K)')
legendCell = cellstr(num2str(fields, '%-d Oe'));
legend(legendCell)
hold off

%%
for i=1:length(fields)
    separatedCpData(i).f=fitSpline(separatedCpData(i).T,-separatedCpData(i).Cp);
%Is fitSpline a Matlab function? Cannot find it in the help.
%    separatedCpData(i).d2=differentiate(separatedCpData(i).f,separatedCpData(i).T);
    plot(separatedCpData(i).f,'deriv1')
    hold on
end
legend(legendCell)
xlabel('Temperature (K)')
ylabel('-dCp/dT (uJ/K^2)')
hold off

%% Plot heat capacity derivative vs temperature and field
figure
% n = length(fields);
% C=[]
% for i=0:n-1
%     C=[C; 1-i/n i/n 0];
% end
% set(0,'defaultaxescolororder',C) %red to green
for i=1:length(fields)
    separatedCpData(i).f=fitSpline(separatedCpData(i).T,separatedCpData(i).Cp);
    separatedCpData(i).d2=differentiate(separatedCpData(i).f,separatedCpData(i).T);
    scatter3(separatedCpData(i).H,separatedCpData(i).T,-separatedCpData(i).d2,48,'filled')
    hold on
end
xlim([0 hmax])
ylim([0 3])
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('-dCp/dT (uJ/K)')
hold off

%% Fit and plot heat capacity derivative vs temperature and field
T2=[];
H2=[];
d2Cp2=[];
for i=1:length(fields)
    T2= [T2; separatedCpData(i).T];
    H2=[ H2; separatedCpData(i).H];
    d2Cp2=[d2Cp2; separatedCpData(i).d2];
end

tstep=0.025;
% hmax=8000;
% hstep=500;
% [Hg,Tg] = meshgrid(0:hstep:hmax,0.365:tstep:3);%for griddata
Hgl=0:hstep:hmax;%for gridfit
Tgl=0.365:tstep:3;%for gridfit
figure
d2Cpgf = gridfit(H2,T2,d2Cp2,Hgl,Tgl);
d2Cpg = griddata(H2,T2,d2Cp2,Hg,Tg);
surf(Hg,Tg,-d1Cpg,'EdgeColor','none')
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('-dCp/dT (uJ/K^3)')
hold on
%scatter3(H2,T2,-d2Cp2,'.k')
xlim([0 hmax])
ylim([0 3])

%% Identify experimental critical temperature at each field
M = ones(1,length(fields));
Tc = ones(1,length(fields));
for i=1:length(fields)
    [M(i),I] = min(separatedCpData(i).d2);
    Tc(i) = separatedCpData(i).T(I);
end

%% Fit: 'Phase boundary'.
[xData, yData] = prepareCurveData( fields, Tc );

% Set up fittype and options.
ft = fittype( 'Tc0/Hc0*x/atanh(x/Hc0)', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', [7 8] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [5000 2.2];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData, excludedPoints );
xlim([0 6000]); ylim([0 2.5])
legend( h, 'Tc vs. fields', 'Excluded Tc vs. fields', 'Phase boundary', 'Location', 'best' );
% Label axes
xlabel fields; ylabel Tc;
grid on
%% Print fit parameter values with error bars
cval = coeffvalues(fitresult);% extract fit parameter values
cft=confint(fitresult);% extract confidence intervals from fit
sprintf("Hc(T=0) = %.1d +- %.0e Oe",fitresult.Hc0,cval(1)-cft(1,1))% print out value of critical field, with error bars
sprintf("Tc(H=0) = %.2f +- %.2f K",fitresult.Tc0,cval(2)-cft(1,2))% print out value of critical temperature, with error bars








