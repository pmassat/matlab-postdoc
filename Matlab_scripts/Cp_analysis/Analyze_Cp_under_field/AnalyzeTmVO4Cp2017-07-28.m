Tc = 2.22;% transition temperature at zero field, in Kelvin units
Hc = 5100;% critical field at zero temperature, in Oersted units

%% Import magnetic field distribution from CSV file
cd 'C:\Users\Pierre\Desktop\Postdoc\Software\COMSOL\TmVO4-LS5200_HC2017-07\TmVO4-LS5200_HC2017-07_COMSOL_results'
% mfdfilename = cell(1,2);
% S = cell(1,2);
mfdfilename{1} = 'TmVO4-LS5200_HC2017-07_sample1-Needle_mesh-p1mm_T=p3-p4-3p1_Hext=all.csv';
mfdfilename{2} = 'TmVO4-LS5200_HC2017-07_sample2-Arya_mesh-20um_T=p3-p4-3p1_Hext=all.csv';
mfdDomNum = {2,3};
for ic=1:length(S)
    S{ic} = importfielddistrib_csv(mfdfilename{ic}, 56, 'domainNum', mfdDomNum{ic});
end

%% Extract values of temperature and external magnetic field from structure header
Smfd = [S{1,:}];% concatenate structures
for sidx=1:length(Smfd)
    % split each header into a cell array of strings, using whitespace and '=' sign as separators
    TBextCell = strsplit(Smfd(sidx).T_Bext,{' ','='});
    % Find the index of the string cell containing "T" (temperature), and
    % convert the following string cell (which contains the value of temperature) into a number
    % Note: circshift(A,1) shifts the index of elements of array A by 1 to the right
    Smfd(sidx).T_K = str2double(TBextCell{circshift(TBextCell=="T",1)});
    % Same with "Bext", which is the external magnetic flux density, in Tesla
    % units, and convert it into a value of magnetic field, in Oersted units
    Smfd(sidx).Hext_Oe = str2double(TBextCell{circshift(TBextCell=="Bext",1)})*10^4;
end

%% Compute probability distribution of fields at a given value of T and Hext
for i=1:length(Smfd)
mfd = Smfd(i).mfd(Smfd(i).mfd>0)/Smfd(i).Hext_Oe;% create distribution from non-zero values
% h = histogram(mfd, 'Normalization', 'pdf');% plot histogram of distribution
[Smfd(i).hc, edges] = histcounts(mfd, 50, 'Normalization', 'pdf');% plot histogram of distribution
% binCenters = h.BinEdges + (h.BinWidth/2);
Smfd(i).binCenters = mean([edges(1:end-1);edges(2:end)],1);
Smfd(i).binWidths = edges(2:end)-edges(1:end-1);
end

%% Plot distribution of fields at a given value of T and Hext
figure;
hold on
needle_index = 2;% there are 2 needle-shaped samples 
start_index = [0,56];% needle 1 data start at index 0+1, needle 2 data at 56+1
param_index = 2;% 1 is constant T, 2 is constant Hext, see param_range
param_range = {[4:8:56], [41:48]};% first range corresponds to a 
% field dependence at constant temp, second range corresponds to a 
% temperature dependence at constant field
rng = start_index(needle_index)  + param_range{param_index};
for i=rng 
T = Smfd(i).T_K;%
Hext = Smfd(i).Hext_Oe;%
p = plot(Smfd(i).binCenters, Smfd(i).hc, '.-', 'DisplayName', sprintf('%.2g, %.2g',T/Tc,Hext/Hc));
end
lgd = legend('show','Location','northwest'); lgd.Title.String = '$T/T_c$, $H_{\mathrm{ext}}/H_c$';
needle_title = sprintf('Distribution of fields in TmVO$_4$ needle %d', needle_index);
param_title = {[', $T=$' sprintf('%.2g K',T)],...
    [', $H_{\mathrm{ext}}=$' sprintf('%.2d Oe',Hext)]};
title([needle_title param_title{param_index}])
xlabel('$H_{\mathrm{in}}/H_{\mathrm{ext}}$')
ylabel('Normalized PDF')

%% Import HC data
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2017-07_TmVO4_Cp_MCE\2017-07-28--31\2017-07-28_Cp'
DATA=ImportTmVO4Cp('TmVO4_Mosaic_2017-07-28.dat');
% Use this data to plot Cp vs (H, T) superimposed on theoretical curves
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
ylabel('$C_p$ (J/K/mol)')

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
    fp = fplot(@(t)31.5*Cp_TFIM(t/2.21,fields(i)/(5.15e3)),[0 3],'LineWidth',2);
% 31.5 is the value of amplitude coefficient extracted from fit using curve
% fitting tool
    clr{rng==i} = get(fp,'Color');
end
for i=rng
    eb{rng==i} = errorbar(separatedCpData(i).Tm,separatedCpData(i).Cpm,separatedCpData(i).CpmFullErr,...
        '.','MarkerSize',18,'DisplayName',num2str(fields(i)/1e4,'%.2f T'),...
        'Color',clr{rng==i},'LineWidth',2);
end
xlabel('Temperature (K)'); ylabel('C$_p$ ($\mu$J/K)');
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
    plot(separatedCpData(i).T,separatedCpData(i).Cp,'.-')
    hold on
end
xlabel('Temperature (K)')
ylabel('Cp (uJ/K)')
legendCell = cellstr(num2str(fields, '%-d Oe'));
legend(legendCell)
hold off

%%
% Commented out on 2019-03-18. Remove if never used.
% x=separatedCpData(4).T;
% y=separatedCpData(4).Cp;

%%
figure
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
[xData, yData] = prepareCurveData( fields', Tc );

% Set up fittype and options.
ft = fittype( 'Tc0/Hc0*x/atanh(x/Hc0)', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', [7 8] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [5000 2.2];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
cft = confint(fitresult);
Hc0err = cft(2,1)-fitresult.Hc0;
Tc0err = cft(2,2)-fitresult.Tc0;

% Plot fit with data.
figure(); hold on
data = plot( xData(~excludedPoints), yData(~excludedPoints), '.', 'DisplayName', 'Data' );
excl = plot( xData(excludedPoints), yData(excludedPoints), 'xk',...
    'MarkerSize', 12, 'DisplayName', 'Excluded' );
xfit = [linspace(1,5000,50),fitresult.Hc0-logspace(2,-14)];
pfit = plot(xfit, fitresult(xfit), '-r', 'DisplayName', 'Fit'); 
xlim([0 6000]); ylim([0 2.5])
% lgd = legend( data, 'Tc vs. H', 'Excluded', 'Fit', 'Location', 'southwest' );
lgd = legend('Location', 'southwest' );

% Label axes
xlabel('$H_{\mathrm{ext}}$ (Oe)'); 
ylabel('$T$ (K)');
title('TmVO$_4$ phase boundary from $dC_p/dT$')
fitparam = annotation('textbox',[0.2 0.5 0.2 0.1],'interpreter','latex',...
    'String',{['Best fit parameters',...
    sprintf('\nTc0 = %.2f(%.1g) K', fitresult.Tc0, Tc0err*1e2),...
    sprintf('\nHc0 = %.3g(%.1g)e+3 Oe', fitresult.Hc0/10^3, Hc0err/10) ]},...
    'LineStyle','-','LineWidth',1,'FontSize',14,...
    'BackgroundColor','w','Color','k','FitBoxToText','on');%
grid on

%% Print fit parameter values with error bars
cval = coeffvalues(fitresult);% extract fit parameter values
cft=confint(fitresult);% extract confidence intervals from fit
sprintf("Hc(T=0) = %.1d +- %.0e Oe",fitresult.Hc0,cval(1)-cft(1,1))% print out value of critical field, with error bars
sprintf("Tc(H=0) = %.2f +- %.2f K",fitresult.Tc0,cval(2)-cft(1,2))% print out value of critical temperature, with error bars

%% Export figure
formatFigure;
% printPNG([todaystr '_TmVO4-2017-07-needle2_mfd@H=4750Oe']);
% printPDF(['2019-06-18_TmVO4-RF-E_fit_Schottky_' strrep(hrstr,'.','p') 'xHc']);







