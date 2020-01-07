%% Data importation
sample = 'TmVO4-RF-E';
filename = 'TmVO4_RF-E_2017-07-14.dat';
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2017-07_TmVO4_Cp_MCE\2017-07-20_Cp\2017-07-20_TmVO4_Cp_analysis'
DATA=ImportTmVO4Cp(filename);% Use this data to plot color map of phase diagram

%% Sample properties
m = 0.25e-3;% mass of sample, in g
M = 283.87;% molar mass of TmVO4, in g/mol
Tc = 2.126;% Value of transition temperature in this sample

%%
Hc0 = 0.49;% value in Tesla units of the critical field at zero temperature
% in the absence of demagnetizing factor
% see data taken on needles of TmVO4-LS5200 in July 2017
rescaling = Hc0/0.72;% rescaling factor, due to demag 
H=[DATA.FieldOersted]*rescaling;
T=[DATA.SampleTempKelvin];
Cp=[DATA.SampHCJmoleK];
CpErr=[DATA.SampHCErrJmoleK];

whichPoints = isfinite(H) & isfinite(T) & isfinite(Cp);
H=H(whichPoints);
T=T(whichPoints);
Cp=Cp(whichPoints).*M/m*1e-6;% factor 1e-6 converts from uJ/K to J/K
CpErr=CpErr(whichPoints).*M/m*1e-6;%

%% Plot 3D scatter of Cp data
[uh,~,X] = unique(round(H,-2));
hmax=max(uh);
tmin = min(T);
tmax = max(T);
figure
scatter3(H,T,Cp,'.')
xlim([0 hmax]);ylim([0 tmax]);
xlabel('Field (Oe)');ylabel('Temperature (K)');
zlabel('Cp (J/K/mol)')

%% Gaussian convolution
tstep=0.025; 
hstep=500*rescaling;
steps = -50:50;
x= tstep*steps;
s=0.05;
d1Gaussian = -exp(-x.^2/(2*s^2)).*x./sqrt(s^6*2*pi);
d2Gaussian = exp(-x.^2/(2*s^2)).*(x.^2 - s^2)/sqrt(s^10*2*pi);
d1Cp = conv2(Cp,d1Gaussian','same');

%% Plot 3D scatter of derivative of Cp
figure
scatter3(H,T,-d1Cp,'.')
xlim([0 hmax]);ylim([0 tmax]);
xlabel('Field (Oe)');ylabel('Temperature (K)');
zlabel('-dCp/dT (J/K/mol)')

%% Fit a 3D surface to the data
% [Hg,Tg] = meshgrid(unique(round(H,1)),unique(round(T,2)));% for use with griddata
[Hg,Tg] = meshgrid(0:hstep:hmax,tmin:tstep:tmax);% for use with gridfit
Hgl = 0:hstep:hmax;% for gridfit
Tgl = tmin:tstep:tmax;% for gridfit
figure
Cpg = gridfit(H,T,Cp,Hgl,Tgl);
% Cpg = griddata(H,T,Cp,Hg,Tg);
surf(Hg,Tg,Cpg)
hold on
scatter3(H,T,Cp,100,'.','.k')
ylim([0 tmax])
zlim([0 1.1*max(Cp)])
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('Cp (J/K/mol)')
%% Plot a 3D surface of the derivative of Cp
% Note: this is not a fit!
figure
d1Cpg = conv2(Cpg,d1Gaussian','same');
d2Cpg = conv2(Cpg,d2Gaussian','same');
surf(Hg./10000,Tg,-d1Cpg,'EdgeColor','none')
xlabel('Field (T)')
ylabel('Temperature (K)')
zlabel('-dCp/dT (J/K$^2$/mol)')
h=colorbar('eastoutside');
h.Label.String = '-dC_p/dT (J/K$^2$/mol)';
ylim([round(tmin,2) 3])
xlim([0 hmax/10^4])
view([0 90])

%% Plot 2D contour of derivative of Cp before averaging data
figure
n = 300;
contourf(Hg./5100,Tg,-d1Cpg,n,'EdgeColor','none');
hold on;
fplot(@(h)h/atanh(h),[0 1.1])
xlabel('H/H$_c$'); ylabel('Temperature (K)');
zlabel('-dCp/dT (J/K$^2$/mol)');
% title(sprintf('n = %i',n));

%% Separate data according to value of field
clear separatedCpData
for i = 1:length(uh)
    wp = abs(H-uh(i))<50;
    separatedCpData(i).H = H(wp);
    separatedCpData(i).T = T(wp);
    separatedCpData(i).Cp = Cp(wp);
    separatedCpData(i).d1Cp = d1Cp(wp);
    separatedCpData(i).CpErr = CpErr(wp);
    
    [separatedCpData(i).T,wo] = sort(separatedCpData(i).T);
    separatedCpData(i).H = separatedCpData(i).H(wo);
    separatedCpData(i).Cp = separatedCpData(i).Cp(wo);
    separatedCpData(i).d1Cp = separatedCpData(i).d1Cp(wo);
    separatedCpData(i).CpErr = separatedCpData(i).CpErr(wo);
end

%% Average datapoints taken at any given temperature and field setpoint
% We want to compute the average of data points that are repetitions of the
% same measurement, i.e. with same temperature and field setpoints
Tsep = 6e-3;% Data points taken within an interval of Tsep are considered to be measured at the same temperature setpoint
% 6mK is an empirical estimate of the dispersion of data points taken consecutively at a given temperature. 
for i = 1:length(uh)
    avgData(i) = averageCpwithH(Tsep,separatedCpData(i).T,separatedCpData(i).Cp,...
        separatedCpData(i).CpErr,separatedCpData(i).H);
end

%% Prepare fit of Cp vs Temperature 
i = 9;
% Use these variables in Curve fitting tool
clear fitT fitCp fitCpErr fitwghts 
fitT = avgData(i).T;
fitCp = avgData(i).Cp/R;
fitCpErr = avgData(i).CpFullErr/R;
fitwghts = 1./fitCpErr;

%% Compute Cp for Gaussian distribution of fields
h = uh(i)/(Hc0*1e4);
rhsgm = 0.08;
sgm = h*rhsgm;
T = linspace(0,1.4,701);
Cpid = zeros(length(T),1);
Cpreal = zeros(length(T),1);
for j=1:length(T)
    Cpid(j) = Cp_TFIM(T(j),h);
    Cpreal(j) = CpTFIM_normpdf(T(j),h,sgm);
end

%% Plot Cp for Gaussian distribution of fields
Tc = 2.15;
figure
plot(fitT,fitCp,'.','DisplayName','data')
hold on;
plot(T*Tc,Cpid,'DisplayName',sprintf('h=%.2f',h));
plot(T*Tc,Cpreal,'DisplayName',sprintf('h=%.2f,r=%.1e',h,rhsgm));
title(sprintf('Cp_TFIM vs CpTFIM_normpdf H=%.0fOe',uh(i)))
legend()

%% Compute additional Cp for Gaussian distribution of fields
h = 0.95;
rhsgm = 0.08;
sgm = h*rhsgm;
for j=1:length(T)
    Cpreal(j) = CpTFIM_normpdf(T(j),h,sgm);
end
plot(T*Tc,Cpreal,'DisplayName',sprintf('h=%.2f,r=%.1e',h,rhsgm));

%% Plot additional Cp for Gaussian distribution of fields
offset = 0.01;
plot(T*Tc,Cpreal+offset,'DisplayName',sprintf('h=%.2f,r=%.1e',h,rhsgm));



%% Compute and plot derivative of averaged data
figure; hold on
for i = 1:length(uh)
    avgData(i).d1Cp = conv2(avgData(i).Cp,d1Gaussian','same');
    plot(avgData(i).T,avgData(i).d1Cp)
end

%% Convert averaged data to table for exportation
avgtbl = table();
fn = fieldnames(avgData);
for i=1:numel(fn)
    avgtbl.(fn{i}) = cell2mat( arrayfun(@(c) c.(fn{i})(:), avgData(1:length(avgData)).', 'Uniform', 0) );
end
%% Add stuff to the table
% tbl.comments(1) = string(filename);
%% Export table of averaged data
%writetable(tbl,'2019-05-02_TmVO4-RF-E_Cp_avgData');

%% Create mesh for 2D contour
[Hgm,Tgm] = meshgrid(0:hstep:hmax,tmin:tstep:tmax);% for use with gridfit
Hgml = 0:hstep:hmax;% for gridfit
Tgml = tmin:tstep:tmax;% for gridfit
Cpgm = gridfit(H,T,Cp,Hgml,Tgml);
d1Cpgm = conv2(Cpgm,d1Gaussian','same');

%% Plot 2D contour of derivative of Cp
% Use this section to plot the phase diagram to combine with MCE traces
% from 'AnalyzeMCEinDR_TmVO4-LS5228-DR-HC180731.m'
figure
n = 300;
contourf(Hgm./5100,Tgm/Tc,-d1Cpgm,n,'EdgeColor','none');
hold on;
fplt = fplot(@(h)h/atanh(h),[0 1.1],'Color','k','LineWidth',1);
xlabel('$H / H_c(T=0)$'); ylabel('$T / T_D(H=0)$');
zlabel('-dCp/dT (J/K$^2$/mol)');
ylim([0.17 1.35]);
% title(sprintf('n = %i',n));
cb=colorbar('north'); cb.Ticks = 0:1000:2000;%make a horizontal colorbar at the top of the plot
cb.Label.String = '$-dC_p/dT$ (J/K$^2$/mol)'; cb.Label.Interpreter = 'latex'; cb.TickLabelInterpreter = 'latex';

%% Export d1Cp matrix to to a tab delimited file
filename = '2019-07-29_TmVO4-RF-E_dCp-dT.txt';
% Header
hdr1={'From file "2017-07_TmVO4-RF-E_Analyze_Cp_under_field.m"'};% First line header: quantities names
fmt1 = repmat('%s\t ', 1, length(hdr1)); fmt1(end:end+1) = '\n';% String formatting for header lines 1 and 2
fid = fopen(filename, 'wt');
fprintf(fid, fmt1, hdr1{:});% header 1
fclose(fid);
% dlmwrite(filename,-d1Cpgm,'-append','Delimiter','\t')

%% Export X and Y axes values of above matrix
Hgmr = Hgm(1,:)./5100;
Tgmr = Tgm(:,1)/Tc;
% Hgmr and Tgmr exported manually on 2019-07-29 by copying from Matlab 
% variables and pasting into worksheet of TmVO4_phase_diagram_Cp_MCE.opju


%% Plot errorbar plot of Hc(T) from xls file
% Need to import table from xls file first, e.g.
% '2019-05-21_TmVO4-LS5228-DR-HC180731_MCE_Hc_vs_T.xlsx'
% figure;% comment out this line when plotting on top of above colormap
ebdown = errorbar(hctbl.Hcrdown,hctbl.Tr,hctbl.dTr,hctbl.dTr,...
    hctbl.dHcrdown,hctbl.dHcrdown,'.m','MarkerSize',18,'LineWidth',2);
hold on;
ebup = errorbar(hctbl.Hcrup,hctbl.Tr,hctbl.dTr,hctbl.dTr,...
    hctbl.dHcrup,hctbl.dHcrup,'.g','MarkerSize',18,'LineWidth',2);
legend([ebup,ebdown],'$H_c^{\mathrm{min}}$','$H_c^{\mathrm{max}}$','Location','northeast');

%% Identify experimental critical temperature at each field
% Then we can plot Tc vs uh and fit using equation f(h) = Tc/Hc*h/atanh(h/Hc)
% From this, we can see that the experimental value of Hc is ~0.72T 
% and hence correct for demag, knowing the value of critical field Hc0 (see beginning of code)
M = ones(1,length(uh));
Tcd1 = ones(1,length(uh));
for i=1:length(uh)
    [M(i),I] = min(avgData(i).d1Cp);
    Tcd1(i) = avgData(i).T(I);
end

%% Prepare MF fit of Cp vs Temperature at given field
j=1;
Tfit = avgData(j).T;
Cpfit = avgData(j).Cp;
CpfitErr = avgData(j).CpFullErr;
wghts = 1./CpfitErr;
Tmaxfit = 2.15;

%% Fit and plot
[fitresult, gof] = fitCpmfvsTemp(Tfit,Cpfit,wghts,Tmaxfit);

%% Fit Cp data at H=0 and extract value of longitudinal strain
[fitresult, gof] = fitCpLFIM(Tfit,Cpfit/R,wghts*R,Tmaxfit)
% Results for Tmaxfit = 2.142; Tc = 2.128  (2.124, 2.132); e = 0.001947  (0.00129, 0.002604);
% Results for Tmaxfit = 2.15;% Tc = 2.126  (2.122, 2.129);% e = 0.001474 (0.001085, 0.001862);

%% Plot averaged data at each field separately
figure; hold on
rng = [1 5 7 9 13];
clr = lines(length(rng));
eb = cell(size(rng));
% for i=1:length(rng)
%     fp = plot(Ttlf(2:end-1)*Tc,Cptlf(:,i));% requires computing Cptlf in 'F_S_Cp_TLFIM_compute.m'
% end
for i=rng
    fp = fplot(@(t)Cp_TFIM(t/2.14,uh(i)/(5.1e3)),[0 3.2],'LineWidth',2);
%     fp = fplot(@(t)Cp_TFIM_offset_strain(t/2.125,uh(i)/(5.1e3),1.5e-3),[0 4],'LineWidth',2);
% Fit parameters on data at H=0: Tc=2.125(3), e=1.5(4)e-3
% Note: the values of amplitude coefficient and Tc extracted from fit 
% in curve fitting tool using Cp_TFIM (no offset strain) are A=7.35 and Tc=2.142K
%     clr{rng==i} = get(fp,'Color');
end
for i=rng
    eb{rng==i} = errorbar(avgData(i).T,avgData(i).Cp/R,avgData(i).CpFullErr/R,...
        '.','MarkerSize',18,'DisplayName',num2str(uh(i)/(Hc0*1e4),'%.2f'),...
        'Color',clr(rng==i,:),'LineWidth',2);
end
xlabel('Temperature (K)'); ylabel('C$_p$/R');%ylabel('C$_p$ (JK$^{-1}$mol$^{-1}$)');
% title('Heat capacity of TmVO4 at various fields')
% title([sprintf('h=H/Hc*%.3g, e=',factor) mat2str(e,2)])
lgd = legend([eb{:}]); lgd.Title.String = '$H/H_c$';
% legendCell = cellstr(num2str(uh, '%-d Oe')); legend(legendCell)
ax = gca; ax.YMinorTick = 'on';% Add minor ticks on Y axis
grid on;%
hold off

%% Export figure 
formatFigure
printPNG('2020-01-07_TmVO4-RF-E_CpTFIM_normpdf_5000Oe')
% printPDF('2019-08-21_TmVO4-RF-E_Cp_vs_T_5H_theory_e-propto-h-cube')

%% Prepare MF fit of Cp vs Temperature under field
index = 15;
H1 = uh(index);
T1 = avgData(index).T;
Cp1 = avgData(index).Cp/R;
Cp1Err = avgData(index).CpFullErr/R;
wghts1 = 1./Cp1Err;
%% Fit and plot
maxTfit = 2;%Kelvin
hrstr = sprintf('%.2f',H1/(Hc0*1e4));
% [fitresult1, gof1] = fitCpTFIM(T1,Cp1,wghts1,2.142,H1/5.1e3);
[fitresult1, gof1] = fitSchTemp(T1,Cp1,wghts1,maxTfit);
fitprms = coeffvalues(fitresult1);
prmerr = coeffvalues(fitresult1)-confint(fitresult1);
gapstr = [sprintf('%.2f',fitprms(1)) '$\pm$' sprintf('%.2f',prmerr(1,1))];
offsetstr = [sprintf('%.3f',fitprms(2)) '$\pm$' sprintf('%.3f',prmerr(1,2))];
title({['Schottky anomaly fit of ' sample],[' at $H/H_c=$' hrstr]});
annfit = annotation('textbox',[0.575 0.175 0.2 0.1],'interpreter','latex',...
    'String',{['$R^2=$ ' sprintf('%.4f',gof1.rsquare)] ['$\Delta=$ ' gapstr 'K']...
    ['offset ' offsetstr]}, 'LineStyle','-','EdgeColor','k',...
    'FitBoxToText','on','LineWidth',1,'BackgroundColor','w','Color','k');% add annotation
formatFigure
annfit.Position(2)=.175;
%% Export figure to pdf
% printPDF(['2019-06-18_TmVO4-RF-E_fit_Schottky_' strrep(hrstr,'.','p') 'xHc']);



%% 3D scatter of Cp(T) at each field separately
figure
for i=1:length(uh)
    scatter3(separatedCpData(i).H,separatedCpData(i).T,separatedCpData(i).Cp,'.',...
        'Displayname',num2str(uh(i)', '%d Oe'))
%     plot(separatedCpData(i).T,separatedCpData(i).Cp,'-+')
    hold on
end
xlabel('H (Oe)');ylabel('T (K)');zlabel('C_p (\muJ/K)');
xlim([0 hmax])
legend()
hold off
%% Plot dCp/dT at each field separately
% Need to average data points measured at each temperature first

figure
for i=1:length(uh)
    scatter3(separatedCpData(i).H,separatedCpData(i).T,-separatedCpData(i).d1Cp,'.',...
        'Displayname',num2str(uh(i)', '%d Oe'))
%     plot(separatedCpData(i).T,separatedCpData(i).Cp,'-+')
    hold on
end
xlabel('H (Oe)');ylabel('T (K)');zlabel('-dC_p/dT (\muJ/K)');
xlim([0 hmax])
legend()
hold off
%% Plot Cp(T) at each field separately
%%
figure
for i=1:length(uh)
    plot(separatedCpData(i).T,separatedCpData(i).Cp,'-+')
    hold on
end
xlabel('Temperature (K)')
ylabel('Cp/T (J/K^2/mol)')
legendCell = cellstr(num2str(uh, '%-d Oe'));
legend(legendCell)
hold off
