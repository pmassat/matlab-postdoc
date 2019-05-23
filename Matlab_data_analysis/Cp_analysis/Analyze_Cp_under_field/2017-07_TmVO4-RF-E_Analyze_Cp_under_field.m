%% Data importation
filename = 'TmVO4_RF-E_2017-07-14.dat';
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2017-07_TmVO4_Cp_MCE\2017-07-20_Cp\2017-07-20_TmVO4_Cp_analysis'
DATA=ImportTmVO4Cp(filename);% Use this data to plot color map of phase diagram
%% Sample properties
m = 0.25e-3;% mass of sample, in g
M = 283.87;% molar mass of TmVO4, in g/mol
R = 8.314;% gas constant, in J/K/mol
%%
Hc0 = 0.51;% value in Tesla units of the critical field at zero temperature
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
%     T2 = repmat(separatedCpData(i).T,1);
%     Tm = zeros(length(separatedCpData(i).T),1);
%     Tsd = zeros(length(separatedCpData(i).T),1);
%     Cpm = zeros(length(separatedCpData(i).Cp),1);
%     Cpsd = zeros(length(separatedCpData(i).Cp),1);
%     CpmErr = zeros(length(separatedCpData(i).Cp),1);
%     for k = 1:length(T2)
%         if T2(k)==0
%             continue
%         elseif length(T2(abs(T2-T2(k))<Tsep))>3
%             Tsep2 = Tsep/2;% reduce the temperature interval
% %             T2(abs(T2-T2(k))<Tsep2)%print out values of temperature which
% %             verify the if statement
%             Tm(k) = mean(T2(abs(T2-T2(k))<Tsep2));
%             Tsd(k) = std(T2(abs(T2-T2(k))<Tsep2));
%             Cpm(k) = mean(separatedCpData(i).Cp(abs(T2-T2(k))<Tsep2));
%             Cpsd(k) = std(separatedCpData(i).Cp(abs(T2-T2(k))<Tsep2));
%             CpmErr(k) = sum(separatedCpData(i).CpErr(abs(T2-T2(k))<Tsep2))/...
%                 sqrt(length(separatedCpData(i).CpErr(abs(T2-T2(k))<Tsep2)));
%             T2(abs(T2-T2(k))<Tsep2)=0;
%         else
%             Tm(k) = mean(T2(abs(T2-T2(k))<Tsep));
%             Tsd(k) = std(T2(abs(T2-T2(k))<Tsep));
%             Cpm(k) = mean(separatedCpData(i).Cp(abs(T2-T2(k))<Tsep));
%             Cpsd(k) = std(separatedCpData(i).Cp(abs(T2-T2(k))<Tsep));
%             CpmErr(k) = sum(separatedCpData(i).CpErr(abs(T2-T2(k))<Tsep))/...
%                 sqrt(length(separatedCpData(i).CpErr(abs(T2-T2(k))<Tsep)));
%             T2(abs(T2-T2(k))<Tsep)=0;
%         end
%     end
%     separatedCpData(i).Hm = separatedCpData(i).H(Tm>0);
%     separatedCpData(i).Tm = Tm(Tm>0);
%     separatedCpData(i).Tsd = Tsd(Tm>0);
%     separatedCpData(i).Cpm = Cpm(Cpm>0);
%     separatedCpData(i).CpmFullErr = Cpsd(Cpm>0) + CpmErr(Cpm>0);
    avgData(i) = averageCpwithH(Tsep,separatedCpData(i).T,separatedCpData(i).Cp,...
        separatedCpData(i).CpErr,separatedCpData(i).H);
end

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
% from 'AnalyzeMCEinDR_TmVO4-LS5228-DR-HC180731_2018-09-19.m'
figure
n = 300;
Tc = 2.14;% Value of transition temperature in this sample
contourf(Hgm./5100,Tgm/Tc,-d1Cpgm,n,'EdgeColor','none');
hold on;
fplot(@(h)h/atanh(h),[0 1.1],'Color','k')
xlabel('H/H$_c$'); ylabel('T/T$_D$');
zlabel('-dCp/dT (J/K$^2$/mol)');
% title(sprintf('n = %i',n));
% h=colorbar('eastoutside');

%% Plot errorbar plot of Hc(T) from xls file
% Need to import table from xls file first
% figure;% comment out this line when plotting on top of above colormap
ebup = errorbar(hctbl.Hcrup,hctbl.Tr,hctbl.dTr,hctbl.dTr,...
    hctbl.dHcrup,hctbl.dHcrup,'.g','MarkerSize',12,'LineWidth',1);
hold on;
ebdown = errorbar(hctbl.Hcrdown,hctbl.Tr,hctbl.dTr,hctbl.dTr,...
    hctbl.dHcrdown,hctbl.dHcrdown,'.r','MarkerSize',12,'LineWidth',1);
legend([ebup,ebdown],'$H_c^{\mathrm{min}}$','$H_c^{\mathrm{max}}$','Location','northeast');

%% Identify experimental critical temperature at each field
% Then we can plot Tc vs uh and fit using equation f(h) = Tc/Hc*h/atanh(h/Hc)
% From this, we can see that the experimental value of Hc is ~0.72T 
% and hence correct for demag, knowing the value of critical field Hc0 (see beginning of code)
M = ones(1,length(uh));
Tc = ones(1,length(uh));
for i=1:length(uh)
    [M(i),I] = min(avgData(i).d1Cp);
    Tc(i) = avgData(i).T(I);
end

%% Plot averaged data at each field separately
figure; hold on
rng = 1:4:length(uh);
clr = cell(size(rng));
eb = cell(size(rng));
for i=rng
    fp = fplot(@(t)Cp_TFIM_offset_strain(t/2.14,0,uh(i)/(5.1e3)),[0 4],'LineWidth',2);
%     fp = fplot(@(t)Cp_TFIM_offset_strain(t/2.125,1.5e-3,uh(i)/(5.1e3)),[0 4],'LineWidth',2);
% Fit parameters on data at H=0: Tc=2.125(3), e=1.5(4)e-3
% Note: the values of amplitude coefficient and Tc extracted from fit 
% in curve fitting tool using the Cp_TFIM (no offset strain) are A=7.35 and Tc=2.142K
    clr{rng==i} = get(fp,'Color');
end
for i=rng
    eb{rng==i} = errorbar(avgData(i).T,avgData(i).Cp/R,avgData(i).CpFullErr/R,...
        '.','MarkerSize',18,'DisplayName',num2str(uh(i)/(Hc0*1e4),'%.2f'),...
        'Color',clr{rng==i},'LineWidth',2);
end
xlabel('Temperature (K)'); ylabel('C$_p$/R');%ylabel('C$_p$ (JK$^{-1}$mol$^{-1}$)');
% title('Heat capacity of TmVO4 at various fields')
lgd = legend([eb{:}]); lgd.Title.String = '$H/H_c$';
% legendCell = cellstr(num2str(uh, '%-d Oe'));
% legend(legendCell)
hold off

%% Export figure 
% printPNG('2019-05-17_TmVO4-RF-E_Cp_vs_T_4H_+fits_No-strain')
% printPDF('2019-05-17_TmVO4-RF-E_Cp_vs_T_4H_+fits')

%% Prepare MF fit of Cp vs Temperature at zero field
T0 = avgData(1).T;
Cp0 = avgData(1).Cp;
Cp0Err = avgData(1).CpFullErr;
wghts = 1./Cp0Err;
%% Fit and plot
[fitresult, gof] = fitCpmfvsTemp(T0,Cp0,wghts,2.142);

%% Prepare MF fit of Cp vs Temperature under field
index = 5;
H1 = uh(index);
T1 = avgData(index).T;
Cp1 = avgData(index).Cp;
Cp1Err = avgData(index).CpFullErr;
wghts1 = 1./Cp1Err;
%% Fit and plot
[fitresult1, gof1] = fitCpTFIM(T1,Cp1,wghts1,2.142,H1/5.1e3);





%% 
%%
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
%% Plot Cp/T(T)
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
%%
whichT = 1.8:0.1:6.1;
clear cstTCpdata
for i=length(whichT):-1:1
    if length(T(abs(T-whichT(i))<0.005)) < 4
        whichT(i) = [];
    end
end

for i=1:length(whichT)
    cstTfilter = abs(T-whichT(i))<0.01;
    cstTCpdata(i).H = H(cstTfilter);
    cstTCpdata(i).T = T(cstTfilter);
    cstTCpdata(i).Cp = Cp(cstTfilter);
    
    [cstTCpdata(i).H,wo] = sort(cstTCpdata(i).H);
    cstTCpdata(i).H = cstTCpdata(i).H(wo);
    cstTCpdata(i).Cp = cstTCpdata(i).Cp(wo);
end
%% 
% For each value of temperature rounded to 0.1K (whichT), identify the closest 
% temperature value measured and store this value and the corresponding Cp value
%% Plot Cp(H) at constant T
%%
index = 1:length(whichT);
for i=index
    plot(cstTCpdata(i).H,cstTCpdata(i).Cp,'-+')
    hold on
end
xlabel('Field (Oe)')
ylabel('Cp (J/K/mol)')
legendCell = cellstr(num2str(whichT(index)', '%1.1f K'));
legend(legendCell)
hold off
%% Plot Cp/T(H) at constant T
%%
index = 1:length(whichT);
for i=index
    plot(cstTCpdata(i).H,cstTCpdata(i).Cp./cstTCpdata(i).T,'-+')
    hold on
end
xlabel('Field (Oe)')
ylabel('Cp/T (J/K^2/mol)')
legendCell = cellstr(num2str(whichT(index)', '%1.1f K'));
legend(legendCell)
hold off
%% 
% Next steps: 
% 
% # save file in TmAsO4 folder
% # plot Cp/T = f(T)
% # Calculate Cp = f(H) at fixed temperature
% 
% 
%% Below this point: code not adapted yet
%%
x=separatedCpData(4).T;
y=separatedCpData(4).Cp;
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
ylabel('-dCp/dT (J/K^2/mol)')
hold off
%%
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
    scatter3(separatedCpData(i).H,separatedCpData(i).T,-separatedCpData(i).d2,'filled')
    hold on
end
xlim([0 hmax])
ylim([0 3])
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('-dCp/dT (J/K/mol)')
hold off
%%
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
% [Hg,Tg] = meshgrid(0:hstep:hmax,tmin:tstep:3);%for griddata
Hgl=0:hstep:hmax;%for gridfit
Tgl=tmin:tstep:3;%for gridfit
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