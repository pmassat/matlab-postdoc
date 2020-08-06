%% Sample properties
m = 0.25e-3;% mass of sample, in g
M = 283.87;% molar mass of TmVO4, in g/mol
Tc0 = 2.126;% Value of transition temperature in this sample, in Kelvin units
Tcnum = 2.15;% in Kelvin units
Hc = 5100;% in Oersted units

%% Import magnetic field distribution from TXT file
% cd 'C:\Users\Pierre\Desktop\Postdoc\Software\COMSOL\TmVO4-RF-E_HC2017-07'
% mfdRFfname = '2020-07-21_TmVO4-RF-E_zero-temp-magnetization_mag-field_distrib.txt';
% Smfd = importfielddistrib(mfdRFfname, 15);

%% Import magnetic field distribution from CSV file
cd 'C:\Users\Pierre\Desktop\Postdoc\Software\COMSOL\TmVO4-RF-E_HC2017-07'
mfdRFfname = cell(1,2);
Srf = cell(1,2);
mfdRFfname{1} = '2020-07-27_TmVO4-RF-E_COMSOL_mfd_mesh=finer-out+max20um-in_T=p001-1K_H=all.csv';
mfdRFfname{2} = '2020-07-27_TmVO4-RF-E_COMSOL_mfd_mesh=finer-out+max20um-in_T=2-3K_H=all.csv';
for ic=1:length(Srf)
    Srf{ic} = importfielddistrib_csv(mfdRFfname{ic}, 28);
end

%% Extract values of temperature and external magnetic field from structure header
Smfd_RF = [Srf{1,:}];% concatenate cell arrays into structure
for sidx=1:length(Smfd_RF)
    % split each header into a cell array of strings, using whitespace and '=' sign as separators
    TBextCell = strsplit(Smfd_RF(sidx).T_Bext,{' ','='});
    % Find the index of the string cell containing "T" (temperature), and
    % convert the following string cell (which contains the value of temperature) into a number
    % Note: circshift(A,1) shifts the index of elements of array A by 1 to the right
    Smfd_RF(sidx).T_K = str2double(TBextCell{circshift(TBextCell=="T",1)});
    % Same with "Bext", which is the external magnetic flux density, in Tesla
    % units, and convert it into a value of magnetic field, in Oersted units
    Smfd_RF(sidx).Hext_Oe = str2double(TBextCell{circshift(TBextCell=="Bext",1)})*10^4;
end

%% Compute probability distribution of fields at a given value of T and Hext
for i=1:length(Smfd_RF)
mfd = Smfd_RF(i).mfd(Smfd_RF(i).mfd>0)/Hc;% create distribution from non-zero values
% h = histogram(mfd, 'Normalization', 'pdf');% plot histogram of distribution
[Smfd_RF(i).hc, edges] = histcounts(mfd, 50, 'Normalization', 'pdf');% plot histogram of distribution
% binCenters = h.BinEdges + (h.BinWidth/2);
Smfd_RF(i).binCenters = mean([edges(1:end-1);edges(2:end)],1);
Smfd_RF(i).binWidths = edges(2:end)-edges(1:end-1);
Smfd_RF(i).label = 'RF-E';
end

%% Plot distribution of fields at a given value of T and Hext
figure;
hold on
param_index = 1;% 1 is constant T, 2 is constant Hext, see param_range
temp_index = 1;% determines temperature of data to plot, taken among [1e-3, 1, 2, 3]K
field_index = 4;% determines field of data to plot from [.1,.2,.3,.4,.45,.5,.55,.6,.62,.63,.65,.67,.7,.8]T
param_range = {temp_index*14+[1:2:14], [field_index:14:56]};% first range corresponds to a 
% field dependence at constant temp, second range corresponds to a 
% temperature dependence at constant field
rng = param_range{param_index};
for i=rng
T = Smfd_RF(i).T_K;%
Hext = Smfd_RF(i).Hext_Oe;%
p = plot(Smfd_RF(i).binCenters, Smfd_RF(i).hc, '.-', 'DisplayName', sprintf('%.2g, %.2g',T/Tcnum,Hext/Hc));
end
lgd = legend('show'); lgd.Title.String = '$T/T_c$, $H_{\mathrm{ext}}/H_c$';
param_title = {[', $T=$' sprintf('%.2g K',T)],...
    [', $H_{\mathrm{ext}}=$' sprintf('%.2d Oe',Hext)]};
title(['Distribution of fields in TmVO$_4$-RF-E' param_title{param_index}])
xlabel('$H_{\mathrm{in}}/H_{\mathrm{c}}$')
ylabel('Normalized PDF')

%% Compute the relative difference between curves at 1000 Oe and 7000 Oe
mfd_diff(113).abs_diff = abs(Smfd_RF(1).hc);%-Smfd(13).hc);
% mfd_diff(113).rel = mfd_diff(113).abs/Smfd(1).binWidths(1);

%% Data importation
sample = 'TmVO4-RF-E';
filename = 'TmVO4_RF-E_2017-07-14.dat';
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2017-07_TmVO4_Cp_MCE\2017-07-20_Cp\2017-07-20_TmVO4_Cp_analysis'
RFDATA=ImportTmVO4Cp(filename);% Use this data to plot color map of phase diagram

%% Assign data to variables
Hc0 = 0.51;% value in Tesla units of the critical field at zero temperature
% in the absence of demagnetizing factor
% see data taken on needles of TmVO4-LS5200 in July 2017
rescaling = Hc0/0.69;% rescaling factor, due to demag; 
% obtained from fit using normal PDF of fields close to Hc; see below
H=[RFDATA.FieldOersted];%*rescaling;
T=[RFDATA.SampleTempKelvin];
Cp=[RFDATA.SampHCJmoleK];
CpErr=[RFDATA.SampHCErrJmoleK];

whichPoints = isfinite(H) & isfinite(T) & isfinite(Cp);
H=H(whichPoints);
T=T(whichPoints);
Cp=Cp(whichPoints).*M/m*1e-6;% factor 1e-6 converts from uJ/K to J/K
CpErr=CpErr(whichPoints).*M/m*1e-6;%
[uhrf,~,X] = unique(round(H,-1));
uhrf(2,:)=[];% When rounding to nearest tens, remove uh=10Oe because it is redundant with uh=0Oe

%% Plot 2D scatter of Cp data at H=0
figure
plot(T(round(H,-1)<20),Cp(round(H,-1)<20),'.')
% xlim([0 hmax]);ylim([0 tmax]);
xlabel('Temperature (K)');
ylabel('Cp (J/K/mol)')

%% Plot 3D scatter of Cp data
hmax=max(uhrf);
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
for i = 1:length(uhrf)
    wp = abs(H-uhrf(i))<50;
    separatedRFCpData(i).H = H(wp);
    separatedRFCpData(i).T = T(wp);
    separatedRFCpData(i).Cp = Cp(wp);
    separatedRFCpData(i).d1Cp = d1Cp(wp);
    separatedRFCpData(i).CpErr = CpErr(wp);
    
    [separatedRFCpData(i).T,wo] = sort(separatedRFCpData(i).T);
    separatedRFCpData(i).H = separatedRFCpData(i).H(wo);
    separatedRFCpData(i).Cp = separatedRFCpData(i).Cp(wo);
    separatedRFCpData(i).d1Cp = separatedRFCpData(i).d1Cp(wo);
    separatedRFCpData(i).CpErr = separatedRFCpData(i).CpErr(wo);
end

%% Average datapoints taken at any given temperature and field setpoint
% We want to compute the average of data points that are repetitions of the
% same measurement, i.e. with same temperature and field setpoints
clear avgRFData
Tsep = 6e-3;% Data points taken within an interval of Tsep are considered to be measured at the same temperature setpoint
% 6mK is an empirical estimate of the dispersion of data points taken consecutively at a given temperature. 
Tp = 27;% Temperature scale of phonons contribution to Cp in TmVO4, in K; see 'TmVO4_Cp_phonons.m'
for i=length(uhrf):-1:1
avgRFData(i) = averageCpwithH2(Tsep,separatedRFCpData(i).T,separatedRFCpData(i).Cp,...
    separatedRFCpData(i).CpErr,separatedRFCpData(i).H);
end
for i=length(uhrf):-1:1
avgRFData(i).Cpel = avgRFData(i).Cp - R*(avgRFData(i).T/Tp).^3;% electronic contribution to Cp, after subtracting phonons contribution
avgRFData(i).Cpelr = avgRFData(i).Cpel/R;
avgRFData(i).CpelrErr = avgRFData(i).CpFullErr/R;
avgRFData(i).uh = unique(round(avgRFData(i).H,-2));
end

%% Plot averaged data 
i = 1;
figure
plot(avgRFData(i).T,avgRFData(i).Cp,'.','DisplayName','$C_p^{\mathrm{full}}$')
hold on
plot(avgRFData(i).T,R*(avgRFData(i).T/Tp).^3,'.','DisplayName','$C_p^{\mathrm{phonons}}$')
plot(avgRFData(i).T,avgRFData(i).Cpel,'.','DisplayName','$C_p^{\mathrm{full}}-C_p^{\mathrm{phonons}}$')
legend('show')
title('TmVO$_4$ heat capacity')
xlabel('T (K)'); ylabel('$C_p$ (J/K/mol)');

%% Compute and plot derivative of averaged data
figure; hold on
for i = 1:length(uhrf)
    avgRFData(i).d1Cp = conv2(avgRFData(i).Cp,d1Gaussian','same');
    plot(avgRFData(i).T,avgRFData(i).d1Cp)
end

%% Convert averaged data to table for exportation
avgtbl = table();
fn = fieldnames(avgRFData);
for i=1:numel(fn)
    avgtbl.(fn{i}) = cell2mat( arrayfun(@(c) c.(fn{i})(:), avgRFData(1:length(avgRFData)).', 'Uniform', 0) );
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
contourf(Hgm./5100,Tgm/Tc0,-d1Cpgm,n,'EdgeColor','none');
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
Tgmr = Tgm(:,1)/Tc0;
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
M = ones(1,length(uhrf));
Tcd1 = ones(1,length(uhrf));
for i=1:length(uhrf)
    [M(i),I] = min(avgRFData(i).d1Cp);
    Tcd1(i) = avgRFData(i).T(I);
end

%% Prepare MF fit of Cp vs Temperature at given field
j=1;
Tfit = avgRFData(j).T;
Cpfit = avgRFData(j).Cp;
CpfitErr = avgRFData(j).CpFullErr;
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
rng = [1:2:13];% [1 5 7 9 13];
clr = lines(length(rng));
eb = cell(size(rng));
% for i=1:length(rng)
%     fp = plot(Ttlf(2:end-1)*Tc,Cptlf(:,i));% requires computing Cptlf in 'F_S_Cp_TLFIM_compute.m'
% end
maxTplot = 3.2;%
for i=rng
    fp = fplot(@(t)Cp_TFIM(t/Tc0,uhrf(i)*rescaling/(Hc0*1e4)),[0 maxTplot],'LineWidth',2);
%     fp = fplot(@(t)Cp_TFIM_offset_strain(t/2.125,uh(i)/(5.1e3),1.5e-3),[0 4],'LineWidth',2);
% Fit parameters on data at H=0: Tc=2.125(3), e=1.5(4)e-3
% Note: the values of amplitude coefficient and Tc extracted from fit 
% in curve fitting tool using Cp_TFIM (no offset strain) are A=7.35 and Tc=2.142K
%     clr{rng==i} = get(fp,'Color');
end
for i=rng
    eb{rng==i} = errorbar(avgRFData(i).T,avgRFData(i).Cpelr,avgRFData(i).CpFullErr/R,...
        '.','MarkerSize',18,'DisplayName',num2str(uhrf(i)/(Hc0*1e4),'%.2f'),...
        'Color',clr(rng==i,:),'LineWidth',2);
end
xlabel('Temperature (K)'); ylabel('C$_p$/R');%ylabel('C$_p$ (JK$^{-1}$mol$^{-1}$)');
xlim([0 maxTplot])
% title('Heat capacity of TmVO4 at various fields')
% title([sprintf('h=H/Hc*%.3g, e=',factor) mat2str(e,2)])
lgd = legend([eb{:}]); lgd.Title.String = '$H_{\mathrm{ext}}/H_c$';
% legendCell = cellstr(num2str(uh, '%-d Oe')); legend(legendCell)
ax = gca; ax.YMinorTick = 'on';% Add minor ticks on Y axis
grid on;%
hold off

%% Export figure 
formatFigure; 
% printPNG('2020-01-08_TmVO4-RF-E_CpTFIM_normpdf_4580Oe')
% printPDF('2019-08-21_TmVO4-RF-E_Cp_vs_T_5H_theory_e-propto-h-cube')

%% Prepare fit of Cp vs Temperature 
i = 9;
% Use these variables in Curve fitting tool
clear fitT fitCp fitCpErr fitwghts 
fitT = avgRFData(i).T;
fitCp = avgRFData(i).Cpelr;
fitCpErr = avgRFData(i).CpFullErr/R;
fitwghts = 1./fitCpErr;
% Tc = 2.15;

%% Compute Cp for Gaussian distribution of fields
clear Cpnum
for i=rng
Cpnum(i).h = uhrf(i)*rescaling/Hc;
Cpnum(i).rhsgm = 0.09;
Cpnum(i).sgm = Cpnum(i).h*Cpnum(i).rhsgm;
Cpnum(i).t_single_h = linspace(0,1.5,601);% reduced temperature, T/Tc
Cpnum(i).single_h = zeros(size(Cpnum(i).t_single_h));
Cpnum(i).t_phenomeno = linspace(0,1.5,301);% reduced temperature, T/Tc
Cpnum(i).normpdf = zeros(size(Cpnum(i).t_phenomeno));
end
i = 1;
Cpnum(i).single_h = Cp_TFIM(Cpnum(i).t_single_h,Cpnum(i).h);
Cpnum(i).normpdf = Cp_LFIM(Cpnum(i).t_phenomeno,1.5e-3);
% Fit parameters on data at H=0: Tc=2.125(3), e=1.5(4)e-3
% Cp_LFIM(h=0)
for i=rng(2:end)
Cpnum(i).single_h = Cp_TFIM(Cpnum(i).t_single_h,Cpnum(i).h);
    for j=2:length(Cpnum(i).t_phenomeno)
    Cpnum(i).normpdf(j) = CpTFIM_normpdf(Cpnum(i).t_phenomeno(j),Cpnum(i).h,Cpnum(i).sgm);
    end
end

%% Compute Cp for COMSOL distribution of fields
for i=1
Cpnum(i).comsolpdf = zeros(size(Cpnum(i).t_phenomeno));
end

%% Gaussian function (see normpdf.m)
mu = .58; 
sigma =.09*mu; 
gauss = @(x) exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);

%% Compute probability distribution of fields at a given value of T and Hext
Hcnum = 5500;
for mfdidx=4:14:56
mfd = Smfd_RF(mfdidx).mfd(Smfd_RF(mfdidx).mfd>0)/Hcnum;% create distribution from non-zero values
% h = histogram(mfd, 'Normalization', 'pdf');% plot histogram of distribution
[Smfd_RF(mfdidx).hc, edges] = histcounts(mfd, 100, 'Normalization', 'pdf');% plot histogram of distribution
% binCenters = h.BinEdges + (h.BinWidth/2);
Smfd_RF(mfdidx).binCenters = mean([edges(1:end-1);edges(2:end)],1);
Smfd_RF(mfdidx).binWidths = edges(2:end)-edges(1:end-1);
end

%% Create table from structure
Tmfd_RF = struct2table(Smfd_RF);% 
utmfd = unique(Tmfd_RF.T_K);
uhmfd = unique(Tmfd_RF.Hext_Oe);

%% For a given dataset, find closest values of temperature and field in COMSOL mfd
i = 5;
Hdata = unique(round(avgRFData(i).H,-2));
[~,mfdhidx] = min(abs(Tmfd_RF.Hext_Oe-Hdata));
h = Tmfd_RF.Hext_Oe(mfdhidx);
tref = 0;
for j=1:length(Cpnum(i).t_phenomeno)
    % Find value of temperature in COMSOL mfd closest to that of actual data
    [~,mfdtidx] = min(abs(Tmfd_RF.T_K/Tcnum-Cpnum(i).t_phenomeno(j)));
    % Improvement note: use sort instead of min, to be able to interpolate...
    t = Tmfd_RF.T_K(mfdtidx);
    if t ~= tref
        sprintf('j=%i, T=%.2gK',j,t)
        tref=t;
    end
    % Same for value of field
    % Find the row in Tmfd that matches both t and h 
    row = find(Tmfd_RF.T_K==t & Tmfd_RF.Hext_Oe==h);
    Cph = zeros(size(Tmfd_RF.binCenters(row,:)));
    for col=1:length(Tmfd_RF.binCenters(row,:))
        Cph(col) = Cp_TFIM(Cpnum(i).t_phenomeno(j),Tmfd_RF.binCenters(row,col));
    end
    % Compute the corresponding value of heat capacity 
%     Cpnum(i).comsolpdf(j) = trapz(Tmfd_RF.binCenters(row,:),Cph.*Tmfd_RF.hc(row,:));
    Cpnum(i).comsolpdf(j) = sum(Tmfd_ndl.binWidths(row,:).*Cph.*...
        Tmfd_ndl.hc(row,:));
end

%% Plot Cp for COMSOL distribution of fields
figure
plot(avgRFData(i).T,avgRFData(i).Cpelr,'.','DisplayName','data')
hold on;
plot(Cpnum(i).t_single_h*Tc0,Cpnum(i).single_h,'DisplayName','MF');
plot(Cpnum(i).t_phenomeno*Tc0,Cpnum(i).comsolpdf,'DisplayName',sprintf('Hc=%.2dOe',Hcnum));
title(['Cp mean-field vs COMSOL pdf $H_{\mathrm{ext}}=$' sprintf('%.0fOe',uhrf(i))]);
lgd = legend();% title(lgd,'TmVO4-RF-E');

%% Plot Cp for Gaussian distribution of fields
figure
plot(avgRFData(i).T,avgRFData(i).Cpelr,'.','DisplayName','data')
hold on;
plot(Cpnum(i).t_single_h*Tc0,Cpnum(i).single_h,'DisplayName','MF');
plot(Cpnum(i).t_phenomeno*Tc0,Cpnum(i).normpdf,'DisplayName','Gauss');
title(['Cp mean-field vs normal pdf $H_{\mathrm{ext}}=$' sprintf('%.0fOe',uhrf(i))]);
lgd = legend();% title(lgd,'TmVO4-RF-E');

%% Plot averaged data at each field separately
figure; hold on
clr = lines(length(rng));
eb = cell(size(rng));
for i=rng
% fp = fplot(@(t)Cp_TFIM(t/Tc0,Cptheo(i).h),[0 3.2],'--','LineWidth',2,'Color',clr(rng==i,:));
plot(Cpnum(i).t_single_h*Tc0,Cpnum(i).single_h,'--','Color',clr(rng==i,:),'DisplayName',sprintf('h=%.2f',Cpnum(i).h));
plot(Cpnum(i).t_phenomeno*Tc0,Cpnum(i).normpdf,'Color',clr(rng==i,:),'DisplayName',...
    sprintf('h=%.2f,r=%.1e',Cpnum(i).h,Cpnum(i).rhsgm));
end
for i=rng
eb{rng==i} = errorbar(avgRFData(i).T,avgRFData(i).Cpelr,avgRFData(i).CpFullErr/R,...
    '.','MarkerSize',18,'DisplayName',num2str(uhrf(i)/(Hc0*1e4),'%.2f'),...
    'Color',clr(rng==i,:),'LineWidth',2);
end
xlabel('$T$ (K)'); ylabel('$C_p/R$');%ylabel('C$_p$ (JK$^{-1}$mol$^{-1}$)');
xlim([0 max(Cpnum(rng(1)).t_single_h*Tc0)]);
lgd = legend([eb{:}]); lgd.Title.String = '$H/H_c$';
ax = gca; ax.YMinorTick = 'on';% Add minor ticks on Y axis
anntheo = annotation('textbox',[0.13 0.83 0.2 0.1],'interpreter','latex',...
    'String',{['$--$ 1-parameter fit'] ['----- 2-parameter fit']...
    }, 'LineStyle','-','EdgeColor','k',...
    'FitBoxToText','on','LineWidth',1,'BackgroundColor','w','Color','k');% add annotation
annnum = annotation('textbox',[0.01 0.01 0.2 0.1],'interpreter','latex',...
    'String',{['(a)']}, 'LineStyle','-','EdgeColor','None',...
    'BackgroundColor','none','Color','k','VerticalAlignment','bottom');% add numbering annotation
% title('Heat capacity of TmVO$_4$')
grid on;%
hold off

%% Compute additional Cp for Gaussian distribution of fields
% Cptheo(i).h = 0.87;
Cpnum(i).rhsgm = 0.1;
Cpnum(i).sgm = Cpnum(i).h*Cpnum(i).rhsgm;
for j=1:length(Cpnum(i).t)
    Cpnum(i).normpdf(j) = CpTFIM_normpdf(Cpnum(i).t(j),Cpnum(i).h,Cpnum(i).sgm);
end
plot(Cpnum(i).t*Tc0,Cpnum(i).normpdf,'DisplayName',sprintf('h=%.2f,r=%.1e',Cpnum(i).h,Cpnum(i).rhsgm));

%% Plot additional Cp for Gaussian distribution of fields
offset = 0.01;
plot(Cpnum(i).t*Tc0,Cpnum(i).normpdf+offset,'DisplayName',sprintf('h=%.2f,r=%.1e',Cpnum(i).h,Cpnum(i).rhsgm));

%% Prepare MF fit of Cp vs Temperature under field
index = 15;
H1 = uhrf(index);
T1 = avgRFData(index).T;
Cp1 = avgRFData(index).Cp/R;
Cp1Err = avgRFData(index).CpFullErr/R;
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

%% Export figure
% formatFigure;
printPNG([todaystr '_TmVO4-RF-E_mfd@4000Oe']);
% printPDF(['2019-06-18_TmVO4-RF-E_fit_Schottky_' strrep(hrstr,'.','p') 'xHc']);



%% 3D scatter of Cp(T) at each field separately
figure
for i=1:length(uhrf)
    scatter3(separatedRFCpData(i).H,separatedRFCpData(i).T,separatedRFCpData(i).Cp,'.',...
        'Displayname',num2str(uhrf(i)', '%d Oe'))
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
for i=1:length(uhrf)
    scatter3(separatedRFCpData(i).H,separatedRFCpData(i).T,-separatedRFCpData(i).d1Cp,'.',...
        'Displayname',num2str(uhrf(i)', '%d Oe'))
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
for i=1:length(uhrf)
    plot(separatedRFCpData(i).T,separatedRFCpData(i).Cp,'-+')
    hold on
end
xlabel('Temperature (K)')
ylabel('Cp/T (J/K^2/mol)')
legendCell = cellstr(num2str(uhrf, '%-d Oe'));
legend(legendCell)
hold off
