%% Sample properties
m = 0.25e-3;% mass of sample, in g
M = 283.87;% molar mass of TmVO4, in g/mol
Tc0rf = 2.126;% Value of transition temperature in this sample, in Kelvin units
Tcnum = 2.15;% in Kelvin units
Hc = 5000;% in Oersted units; see data taken on needles of TmVO4-LS5200 in July 2017
e = 1.5e-3;% constant longitudinal field
% Results for Tmaxfit = 2.15;% Tc = 2.126  (2.122, 2.129);% e = 0.001474 (0.001085, 0.001862);

%% Import magnetic field distribution from TXT file
% cd 'C:\Users\Pierre\Desktop\Postdoc\Software\COMSOL\TmVO4-RF-E_HC2017-07'
% mfdRFfname = '2020-07-21_TmVO4-RF-E_zero-temp-magnetization_mag-field_distrib.txt';
% Smfd = importfielddistrib(mfdRFfname, 15);

%% Import magnetic field distribution from CSV file
% cd 'C:\Users\Pierre\Desktop\Postdoc\Software\COMSOL\TmVO4-RF-E_HC2017-07\TmVO4-RF-E_magnetic_field_distribution_data'
cd 'C:\Users\Pierre\Desktop\Postdoc\Software\COMSOL\TmVO4_basic_shapes\TmVO4_basic_shapes_MFD_data'
% mfdRFfname = cell(1,2);
% Srf = cell(1,2);
% % mfdRFf{1} = {'2020-07-27_TmVO4-RF-E_COMSOL_mfd_mesh=finer-out+max20um-in_T=p001-1K_H=all.csv',28};
% mfdRFf{2} = {'2020-07-27_TmVO4-RF-E_COMSOL_mfd_mesh=finer-out+max20um-in_T=2-3K_H=all.csv',28};
% mfdRFf{3} = {'2020-08-11_TmVO4-RF-E_COMSOL_mfd_mesh=max20um-in_T=p5-1-1p5K_H=p1-p1-p8.csv',24};
% mfdRFf{1} = {'2020-08-25_TmVO4-RF-E_COMSOL_mfd_mesh=max20um-in_T=p3-p4-3p1K_H=p1-p1-p8T.csv',64};
Ntemp = 6;
Nfields = 1;
mfdRF{1} = {'2020-09-08_TmVO4_half-diamond+front-cut_V=Vsample_Hzm=0p71513_mesh=max20um-in_comsol-mfd_T=p3-3K-6values_H=5kOe.csv',Ntemp*Nfields};
for ic=length(mfdRF):-1:1
%     Srf{ic} = importfielddistrib_csv(mfdRF{ic}{1}, mfdRF{ic}{2}, 'compNum', 2);
    Srf{ic} = importfielddistrib_csv(mfdRF{ic}{1}, mfdRF{ic}{2});
end

%% Extract values of temperature and external magnetic field from structure header
Smfd_RF = [Srf{1,:}];% concatenate cell arrays into structure
for sidx=1:length(Smfd_RF)
    % split each header into a cell array of strings, using whitespace and '=' sign as separators
    TBextCell = strsplit(Smfd_RF(sidx).T_Bext,{' ','='});
    
    % Find the index of the string cell containing "T" (temperature), and
    % convert the following string cell (which contains the value of temperature) into a number
    % Note: circshift(A,1) shifts the index of elements of array A by 1 to the right
    try
        Smfd_RF(sidx).T_K = str2double(TBextCell{circshift(TBextCell=="T",1)});
    catch
        user_temp = 1;% K
        Smfd_RF(sidx).T_K = user_temp;
        warning(sprintf('No temperature values found. Using user value: %.2g K',user_temp))
    end
    
    % Same with "Bext", which is the external magnetic flux density, in Tesla
    % units, and convert it into a value of magnetic field, in Oersted units
    try
        Smfd_RF(sidx).Hext_Oe = str2double(TBextCell{circshift(TBextCell=="Bext",1)})*10^4;
    catch
        user_mag = 5000;% Oe
        Smfd_RF(sidx).Hext_Oe = user_mag;
        warning(sprintf('No magnetic field values found. Using user value: %i Oe',user_mag))
    end
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
% Ntemp = 7;
param_index = 1;% 1 is constant T, 2 is constant Hext, see param_range
temp_index = 1;% determines temperature of data to plot, taken among [.3:.4:3.1]K
field_index = 1;% determines field of data to plot from [.1,.2,.3,.4,.45,.5,.55,.6,.62,.63,.65,.67,.7,.8]T
param_range = {(temp_index-1)*Ntemp+[1:1:Ntemp], field_index+[0:Ntemp:56]};% first range corresponds to a 
% field dependence at constant temp, second range corresponds to a 
% temperature dependence at constant field
rng = param_range{param_index};
for i=rng
T = Smfd_RF(i).T_K;%
Hext = Smfd_RF(i).Hext_Oe;%
p = plot(Smfd_RF(i).binCenters, Smfd_RF(i).hc, '.-', 'DisplayName', sprintf('%.2g, %.2g',T/Tcnum,Hext/Hc));
end
lgd = legend('show'); lgd.Title.String = '$T/T_{c,0}$, $H_{\mathrm{ext}}/H_{c,0}$';
param_title = {[', $T=$ ' sprintf('%.2g K',T)],...
    [', $H_{\mathrm{ext}}=$ ' sprintf('%.2d Oe',Hext)]};
title(['Distribution of fields in TmVO$_4$-RF-E' param_title{param_index}])
xlabel('$H_{\mathrm{in}}/H_{\mathrm{c}}$')
ylabel('Normalized PDF')

%% Compute the average value of ratio of internal to external magnetic field 
for i=1:length(Smfd_RF)
Smfd_RF(i).Hinm_Oe = round(...
    sum(...
    Smfd_RF(i).binCenters.*Smfd_RF(i).hc.*Smfd_RF(i).binWidths...
    )*Hc...
    );
end

%% Check results for given range of temperatures/field
% idx = 1; rng = (idx-1)*8+[1:8]; [Smfd_RF(rng).T_K;Smfd_RF(rng).Hinm_Oe]

% % Check that the ratio of Hinm/Hext is roughly constant inside the
% % ordered phase, thus allowing to use the value at 0.3 K and 1000 Oe as a proxy
% for i=rng
%     sprintf('%d',Smfd_RF(i).Hinm/Smfd_RF(i).Hext_Oe)
% end
rescaling = Smfd_RF(1).Hinm_Oe/Smfd_RF(1).Hext_Oe;%Hc0/0.69;% rescaling factor, due to demag; 

%% Compute the relative difference between curves at 1000 Oe and 7000 Oe
mfd_diff(113).abs_diff = abs(Smfd_RF(1).hc);%-Smfd(13).hc);
% mfd_diff(113).rel = mfd_diff(113).abs/Smfd(1).binWidths(1);

%% Data importation
sample = 'TmVO4-RF-E';
filename = 'TmVO4_RF-E_2017-07-14.dat';
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2017-07_TmVO4_Cp_MCE\2017-07-20_Cp\2017-07-20_TmVO4_Cp_analysis'
RFDATA=ImportTmVO4Cp(filename);% Use this data to plot color map of phase diagram

%% Assign data to variables
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
contourf(Hgm./5100,Tgm/Tc0rf,-d1Cpgm,n,'EdgeColor','none');
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
Tgmr = Tgm(:,1)/Tc0rf;
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

%% Prepare fit of Cp vs Temperature 
i = 9;
% Use these variables in Curve fitting tool
clear fitT fitCp fitCpErr fitwghts 
fitT = avgRFData(i).T;
fitCp = avgRFData(i).Cpelr;
fitCpErr = avgRFData(i).CpFullErr/R;
fitwghts = 1./fitCpErr;
% Tc = 2.15;

%% Compute Cp for distribution of fields
clear CpnumRF
rngNum = 7;
for i=rngNum
CpnumRF(i).h = uhrf(i)*rescaling/Hc;
% CpnumRF(i).rhsgm = 0.09;
% CpnumRF(i).sgm = CpnumRF(i).h*CpnumRF(i).rhsgm;
CpnumRF(i).t_single_h = linspace(0,1.5,601);% reduced temperature, T/Tc
CpnumRF(i).single_h_no_e = zeros(size(CpnumRF(i).t_single_h));
CpnumRF(i).single_h_w_e = zeros(size(CpnumRF(i).t_single_h));
CpnumRF(i).t_h_dist_noe = [linspace(0.01,.5,50) linspace(0.505,1.1,120) linspace(1.11,1.4,30)];% reduced temperature, T/Tc
CpnumRF(i).normpdf = zeros(size(CpnumRF(i).t_h_dist_noe));
CpnumRF(i).comsolpdf_no_e = zeros(size(CpnumRF(i).t_h_dist_noe));
CpnumRF(i).comsolpdf_w_e = zeros(size(CpnumRF(i).t_h_dist_noe));
end
% Fit parameters on data at H=0: Tc=2.125(3), e=1.5(4)e-3
%     for j=2:length(CpnumRF(i).t_h_dist_noe)
%     CpnumRF(i).normpdf(j) = CpTFIM_normpdf(CpnumRF(i).t_h_dist_noe(j),CpnumRF(i).h,CpnumRF(i).sgm);
%     end

%% Compute Cp_TLFIM at single value of field
tic
for i=rngNum%(2:end)
CpnumRF(i).single_h_no_e = Cp_TFIM(CpnumRF(i).t_single_h,CpnumRF(i).h);
[~,~,CpnumRF(i).single_h_w_e] = FSCp_TLFIM(CpnumRF(i).t_single_h,CpnumRF(i).h,e);
CpnumRF(i).single_h_w_e = CpnumRF(i).single_h_w_e';
end
toc 

%% Plot averaged data at each field separately
figure; hold on
rng = rngNum;% [1 5 7 9 13];
clr = lines(length(rng));
eb = cell(size(rng));
maxTplot = 3.2;%
for i=rng
    fp = plot(CpnumRF(i).t_single_h(2:end-1)*Tc0rf,CpnumRF(i).single_h_w_e,'DisplayName','MF');
%     fp = fplot(@(t) Cp_TFIM(t/Tc0rf,uhrf(i)*rescaling/(Hc0*1e4)),[0 maxTplot],'LineWidth',2);
% Fit parameters on data at H=0: Tc=2.125(3), e=1.5(4)e-3
% Note: the values of amplitude coefficient and Tc extracted from fit 
% in curve fitting tool using Cp_TFIM (no offset strain) are A=7.35 and Tc=2.142K
%     clr{rng==i} = get(fp,'Color');
end
for i=rng
    eb{rng==i} = errorbar(avgRFData(i).T,avgRFData(i).Cpelr,avgRFData(i).CpFullErr/R,...
        '.','MarkerSize',18,'DisplayName',num2str(uhrf(i)/Hc,'%.2f'),...
        'Color',clr(rng==i,:),'LineWidth',2);
end
xlabel('Temperature (K)'); ylabel('C$_p$/R');%ylabel('C$_p$ (JK$^{-1}$mol$^{-1}$)');
xlim([0 maxTplot])
title('Heat capacity of TmVO4-RF-E')
lgd = legend([eb{:}]); lgd.Title.String = '$H_{\mathrm{ext}}/H_c$';
ax = gca; ax.YMinorTick = 'on';% Add minor ticks on Y axis
grid on;%
hold off

%% Gaussian function (see normpdf.m)
mu = .58; 
sigma =.09*mu; 
gauss = @(x) exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);

%% Compute probability distribution of fields at a given value of T and Hext
Hcnum = Hc;
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
utmfdrf = unique(Tmfd_RF.T_K);
uhmfdrf = unique(Tmfd_RF.Hext_Oe);

%% Outdated; For a given dataset, find closest values of temperature and field in COMSOL mfd
i = 5;
Hdata = unique(round(avgRFData(i).H,-2));
[~,mfdhidx] = min(abs(Tmfd_RF.Hext_Oe-Hdata));
h = Tmfd_RF.Hext_Oe(mfdhidx);
tref = 0;
for j=1:length(CpnumRF(i).t_h_dist_noe)
    % Find value of temperature in COMSOL mfd closest to that of actual data
    [~,mfdtidx] = min(abs(Tmfd_RF.T_K/Tcnum-CpnumRF(i).t_h_dist_noe(j)));
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
        Cph(col) = Cp_TFIM(CpnumRF(i).t_h_dist_noe(j),Tmfd_RF.binCenters(row,col));
    end
    % Compute the corresponding value of heat capacity 
%     CpnumRF(i).comsolpdf(j) = trapz(Tmfd_RF.binCenters(row,:),Cph.*Tmfd_RF.hc(row,:));
    CpnumRF(i).comsolpdf(j) = sum(Tmfd_ndl.binWidths(row,:).*Cph.*...
        Tmfd_ndl.hc(row,:));
end

%% For a given dataset, compute Cp with and without longitudinal field
% for the two closest values of temperature and field in COMSOL mfd
tic% start clock to measure computation time

i=7;%rngNum(3:end)
    
    Hext_data = unique(round(avgRFData(i).H,-1));
    [~,mfdhidx] = min(abs(Tmfd_RF.Hext_Oe-Hext_data));
    Hext_mfd = Tmfd_RF.Hext_Oe(mfdhidx);
    trange = CpnumRF(i).t_h_dist_noe(1:end);
    tref = zeros(2,length(trange)+1);% Not even necessary given the the function find_mfd_temp_rows() includes [0 0] as default values for tref
    % For the computations with longitudinal field, need to compute free
    % energy first, since there is no analytical formula for the heat capacity
    Fwewt = zeros(2,length(trange));% Initialize array of weighted partial free energies
    Cpnoewt = zeros(2,length(trange));% Initialize array of weighted partial heat capacities (w/o longitudinal strain)
    % t = zeros(2,length(CpnumRF(i).t_h_dist_noe));% array that will contain two closest values of temperature
    wt = zeros(2,length(trange));% array that will contain weights attributed to closest temperatures
    rf_rows = zeros(2,length(trange));
    Cphnoe = zeros(size(Tmfd_RF.binCenters(1,:)));

    for jt=1:length(trange)
    %     [rows, wt(:,jt)] = find_mfd_temp_rows(trange(jt),...
    %         utmfd1, Tc0rf, Tmfd_RF, Hext_mfd );
        % add tref as input and output of function find_mfd_temp_rows() to check 
        % that the function finds the correct temperatures to compute CpnumRF
        [rows, wt(:,jt), tref(:,jt+1)] = find_mfd_temp_rows(trange(jt),...
            utmfdrf, Tc0rf, Tmfd_RF, Hext_mfd, tref(:,jt), 'printTref', true );

        rf_rows(:,jt) = rows;

        % Compute heat capacity at each value of internal (transverse) 
        % magnetic field for both nearest temperatures
        for jr=1:length(rows)

            % Compute free energy, including longitudinal field
            % Note: it is important to compute this array line by line using
            % the for loop over jr, because the computation over the entire
            % array at once induces computation errors, which are hard to
            % spot at first, because they only arise at the 5th significant
            % digit, but this is enough to make a big difference in the computation
            % of the second derivative, i.e. the heat capacity
            Fhwe(jr,:) = FSCp_TLFIM(trange(jt),Tmfd_RF.binCenters(rows(jr),:),e);

            % Compute Cp without longitudinal field, for each value of internal field
            for col=1:length(Tmfd_RF.binCenters(rows,:))
                Cphnoe(jr,col) = Cp_TFIM(trange(jt),...
                    Tmfd_RF.binCenters(rows(jr),col));
            end

            % Sum over all internal fields to get the total heat capacity at
            % the temperature associated with jt
    %        Cpnoe(jr,jt) = wt(jr,jt)*sum(...
            Cpnoewt(jr,jt) = wt(jr,jt).*sum(...
                Tmfd_RF.binWidths(rows(jr),:).*...
                Cphnoe(jr,:).*...
                Tmfd_RF.hc(rows(jr),:)...
                );
        end

        % Sum the free energy over all internal fields to get the total free energy at
        % the temperature associated with jt
    %     Fwewt(:,jt) = sum(...
        Fwewt(:,jt) = sum(...% Use this if trange does not cover the full range of CpnumRF(i).t_h_dist_noe
            Tmfd_RF.binWidths(rows,:).*...
            Fhwe.*...
            Tmfd_RF.hc(rows,:),...
            2);

    end

    toc
    
    %%
    % CpnumRF(i).t_h_dist_noe = linspace(0.1,1.5,150);% reduced temperature, T/Tc
    Cpnoe = sum(Cpnoewt);

    %% Compute entropy and heat capacity, with longitudinal fields, at each value of internal field
    twe = trange;
    dt = diff(twe);
    Swewt = -diff(Fwewt,1,2)./dt;% entropy
    dtm = 0.5*(dt(1:end-1)+dt(2:end));
    Cpwewt = twe(2:end-1).*diff(Swewt,1,2)./dtm;

    %% Identify and remove spikes
    % Spikes result from numerical aberrations, most likely due to
    % significant changes in the MFD at "high" temperatures, i.e. close to
    % and above Tc(H=0)
    
    % Identify spikes
    spikes = abs(Cpwewt)>1.2*max(avgRFData(i).Cpelr) | Cpwewt<0;% 
%     Cpwewt(abs(Cpwewt)>1.2*max(avgRFData(i).Cpelr)) = NaN;% Remove numerical aberrations
%     Cpwewt(Cpwewt<0) = NaN;% Remove numerical aberrations

    % Convert spikes appearing in both computed Cp curves into NaN for both curves
    Cpwewt(:,spikes(1,:)&spikes(2,:)) = NaN;
    
    % For those spikes appearing in Cp of second closest temperature
    % only, replace them with the values of the closest temperature
    Cpwewt(2,~spikes(1,:)&spikes(2,:)) = Cpwewt(1,~spikes(1,:)&spikes(2,:));

    %% Compute free energy as sum of weighted partial free energies
    Cpwe = sum(wt(:,2:end-1).*Cpwewt);

    %% Store results into CpnumRF structure
    CpnumRF(i).t_h_dist_w_e = twe(2:end-1);
    CpnumRF(i).comsolpdf_no_e = Cpnoe;
    CpnumRF(i).comsolpdf_w_e = Cpwe;

    %% Prepare data for plotting by removing NaN datapoints
    sel_no_e = ~isnan(CpnumRF(i).comsolpdf_no_e);
    sel_w_e = ~isnan(CpnumRF(i).comsolpdf_w_e);
    CpnumRF(i).t_h_dist_noe = CpnumRF(i).t_h_dist_noe(sel_no_e);
    CpnumRF(i).t_h_dist_w_e = CpnumRF(i).t_h_dist_w_e(sel_w_e);
    CpnumRF(i).comsolpdf_no_e = CpnumRF(i).comsolpdf_no_e(sel_no_e);
    CpnumRF(i).comsolpdf_w_e  = CpnumRF(i).comsolpdf_w_e(sel_w_e);

% end

% toc 

%% Plot Cp for COMSOL distribution of fields
for i=7
single_str = ['$H_{\mathrm{ext}}\times$' sprintf('%.2g',rescaling)];
figure
% plot(avgRFData(i).T,avgRFData(i).Cpelr,'.','DisplayName','data')
eb{rngNum==i} = errorbar(avgRFData(i).T,avgRFData(i).Cpelr,avgRFData(i).CpelrErr,...
    '.','MarkerSize',18,'DisplayName',['Data at $H=H_{\mathrm{ext}}$'],...
    'LineWidth',2);
hold on;
no_e_str = ' ($e=0$)';
w_e_str = [' ($e=$ ' sprintf(' %.2g)',e)];

% Plot results at single value of magnetic field
plot(CpnumRF(i).t_single_h*Tc0rf,CpnumRF(i).single_h_no_e,...
    'DisplayName',[single_str no_e_str]);% no longitudinal field
plot(CpnumRF(i).t_single_h(2:end-1)*Tc0rf,CpnumRF(i).single_h_w_e,...
    'DisplayName',[single_str w_e_str]);% with longitudinal field

% Plot results for distribution of magnetic fields as computed with COMSOL
comsol_str = 'comsol pdf';
plot(CpnumRF(i).t_h_dist_noe*Tc0rf,CpnumRF(i).comsolpdf_no_e,...
    'DisplayName',[comsol_str no_e_str]);
plot(CpnumRF(i).t_h_dist_w_e*Tc0rf,CpnumRF(i).comsolpdf_w_e,...
    'DisplayName',[comsol_str w_e_str]);
title(['$C_p$ single field vs COMSOL PDF $H_{\mathrm{ext}}=$ ' sprintf('%.0f Oe',uhrf(i))]);
lgd = legend('Location','northeast');% title(lgd,'TmVO4-Ndl-E');
lgd.FontSize = 14;

xlabel('$T$ (K)')
ylabel('$C_p/R$')
end

%% Export figure
% cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2017-07_TmVO4_Cp_MCE\2017-07-20_Cp\2017-07-20_TmVO4_Cp_analysis'
formatFigure;
% printPNG([todaystr '_TmVO4_half-diamond+front-cut_Cp_fits_Hext=5kOe_Comsol_V=Vsample_Hzm=p71513']);
% printPDF(['2019-06-18_TmVO4-RF-E_fit_Schottky_' strrep(hrstr,'.','p') 'xHc']);










%% Plot Cp for Gaussian distribution of fields
figure
plot(avgRFData(i).T,avgRFData(i).Cpelr,'.','DisplayName','data')
hold on;
plot(CpnumRF(i).t_single_h*Tc0rf,CpnumRF(i).single_h_no_e,'DisplayName','MF');
plot(CpnumRF(i).t_h_dist_noe*Tc0rf,CpnumRF(i).normpdf,'DisplayName','Gauss');
title(['Cp mean-field vs normal pdf $H_{\mathrm{ext}}=$' sprintf('%.0fOe',uhrf(i))]);
lgd = legend();% title(lgd,'TmVO4-RF-E');

%% Plot averaged data at each field separately
figure; hold on
clr = lines(length(rng));
eb = cell(size(rng));
for i=rng
% fp = fplot(@(t)Cp_TFIM(t/Tc0rf,Cptheo(i).h),[0 3.2],'--','LineWidth',2,'Color',clr(rng==i,:));
plot(CpnumRF(i).t_single_h*Tc0rf,CpnumRF(i).single_h_no_e,'--','Color',clr(rng==i,:),'DisplayName',sprintf('h=%.2f',CpnumRF(i).h));
plot(CpnumRF(i).t_h_dist_noe*Tc0rf,CpnumRF(i).normpdf,'Color',clr(rng==i,:),'DisplayName',...
    sprintf('h=%.2f,r=%.1e',CpnumRF(i).h,CpnumRF(i).rhsgm));
end
for i=rng
eb{rng==i} = errorbar(avgRFData(i).T,avgRFData(i).Cpelr,avgRFData(i).CpFullErr/R,...
    '.','MarkerSize',18,'DisplayName',num2str(uhrf(i)/(Hc0*1e4),'%.2f'),...
    'Color',clr(rng==i,:),'LineWidth',2);
end
xlabel('$T$ (K)'); ylabel('$C_p/R$');%ylabel('C$_p$ (JK$^{-1}$mol$^{-1}$)');
xlim([0 max(CpnumRF(rng(1)).t_single_h*Tc0rf)]);
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
CpnumRF(i).rhsgm = 0.1;
CpnumRF(i).sgm = CpnumRF(i).h*CpnumRF(i).rhsgm;
for j=1:length(CpnumRF(i).t)
    CpnumRF(i).normpdf(j) = CpTFIM_normpdf(CpnumRF(i).t(j),CpnumRF(i).h,CpnumRF(i).sgm);
end
plot(CpnumRF(i).t*Tc0rf,CpnumRF(i).normpdf,'DisplayName',sprintf('h=%.2f,r=%.1e',CpnumRF(i).h,CpnumRF(i).rhsgm));

%% Plot additional Cp for Gaussian distribution of fields
offset = 0.01;
plot(CpnumRF(i).t*Tc0rf,CpnumRF(i).normpdf+offset,'DisplayName',sprintf('h=%.2f,r=%.1e',CpnumRF(i).h,CpnumRF(i).rhsgm));

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
