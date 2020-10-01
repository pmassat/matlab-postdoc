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

%% Identify experimental critical temperature at each field
% Then we can plot Tc vs uh and fit using equation f(h) = Tc/Hc*h/atanh(h/Hc)
% From this, we can see that the experimental value of Hc is ~0.72T 
% and hence correct for demag, knowing the value of critical field Hc0 (see beginning of code)
M = ones(1,length(uhDy2));
Tcd1 = ones(1,length(uhDy2));
for i=1:length(uhDy2)
    [M(i),I] = min(avgCpDy2Data(i).d1Cp);
    Tcd1(i) = avgCpDy2Data(i).T(I);
end

%% Prepare MF fit of Cp vs Temperature at given field
j=1;
Tfit = avgCpDy2Data(j).T;
Cpfit = avgCpDy2Data(j).Cp;
CpfitErr = avgCpDy2Data(j).CpFullErr;
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
fitT = avgCpDy2Data(i).T;
fitCp = avgCpDy2Data(i).Cpelr;
fitCpErr = avgCpDy2Data(i).CpFullErr/R;
fitwghts = 1./fitCpErr;
% Tc = 2.15;

%% Compute Cp for distribution of fields
clear CpnumRF
rngNum = 1:length(uhDy2);
for i=rngNum
CpnumRF(i).h = uhDy2(i)*rescaling/Hc;
% CpnumRF(i).rhsgm = 0.09;
% CpnumRF(i).sgm = CpnumRF(i).h*CpnumRF(i).rhsgm;
CpnumRF(i).t_single_h_no_e = linspace(0,1.5,601);% reduced temperature, T/Tc
CpnumRF(i).single_h_no_e = zeros(size(CpnumRF(i).t_single_h));
CpnumRF(i).t_single_h_w_e = CpnumRF(i).t_single_h_no_e(2:end-1);%
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
rng = find(ismember(uhDy2,10^3*[1:8]))';%[1 5 7 9 13];
rngPlot = rng(3:end-1);

plot_Cp_avg_w_fits(rngPlot, avgCpDy2Data, CpnumRF, Tc0rf, uhDy2/Hc)

% Fit parameters on data at H=0: Tc=2.125(3), e=1.5(4)e-3
% Note: the values of amplitude coefficient and Tc extracted from fit 
% in curve fitting tool using Cp_TFIM (no offset strain) are A=7.35 and Tc=2.142K

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
uhmfdrf = unique(Tmfd_RF.Hext_Oe);

%% Outdated; For a given dataset, find closest values of temperature and field in COMSOL mfd
i = 5;
Hdata = unique(round(avgCpDy2Data(i).H,-2));
[~,mfdhidx] = min(abs(Tmfd_RF.Hext_Oe-Hdata));
h = Tmfd_RF.Hext_Oe(mfdhidx);
tref = 0;
for j=1:length(CpnumRF(i).t_h_dist_noe)
    % Find value of temperature in COMSOL mfd closest to that of actual data
    [~,mfdtidx] = min(abs(Tmfd_RF.T_K/Tc0rf-CpnumRF(i).t_h_dist_noe(j)));
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
rngMFD = find(ismember(uhDy2,10^3*[1:8]));

for idx=length(rngMFD);%rngNum(3:end)
% idx = 8
    i = rngMFD(idx)
    
    Hext_data = unique(round(avgCpDy2Data(i).H,-1));
    [~,mfdhidx] = min(abs(Tmfd_RF.Hext_Oe-Hext_data));
    Hext_mfd = Tmfd_RF.Hext_Oe(mfdhidx)
    utmfdrf = unique(Tmfd_RF.T_K(Tmfd_RF.Hext_Oe==Hext_mfd));
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

%     toc
    
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
    spikes = abs(Cpwewt)>1.2*max(avgCpDy2Data(i).Cpelr) | Cpwewt<0;% 
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

end

toc 

%% Plot Cp for COMSOL distribution of fields

for idx=1%:length(rngMFD)
    i = rngMFD(idx);
    single_str = ['$H_{\mathrm{ext}}\times$' sprintf('%.2g',rescaling)];
    figure
    % plot(avgRFData(i).T,avgRFData(i).Cpelr,'.','DisplayName','data')
    eb{rngNum==i} = errorbar(avgCpDy2Data(i).T,avgCpDy2Data(i).Cpelr,avgCpDy2Data(i).CpelrErr,...
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
    title(['$C_p$ single field vs COMSOL PDF $H_{\mathrm{ext}}=$ ' sprintf('%.0f Oe',uhDy2(i))]);
    lgd = legend('Location','northwest');% title(lgd,'TmVO4-Ndl-E');
    lgd.FontSize = 14;

    xlabel('$T$ (K)')
    ylabel('$C_p/R$')
    
    xlim([0 3.2])
end

%% Prepare exportation of figure/data
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2017-07_TmVO4_Cp_MCE\2017-07-20_Cp\2017-07-20_TmVO4_Cp_analysis'
mfd_ID = '_TmVO4_half-diamond_V=Vsample_Hzm=0p71513';

%% Create file name for field i
for i=rng
fname_field{i} = sprintf('%s_Cpnum_Hext=%dkOe', mfd_ID, uhDy2(i)/10^3);
export_dat = [todaystr fname_field{i} '.dat'];
import_dat = ['2020-09-09' fname_field{i} '.dat'];

%% Import results from .dat file
fileID = fopen(import_dat,'r');
A = textscan(fileID, '%f%f', 'delimiter', ',', 'HeaderLines', 7);
CpnumRF(i).t_h_dist_w_e = A{1}';
CpnumRF(i).comsolpdf_w_e = A{2}';
fclose(fileID);

end

%% Plot data and computed Cp
plot_Cp_avg_w_fits(rngPlot, avgCpDy2Data, CpnumRF, Tc0rf, uhDy2/Hc,...
    'TnumStr', 't_h_dist_w_e', 'CpnumStr', 'comsolpdf_w_e')

%% Export figure
% Place cursor on this line to run loop containing subsections
% for idx=1:4%length(rngMFD)
% i=rngMFD(idx)

% formatFigure;
% printPNG([todaystr mfd_ID '_Cp_vs_T_H'])
% printPDF(['2019-06-18_TmVO4-RF-E_fit_Schottky_' strrep(hrstr,'.','p') 'xHc']);
% end

%% Export results to .dat file
% Note: in order to run a for loop that contains subsections like this
% one, the cursor has to be outside the loop when executing the run command
% save_Cpnum_data(CpnumRF(i), uhrf(i)/10^3, mfd_ID, e, Tc0rf, dat_name);
    












%% Plot Cp for Gaussian distribution of fields
figure
plot(avgCpDy2Data(i).T,avgCpDy2Data(i).Cpelr,'.','DisplayName','data')
hold on;
plot(CpnumRF(i).t_single_h*Tc0rf,CpnumRF(i).single_h_no_e,'DisplayName','MF');
plot(CpnumRF(i).t_h_dist_noe*Tc0rf,CpnumRF(i).normpdf,'DisplayName','Gauss');
title(['Cp mean-field vs normal pdf $H_{\mathrm{ext}}=$' sprintf('%.0fOe',uhDy2(i))]);
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
eb{rng==i} = errorbar(avgCpDy2Data(i).T,avgCpDy2Data(i).Cpelr,avgCpDy2Data(i).CpFullErr/R,...
    '.','MarkerSize',18,'DisplayName',num2str(uhDy2(i)/(Hc0*1e4),'%.2f'),...
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
H1 = uhDy2(index);
T1 = avgCpDy2Data(index).T;
Cp1 = avgCpDy2Data(index).Cp/R;
Cp1Err = avgCpDy2Data(index).CpFullErr/R;
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
for i=1:length(uhDy2)
    scatter3(separatedDy2CpData(i).H,separatedDy2CpData(i).T,separatedDy2CpData(i).Cp,'.',...
        'Displayname',num2str(uhDy2(i)', '%d Oe'))
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
for i=1:length(uhDy2)
    scatter3(separatedDy2CpData(i).H,separatedDy2CpData(i).T,-separatedDy2CpData(i).d1Cp,'.',...
        'Displayname',num2str(uhDy2(i)', '%d Oe'))
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
for i=1:length(uhDy2)
    plot(separatedDy2CpData(i).T,separatedDy2CpData(i).Cp,'-+')
    hold on
end
xlabel('Temperature (K)')
ylabel('Cp/T (J/K^2/mol)')
legendCell = cellstr(num2str(uhDy2, '%-d Oe'));
legend(legendCell)
hold off
