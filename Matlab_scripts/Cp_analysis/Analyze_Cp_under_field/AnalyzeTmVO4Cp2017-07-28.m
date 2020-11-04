Tc_ndl = 2.22;% transition temperature at zero field, in Kelvin units
Hc = 5000;% critical field at zero temperature, in Oersted units
cd 'C:\Users\Pierre\Desktop\Postdoc\Software\COMSOL\TmVO4-LS5200_HC2017-07\TmVO4-LS5200_HC2017-07_COMSOL_results'

%% Import magnetic field distribution from CSV file
clear mfdNdlFname Sndl
mfdNdlFname{1} = '2020-11-03_TmVO4-LS5200_HC17-VII_T=p3-p4-3p1_Hext=all.csv';

for ic=1:length(mfdNdlFname)
    Sndl{ic} = importfielddistrib_csv(mfdNdlFname{ic}, 56);
end

%% Extract values of temperature and external magnetic field from structure header
Smfd_ndl = [Sndl{1,:}];% concatenate structures
for sidx=1:length(Smfd_ndl)
    % split each header into a cell array of strings, using whitespace and '=' sign as separators
    TBextCell = strsplit(Smfd_ndl(sidx).T_Bext,{' ','='});
    % Find the index of the string cell containing "T" (temperature), and
    % convert the following string cell (which contains the value of temperature) into a number
    % Note: circshift(A,1) shifts the index of elements of array A by 1 to the right
    Smfd_ndl(sidx).T_K = str2double(TBextCell{circshift(TBextCell=="T",1)});
    % Same with "Bext", which is the external magnetic flux density, in Tesla
    % units, and convert it into a value of magnetic field, in Oersted units
    Smfd_ndl(sidx).Hext_Oe = str2double(TBextCell{circshift(TBextCell=="Bext",1)})*10^4;
%     if sidx<=56
%         Smfd_ndl(sidx).label = 'needle1';
%     else
%         Smfd_ndl(sidx).label = 'needle2';
%     end
end

%% Compute probability distribution of fields at a given value of T and Hext
for i=1:length(Smfd_ndl)
mfd = Smfd_ndl(i).mfd(Smfd_ndl(i).mfd>0)/Smfd_ndl(i).Hext_Oe;% create distribution from non-zero values; also possible to normalize wrt Hc
% h = histogram(mfd, 'Normalization', 'pdf');% plot histogram of distribution
[Smfd_ndl(i).hc, edges] = histcounts(mfd, 50, 'Normalization', 'pdf');% plot histogram of distribution
% binCenters = h.BinEdges + (h.BinWidth/2);
Smfd_ndl(i).binCenters = mean([edges(1:end-1);edges(2:end)],1);
Smfd_ndl(i).binWidths = edges(2:end)-edges(1:end-1);
end

%% Plot distribution of fields at a given value of T and Hext
figure; hold on
% needle_index = 1;% there are 2 needle-shaped samples 
% start_index = [0,56];% needle 1 data start at index 0+1, needle 2 data at 56+1
param_index = 2;% 1 is constant T, 2 is constant Hext, see param_range
param_range = {[1:8*2:56], 4*8+[1,4,6,8]};% first range corresponds to a 
% field dependence at constant temp, second range corresponds to a 
% temperature dependence at constant field
rng = param_range{param_index};% + start_index(needle_index);

for i=rng 
    Tndl = Smfd_ndl(i).T_K;%
    Hext = Smfd_ndl(i).Hext_Oe;%
    lgd_str = legend_string(param_index, Hext/Hc, Tndl/Tc_ndl)
    p = plot(Smfd_ndl(i).binCenters, Smfd_ndl(i).hc, '.-', 'DisplayName', lgd_str);
end
xlim([.6966 1.05])

lgd = legend('show','Location','northwest'); 
if param_index==1
    lgd.Title.String = '$H_{\mathrm{ext}}/H_c$';
    ann_str = ['$T/T_c=$' sprintf(' %.2g', Tndl/Tc_ndl)];
elseif param_index==2
    lgd.Title.String = '$T/T_c$';
    ann_str = ['$H_{\mathrm{ext}} =$' sprintf(' %.3g kOe', Hext/1e3)];
end
anndist = annotation('textbox',[0.35 0.8 0.2 0.1], 'interpreter','latex',...
    'String', ann_str, 'EdgeColor','none', 'FitBoxToText', 'on',...
    'BackgroundColor','none', 'Color','k');% add annotation
% needle_title = sprintf('Field distribution in TmVO$_4$ needles');
% param_title = {[', $T=$' sprintf('%.2g K',Tndl)],...
%     [', $H_{\mathrm{ext}}=$ ' sprintf('%.2d Oe',Hext)]};
% title([needle_title param_title{param_index}])
xlabel('$H_{\mathrm{in}}/H_{\mathrm{ext}}$')
ylabel('Probability density')
ax = gca; ax.XTick = .7:.1:1;
grid on

%% Export figure
% formatFigure
% printPNG([todaystr '_TmVO4-needles-2017-07_mfd@H=4kOe'])
% printPNG([todaystr '_TmVO4-needles-2017-07_mfd@Tr=1p4'])

%% Compute the average value of ratio of internal to external magnetic field 
for i=1:length(Smfd_ndl)
Smfd_ndl(i).Hinm_Oe = round(...
    sum(...
    Smfd_ndl(i).binCenters.*Smfd_ndl(i).hc.*Smfd_ndl(i).binWidths...
    )*Smfd_ndl(i).Hext_Oe...
    );
end

%% Check results for given range of temperatures/field
rng = param_range{1};% define range

% Check that the ratio of Hinm/Hext is roughly constant inside the
% ordered phase, thus allowing to use the value at 0.3 K and 1000 Oe as a proxy
T = [Smfd_ndl(rng).T_K]
Hext = [Smfd_ndl(rng).Hext_Oe]
Hinm = [Smfd_ndl(rng).Hinm_Oe]
h = Hinm./Hext

rescaling = Smfd_ndl(1).Hinm_Oe/Smfd_ndl(1).Hext_Oe;%Hc0/0.69;% rescaling factor, due to demag; 


%% Import HC data
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2017-07_TmVO4_Cp_MCE\2017-07-28--31\2017-07-28_Cp'
NdlData=ImportTmVO4Cp('TmVO4_Mosaic_2017-07-28.dat');
% Use this data to plot Cp vs (H, T) superimposed on theoretical curves
% Average data points at each temperature + field 

%% Store data into variables
% H=[DATA.FieldOersted; DATA(2).FieldOersted];
% T=[DATA.SampleTempKelvin; DATA(2).SampleTempKelvin];
% Cp=[DATA.SampHCJmoleK; DATA(2).SampHCJmoleK];
Hndl=[NdlData.FieldOersted];
Tndl=[NdlData.SampleTempKelvin];
CpNdl=[NdlData.SampHCJmoleK];
CpErrNdl=[NdlData.SampHCErrJmoleK];

whichPoints = isfinite(Hndl) & isfinite(Tndl) & isfinite(CpNdl);
Hndl=Hndl(whichPoints);
Tndl=Tndl(whichPoints);
CpNdl=CpNdl(whichPoints);
CpErrNdl=CpErrNdl(whichPoints);

hmax=5500;

%% Plot 3D scatter
scatter3(Hndl,Tndl,CpNdl,'.')
xlim([0 hmax])
ylim([0 3])
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('Cp (uJ/K)')

%% Prepare fit of 3D data
tstep=0.05;
hstep=500;
[Hg,Tg] = meshgrid(0:hstep:hmax,0.365:tstep:3);% for griddata
Hgl=0:hstep:hmax;% for gridfit
Tgl=0.365:tstep:3;% for gridfit

%% Plot fit of 3D data
figure
Cpg = gridfit(Hndl,Tndl,CpNdl,Hgl,Tgl);
%Cpg = griddata(H,T,Cp,Hg,Tg);
surf(Hg,Tg,Cpg)
hold on
scatter3(Hndl,Tndl,CpNdl,100,'.','.k')
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('Cp (uJ/K)')

%% Gaussian smoothing
steps = -50:50;
x= tstep*steps;
s=0.05;
d1Gaussian = -exp(-x.^2/(2*s^2)).*x./sqrt(s^6*2*pi);
d2Gaussian = exp(-x.^2/(2*s^2)).*(x.^2 - s^2)/sqrt(s^10*2*pi);

d1Cpg = conv2(Cpg,d1Gaussian','same');
d2Cpg = conv2(Cpg,d2Gaussian','same');

%% Plot smoothed 3D surface
figure
surf(Hg,Tg,-d1Cpg,'EdgeColor','none');
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('-dCp/dT (uJ/K^2)')


%% Separate data according to value of magnetic field
clear separatedCpData

fieldsNdl = unique(round(Hndl,-1));%[10,linspace(1000,4000,4),4500,4750,5000];

for i = 1:length(fieldsNdl)
    wp = abs(Hndl-fieldsNdl(i))<50;
    separatedNdlCpData(i).H = Hndl(wp);
    separatedNdlCpData(i).T = Tndl(wp);
    separatedNdlCpData(i).Cp = CpNdl(wp);
    separatedNdlCpData(i).CpErr = CpErrNdl(wp);

    [separatedNdlCpData(i).T,wo] = sort(separatedNdlCpData(i).T);
    separatedNdlCpData(i).H = separatedNdlCpData(i).H(wo);
    separatedNdlCpData(i).Cp = separatedNdlCpData(i).Cp(wo);
    separatedNdlCpData(i).CpErr = separatedNdlCpData(i).CpErr(wo);
end

%% Average datapoints taken at any given temperature and field setpoint
% We want to compute the average of data points that are repetitions of the
% same measurement, i.e. with same temperature and field setpoints
clear avgNdlData
Tsep = 6e-3;% Data points taken within a relative interval of Tsep are considered to be measured at the same temperature setpoint
Tp = 27;% Temperature scale of phonons contribution to Cp in TmVO4, in K; see 'TmVO4_Cp_phonons.m'
A = 31.5;% 31.5 is the value of amplitude coefficient extracted from fit using curve fitting tool
for i=length(fieldsNdl):-1:1
avgNdlData(i) = averageCpwithH2(Tsep,separatedNdlCpData(i).T,separatedNdlCpData(i).Cp,...
    separatedNdlCpData(i).CpErr,separatedNdlCpData(i).H);
end
for i=length(fieldsNdl):-1:1
avgNdlData(i).Cpel = avgNdlData(i).Cp - R*(avgNdlData(i).T/Tp).^3;% electronic contribution to Cp, after subtracting phonons contribution
avgNdlData(i).Cpelr = avgNdlData(i).Cpel/A;
avgNdlData(i).CpelrErr = avgNdlData(i).CpFullErr/A;
avgNdlData(i).uh = unique(round(avgNdlData(i).H,-1));
end

%% Plot averaged data
figure; hold on
for i=1:length(fieldsNdl)
    errorbar(avgNdlData(i).T,avgNdlData(i).Cpel,avgNdlData(i).CpFullErr,'.','MarkerSize',18)
end
xlabel('Temperature (K)');
ylabel('$C_p$ (J/K/mol)')

%% Prepare MF fit of Cp vs Temperature under field
% Useful to use the curve fitting tool for quick look at the data
index = 1;
H1 = fieldsNdl(index);
T1 = avgNdlData(index).T;
Cp1 = avgNdlData(index).Cpelr;
Cp1Err = avgNdlData(index).CpelrErr;
wghts1 = 1./Cp1Err;

%% Longitudinal field value extracted from curve fitting Cp data at 10 Oe
% Cp_LFIM(T/Tc,e) fit of data at 10 Oe (excluding datapoints at T=2.252K and above) yields:
%        Tc =       2.199  (2.196, 2.203)
%        e =    0.001087  (0.0008557, 0.001319)
e = 1.1e-3;
Tc0_ndl = 2.2;

%% Create table from structure
Tmfd_ndl1 = struct2table(Smfd_ndl(1:56));% 
utmfd1 = unique(Tmfd_ndl1.T_K);
uhmfd1 = unique(Tmfd_ndl1.Hext_Oe);

%% Compute Cp for Gaussian distribution of fields
% clear Cpnum
% Hcnum=4900
rescaling = .97;
rngNum = 1:length(fieldsNdl);
for i=rngNum
Cpnum(i).h = fieldsNdl(i)/Hc*rescaling;
Cpnum(i).t_single_h = linspace(0,1.4,561);% reduced temperature, T/Tc
Cpnum(i).single_h_no_e = zeros(size(Cpnum(i).t_single_h));
Cpnum(i).single_h_w_e = zeros(size(Cpnum(i).t_single_h));
Cpnum(i).t_h_dist_no_e = [linspace(0.01,.5,50) linspace(0.505,1.1,120) linspace(1.11,1.4,30)];% reduced temperature, T/Tc
Cpnum(i).comsolpdf_no_e = zeros(size(Cpnum(i).t_h_dist_no_e));
Cpnum(i).comsolpdf_w_e = zeros(size(Cpnum(i).t_h_dist_no_e));
end

%% Compute Cp_TLFIM at single value of field
tic
for i=rngNum
Cpnum(i).single_h_no_e = Cp_TFIM(Cpnum(i).t_single_h,Cpnum(i).h);
[~,~,Cpnum(i).single_h_w_e] = FSCp_TLFIM(Cpnum(i).t_single_h,Cpnum(i).h,e);
Cpnum(i).single_h_w_e = Cpnum(i).single_h_w_e';
end
toc

%% Plot averaged data at each field separately
figure; hold on
rngAvg = [1,3:6,8];%1:length(fieldsNdl);
clr = cell(size(rngAvg));
eb = cell(size(rngAvg));
for i=rngAvg
%     fp = fplot(@(t) Cp_TFIM(t/Tc_ndl,fieldsNdl(i)/Hc),[0 3],'LineWidth',2);
    fp = plot(Cpnum(i).t_single_h(2:end-1)*Tc0_ndl,Cpnum(i).single_h_w_e,'DisplayName','MF');
    clr{rngAvg==i} = get(fp,'Color');
end
for i=rngAvg
    eb{rngAvg==i} = errorbar(avgNdlData(i).T,avgNdlData(i).Cpelr,avgNdlData(i).CpelrErr,...
        '.','MarkerSize',18,'DisplayName',sprintf('%d',fieldsNdl(i)),...
        'Color',clr{rngAvg==i},'LineWidth',2);
end
xlabel('$T$ (K)'); 
ylabel('$C_p/R$');
% title(['Heat capacity of needles of TmVO$_4$, $e=$ ', sprintf('%.2g',e),...
%     ', $H\times$' sprintf('%.2g',rescaling)])
lgd = legend([eb{:}],'Location','northeast');
lgd.Title.String = '$H$ (Oe)';
hold off
annfield = annotation('textbox',[0.15 0.8 0.2 0.1], 'interpreter','latex',...
    'String', ['Fits at $H\times$' sprintf('%.2g',rescaling)], 'EdgeColor','none',...
    'FitBoxToText','on', 'BackgroundColor','none', 'Color','k');% add annotation
xlim([0 3.1]);

%% 
% formatFigure;
% printPNG([todaystr mfd_ID '_Cp_vs_T_H'])
printPDF([todaystr '_TmVO4-needles_Cp+single-h-fits_e=1p1e-3']);




%% Garbage: For a given dataset, find closest values of temperature and field in COMSOL mfd
i = 6;
Cpnum(i).comsolpdf_no_e = zeros(size(Cpnum(i).t_h_dist_no_e));
[~,mfdhidx] = min(abs(Tmfd_ndl1.Hext_Oe-fieldsNdl(i)));
Hext_mfd = Tmfd_ndl1.Hext_Oe(mfdhidx);
tref = [0,0];
t = zeros(2,length(Cpnum(i).t_h_dist_no_e));
wt = zeros(2,length(Cpnum(i).t_h_dist_no_e));
Trcomp = zeros(length(Tmfd_ndl1.T_K),length(Cpnum(i).t_h_dist_no_e));
for j=1:length(Cpnum(i).t_h_dist_no_e)
    % Find values of temperature in COMSOL mfd closest to that of interest
    Trcomp(:,j) = Tmfd_ndl1.T_K/Tc_ndl-Cpnum(i).t_h_dist_no_e(j);
    [abst,mfdtidx] = unique(abs(Trcomp(:,j)));

    % if the 2 closest temperatures are both above or below the one of
    % interest, just use the single closest, otherwise use both
    if sign(Trcomp(mfdtidx(1),j))==sign(Trcomp(mfdtidx(2),j))
        t(:,j) = Tmfd_ndl1.T_K(mfdtidx(1));
        wt(:,j) = [1,0];
    else
        t(:,j) = Tmfd_ndl1.T_K(mfdtidx(1:2));
        wt(:,j) = 1-abst(1:2)/sum(abst(1:2));
    end

    if ~all(t(:,j)==tref)
        sprintf('j=%i, T=%.2gK, Tref=[%.2g,%.2g]K',j,Cpnum(i).t_h_dist_no_e(j)*Tc_ndl,t(:,j))
        tref=t(:,j);
    end
    % Same for value of field
    % Find the rows in Tmfd that matches both t and h 
    rows = find(ismember(Tmfd_ndl1.T_K,t(:,j)) & Tmfd_ndl1.Hext_Oe==Hext_mfd);
    ndl_temps = [rows(rows<=56),rows(rows>57)];
    [ntemps,nsamples] = size(ndl_temps);
    Cphnoe = zeros(size(Tmfd_ndl1.binCenters(rows,:)));
    for ndl_idx=1:nsamples
        for temp_idx=1:ntemps
            for Hin=1:length(Tmfd_ndl1.binCenters(ndl_temps(1,ndl_idx),:))
                Cphnoe(ndl_idx,Hin) = Cphnoe(ndl_idx,Hin) +...
                    Cp_TFIM(...
                    Cpnum(i).t_h_dist_no_e(j),...
                    Tmfd_ndl1.binCenters(ndl_temps(temp_idx,ndl_idx),Hin)...
                    ).*...
                    wt(temp_idx,j);
            end
        % Compute the corresponding value of heat capacity
        Cpnum(i).comsolpdf_no_e(j) = Cpnum(i).comsolpdf_no_e(j) +...
            sum(...
            Cphnoe(ndl_idx,:).*...
            Tmfd_ndl1.hc(ndl_temps(temp_idx,ndl_idx),:).*...
            Tmfd_ndl1.binWidths(ndl_temps(temp_idx,ndl_idx),:)...
            )./...
            prod(size(ndl_temps));
        end
    end
%     Cph_ndl = mean(Cph,1);

end

%% For a given dataset, find closest values of temperature and field in COMSOL mfd
tic

% for i=2%rngNum(3:end)
    
    Hext_data = unique(round(avgNdlData(i).H,-1));
    [~,mfdhidx] = min(abs(Tmfd_ndl1.Hext_Oe-Hext_data));
    Hext_mfd = Tmfd_ndl1.Hext_Oe(mfdhidx);
    trange = Cpnum(i).t_h_dist_no_e;
    tref = zeros(2,length(trange)+1);% Not even necessary given the the function find_mfd_temp_rows() includes [0 0] as default values for tref
    % For the computations with longitudinal field, need to compute free
    % energy first, since there is no analytical formula for the heat capacity
    Fwewt = zeros(2,length(trange));% Initialize array of weighted partial free energies
    Cpnoewt = zeros(2,length(trange));% Initialize array of weighted partial heat capacities (w/o longitudinal strain)
    % t = zeros(2,length(Cpnum(i).t_h_dist_no_e));% array that will contain two closest values of temperature
    wt = zeros(2,length(trange));% array that will contain weights attributed to closest temperatures
    ndl_rows = zeros(2,length(trange));
    Cphnoe = zeros(size(Tmfd_ndl1.binCenters(1,:)));

    for jt=1:length(trange)
    %     [ndl1_rows, wt(:,jt)] = find_mfd_temp_rows(trange(jt),...
    %         utmfd1, Tc0_ndl, Tmfd_ndl1, Hext_mfd );
        % add tref as input and output of function find_mfd_temp_rows() to check 
        % that the function finds the correct temperatures to compute Cpnum
        [ndl1_rows, wt(:,jt), tref(:,jt+1)] = find_mfd_temp_rows(trange(jt),...
            utmfd1, Tc0_ndl, Tmfd_ndl1, Hext_mfd, tref(:,jt), 'printTref', true );

        ndl_rows(:,jt) = ndl1_rows;

        % Compute heat capacity at each value of internal (transverse) 
        % magnetic field for both nearest temperatures
        for jr=1:length(ndl1_rows)

            % Compute free energy, including longitudinal field
            % Note: it is important to compute this array line by line using
            % the for loop over jr, because the computation over the entire
            % array at once induces computation errors, which are hard to
            % spot at first, because they only arise at the 5th significant
            % digit, but this is enough to make a big difference in the computation
            % of the second derivative, i.e. the heat capacity
            Fhwe(jr,:) = FSCp_TLFIM(trange(jt),Tmfd_ndl1.binCenters(ndl1_rows(jr),:),e);

            % Compute Cp without longitudinal field, for each value of internal field
            for col=1:length(Tmfd_ndl1.binCenters(ndl1_rows,:))
                Cphnoe(jr,col) = Cp_TFIM(trange(jt),...
                    Tmfd_ndl1.binCenters(ndl1_rows(jr),col));
            end

            % Sum over all internal fields to get the total heat capacity at
            % the temperature associated with jt
    %        Cpnoe(jr,jt) = wt(jr,jt)*sum(...
            Cpnoewt(jr,jt) = wt(jr,jt).*sum(...
                Tmfd_ndl1.binWidths(ndl1_rows(jr),:).*...
                Cphnoe(jr,:).*...
                Tmfd_ndl1.hc(ndl1_rows(jr),:)...
                );
        end

        % Sum the free energy over all internal fields to get the total free energy at
        % the temperature associated with jt
    %     Fwewt(:,jt) = sum(...
        Fwewt(:,jt) = sum(...% Use this if trange does not cover the full range of Cpnum(i).t_h_dist_no_e
            Tmfd_ndl1.binWidths(ndl1_rows,:).*...
            Fhwe.*...
            Tmfd_ndl1.hc(ndl1_rows,:),...
            2);

    end

    toc
    
    %%
    % Cpnum(i).t_h_dist_no_e = linspace(0.1,1.5,150);% reduced temperature, T/Tc
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

    %% Store results into Cpnum structure
    Cpnum(i).t_h_dist_w_e = twe(2:end-1)';
    Cpnum(i).comsolpdf_no_e = Cpnoe;
    Cpnum(i).comsolpdf_w_e = Cpwe;

    %% Prepare data for plotting by removing NaN datapoints
    sel_no_e = ~isnan(Cpnum(i).comsolpdf_no_e);
    sel_w_e = ~isnan(Cpnum(i).comsolpdf_w_e);
    Cpnum(i).t_h_dist_no_e = Cpnum(i).t_h_dist_no_e(sel_no_e);
    Cpnum(i).t_h_dist_w_e = Cpnum(i).t_h_dist_w_e(sel_w_e);
    Cpnum(i).comsolpdf_no_e = Cpnum(i).comsolpdf_no_e(sel_no_e);
    Cpnum(i).comsolpdf_w_e  = Cpnum(i).comsolpdf_w_e(sel_w_e);

% end

% toc 

%% Plot Cp for COMSOL distribution of fields
i=7;
single_str = '$H=H_{\mathrm{ext}}$';
figure
% plot(avgNdlData(i).T,avgNdlData(i).Cpelr,'.','DisplayName','data')
eb{rngAvg==i} = errorbar(avgNdlData(i).T,avgNdlData(i).Cpelr,avgNdlData(i).CpelrErr,...
    '.','MarkerSize',18,'DisplayName',['Data at ' single_str],...
    'LineWidth',2);
hold on;
no_e_str = ' ($e=0$)';
w_e_str = [' ($e=$ ' sprintf(' %.2g)',e)];

% Plot results at single value of magnetic field
plot(Cpnum(i).t_single_h*Tc0_ndl,Cpnum(i).single_h_no_e,...
    'DisplayName',[single_str no_e_str]);% no longitudinal field
plot(Cpnum(i).t_single_h(2:end-1)*Tc0_ndl,Cpnum(i).single_h_w_e,...
    'DisplayName',[single_str w_e_str]);% with longitudinal field

% Plot results for distribution of magnetic fields as computed with COMSOL
comsol_str = 'comsol pdf';
plot(Cpnum(i).t_h_dist_no_e*Tc0_ndl,Cpnum(i).comsolpdf_no_e,...
    'DisplayName',[comsol_str no_e_str]);
plot(Cpnum(i).t_h_dist_w_e*Tc0_ndl,Cpnum(i).comsolpdf_w_e,...
    'DisplayName',[comsol_str w_e_str]);
title(['$C_p$ single field vs COMSOL PDF $H_{\mathrm{ext}}=$ ' sprintf('%.0f Oe',fieldsNdl(i))]);
lgd = legend('Location','best');% title(lgd,'TmVO4-Ndl-E');

xlabel('$T$ (K)')
ylabel('$C_p/R$')



%% Export figure
formatFigure;
% printPNG([todaystr...
%     sprintf('_TmVO4-2017-07-needles_Cp_vs_T_@%iOe_single-field_vs_comsol-pdf_fits',...
%     fieldsNdl(i))]);
% printPDF(['2019-06-18_TmVO4-RF-E_fit_Schottky_' strrep(hrstr,'.','p') 'xHc']);

%% Export results to .dat file
cd '2020-08_Full_Cpnum_COMSOL_PDF_w_e_data'
for i=3:8
    A = cat(2,Cpnum(i).t_h_dist_w_e, Cpnum(i).comsolpdf_w_e');
    col_headers = ['T/Tc, Cp/R'];
    file_description = {...
        [sprintf('Heat capacity at Hext = %d Oe, computed numerically ',fieldsNdl(i))...
        'for the largest needle of the TmVO4 "mosaic" measured on 2017-07-28.'];
        ['The computation took into account the distribution of magnetic fields '...
        'computed using COMSOL on 2020-07-31, as well as the rounding '...
        'of the transition due to a constant longitudinal field (strain) '...
        'of e = 1.1e-3 as obtained from fitting the data at zero field.'];
        ['The computed heat capacity is a weighted sum of the heat capacity '...
        'at the two closest temperatures at which the distribution of magnetic '...
        'fields was computed in COMSOL'];
        ['The zero field transition temperature obtained from the same fit '...
        'was Tc(H=0) = 2.20K.'];
        'Any missing data is due to removal of NaN and computational aberrations.'};

    fname = [todaystr...
        sprintf('_TmVO4-2017-07-needles_Cpnum_comsol_pdf_@%dOe.dat',...
        fieldsNdl(i))];
    fileID = fopen(fname,'w');
    fprintf(fileID, '# %s\n', file_description{:});
    fprintf(fileID, '%s\n', col_headers);
    fclose(fileID);
    dlmwrite(fname, A, '-append')
end




%% Older data analysis: dCp/dT, phase boundary

%% Gaussian convolution and computation of derivative of raw Cp
tstep=0.025; 
steps = -50:50;
x= tstep*steps;
s=0.05;
d1Gaussian = -exp(-x.^2/(2*s^2)).*x./sqrt(s^6*2*pi);
d2Gaussian = exp(-x.^2/(2*s^2)).*(x.^2 - s^2)/sqrt(s^10*2*pi);
d1Cp = conv2(CpNdl,d1Gaussian','same');

%% Compute and plot derivative of averaged data
figure; hold on
for i = 1:length(fieldsNdl)
    separatedNdlCpData(i).d1Cpm = conv2(avgNdlData(i).Cp,d1Gaussian','same');
    plot(avgNdlData(i).T,separatedNdlCpData(i).d1Cpm)
end

%% Identify experimental critical temperature at each field
% Then we can plot Tc vs fields and fit using equation f(h) = Tc/Hc*h/atanh(h/Hc)
% From this, we can see that the experimental value of Hc is ~0.72T 
% and hence correct for demag, knowing the value of critical field Hc0 (see beginning of code)
M = ones(1,length(fieldsNdl));
Tc_ndl = ones(1,length(fieldsNdl));
for i=1:length(fieldsNdl)
    [M(i),I] = min(separatedNdlCpData(i).d1Cpm);
    Tc_ndl(i) = avgNdlData(i).T(I);
end


%% Plot data at each field separately
% Probably useless now that the average is computed. 2019-03-18
figure
for i=1:length(fieldsNdl)
    plot(separatedNdlCpData(i).T,separatedNdlCpData(i).Cp,'.-')
    hold on
end
xlabel('Temperature (K)')
ylabel('Cp (uJ/K)')
legendCell = cellstr(num2str(fieldsNdl, '%-d Oe'));
legend(legendCell)
hold off

%% Fit data with smoothing spline
figure
for i=1:length(fieldsNdl)
    separatedNdlCpData(i).f=fitSpline(separatedNdlCpData(i).T,-separatedNdlCpData(i).Cp);
% fitSpline is a custom function? Press ctrl+D while cursor is on the name for more info
%    separatedCpData(i).d2=differentiate(separatedCpData(i).f,separatedCpData(i).T);
    plot(separatedNdlCpData(i).f,'deriv1')
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
for i=1:length(fieldsNdl)
    separatedNdlCpData(i).f=fitSpline(separatedNdlCpData(i).T,separatedNdlCpData(i).Cp);
    separatedNdlCpData(i).d2=differentiate(separatedNdlCpData(i).f,separatedNdlCpData(i).T);
    scatter3(separatedNdlCpData(i).H,separatedNdlCpData(i).T,-separatedNdlCpData(i).d2,48,'filled')
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
for i=1:length(fieldsNdl)
    T2= [T2; separatedNdlCpData(i).T];
    H2=[ H2; separatedNdlCpData(i).H];
    d2Cp2=[d2Cp2; separatedNdlCpData(i).d2];
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
M = ones(1,length(fieldsNdl));
Tc_ndl = ones(1,length(fieldsNdl));
for i=1:length(fieldsNdl)
    [M(i),I] = min(separatedNdlCpData(i).d2);
    Tc_ndl(i) = separatedNdlCpData(i).T(I);
end

%% Fit: 'Phase boundary'.
[xData, yData] = prepareCurveData( fieldsNdl', Tc_ndl );

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







