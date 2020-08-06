Tc_ndl = 2.22;% transition temperature at zero field, in Kelvin units
Hc = 5000;% critical field at zero temperature, in Oersted units

%% Import magnetic field distribution from CSV file
cd 'C:\Users\Pierre\Desktop\Postdoc\Software\COMSOL\TmVO4-LS5200_HC2017-07\TmVO4-LS5200_HC2017-07_COMSOL_results'
% mfdfilename = cell(1,2);
% S = cell(1,2);
mfdNdlFname{1} = 'TmVO4-LS5200_HC2017-07_sample1-Needle_mesh-p1mm_T=p3-p4-3p1_Hext=all.csv';
mfdNdlFname{2} = 'TmVO4-LS5200_HC2017-07_sample2-Arya_mesh-20um_T=p3-p4-3p1_Hext=all.csv';
mfdDomNum = {2,3};
for ic=1:length(mfdNdlFname)
    Sndl{ic} = importfielddistrib_csv(mfdNdlFname{ic}, 56, 'domainNum', mfdDomNum{ic});
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
    if sidx<=56
        Smfd_ndl(sidx).label = 'needle1';
    else
        Smfd_ndl(sidx).label = 'needle2';
    end
end

%% Compute probability distribution of fields at a given value of T and Hext
for i=1:length(Smfd_ndl)
mfd = Smfd_ndl(i).mfd(Smfd_ndl(i).mfd>0)/Hc;% create distribution from non-zero values
% h = histogram(mfd, 'Normalization', 'pdf');% plot histogram of distribution
[Smfd_ndl(i).hc, edges] = histcounts(mfd, 50, 'Normalization', 'pdf');% plot histogram of distribution
% binCenters = h.BinEdges + (h.BinWidth/2);
Smfd_ndl(i).binCenters = mean([edges(1:end-1);edges(2:end)],1);
Smfd_ndl(i).binWidths = edges(2:end)-edges(1:end-1);
end

%% Plot distribution of fields at a given value of T and Hext
figure;
hold on
needle_index = 1;% there are 2 needle-shaped samples 
start_index = [0,56];% needle 1 data start at index 0+1, needle 2 data at 56+1
param_index = 2;% 1 is constant T, 2 is constant Hext, see param_range
param_range = {[3:8:56], 33+[0:7]};% first range corresponds to a 
% field dependence at constant temp, second range corresponds to a 
% temperature dependence at constant field
rng = start_index(needle_index)  + param_range{param_index};
for i=rng 
Tndl = Smfd_ndl(i).T_K;%
Hext = Smfd_ndl(i).Hext_Oe;%
p = plot(Smfd_ndl(i).binCenters, Smfd_ndl(i).hc, '.-', 'DisplayName', sprintf('%.2g, %.2g',Tndl/Tc_ndl,Hext/Hc));
end
lgd = legend('show','Location','northwest'); lgd.Title.String = '$T/T_c$, $H_{\mathrm{ext}}/H_c$';
needle_title = sprintf('Distribution of fields in TmVO$_4$ needle %d', needle_index);
param_title = {[', $T=$' sprintf('%.2g K',Tndl)],...
    [', $H_{\mathrm{ext}}=$' sprintf('%.2d Oe',Hext)]};
title([needle_title param_title{param_index}])
xlabel('$H_{\mathrm{in}}/H_{c}$')
ylabel('Normalized PDF')

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


%%
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
for i=1:2:length(fieldsNdl)
    errorbar(avgNdlData(i).T,avgNdlData(i).Cpel,avgNdlData(i).CpFullErr,'.','MarkerSize',18)
end
xlabel('Temperature (K)');
ylabel('$C_p$ (J/K/mol)')

%% Prepare MF fit of Cp vs Temperature under field
% Useful to use the curve fitting tool for quick look at the data
index = 1;
H1 = fieldsNdl(index);
T1 = separatedNdlCpData(index).Tm;
Cp1 = separatedNdlCpData(index).Cpm;
Cp1Err = separatedNdlCpData(index).CpmFullErr;
wghts1 = 1./Cp1Err;

%% Plot averaged data at each field separately
figure; hold on
rngAvg = 1:length(fieldsNdl)-1;
clr = cell(size(rngAvg));
eb = cell(size(rngAvg));
for i=rngAvg
    fp = fplot(@(t) Cp_TFIM(t/Tc_ndl,fieldsNdl(i)/Hc),[0 3],'LineWidth',2);
    clr{rngAvg==i} = get(fp,'Color');
end
for i=rngAvg
    eb{rngAvg==i} = errorbar(avgNdlData(i).T,avgNdlData(i).Cpelr,avgNdlData(i).CpelrErr,...
        '.','MarkerSize',18,'DisplayName',num2str(fieldsNdl(i)/1e4,'%.2f T'),...
        'Color',clr{rngAvg==i},'LineWidth',2);
end
xlabel('Temperature (K)'); 
ylabel('C$_p/R$');
title('Heat capacity of needles of TmVO4 (no demag) at various fields')
legend([eb{:}]);
hold off



%% Create table from structure
Tmfd_ndl = struct2table(Smfd_ndl);% 
utmfd = unique(Tmfd_ndl.T_K);
uhmfd = unique(Tmfd_ndl.Hext_Oe);

%% Compute Cp for Gaussian distribution of fields
clear Cpnum
% Hcnum=4900
rngNum=1:8;
for i=rngNum
Cpnum(i).h = fieldsNdl(i)/Hc;
Cpnum(i).t_single_h = linspace(0,1.5,601);% reduced temperature, T/Tc
Cpnum(i).single_h = zeros(size(Cpnum(i).t_single_h));
Cpnum(i).t_phenomeno = linspace(0,1.5,301);% reduced temperature, T/Tc
Cpnum(i).comsolpdf = zeros(size(Cpnum(i).t_phenomeno));
end

for i=rngNum
Cpnum(i).single_h = Cp_TFIM(Cpnum(i).t_single_h,Cpnum(i).h);
end

%% Garbage: For a given dataset, find closest values of temperature and field in COMSOL mfd
i = 6;
Cpnum(i).comsolpdf = zeros(size(Cpnum(i).t_phenomeno));
[~,mfdhidx] = min(abs(Tmfd_ndl.Hext_Oe-fieldsNdl(i)));
h = Tmfd_ndl.Hext_Oe(mfdhidx);
tref = [0,0];
t = zeros(2,length(Cpnum(i).t_phenomeno));
wt = zeros(2,length(Cpnum(i).t_phenomeno));
Trcomp = zeros(length(Tmfd_ndl.T_K),length(Cpnum(i).t_phenomeno));
for j=1:length(Cpnum(i).t_phenomeno)
    % Find values of temperature in COMSOL mfd closest to that of interest
    Trcomp(:,j) = Tmfd_ndl.T_K/Tc_ndl-Cpnum(i).t_phenomeno(j);
    [abst,mfdtidx] = unique(abs(Trcomp(:,j)));

    % if the 2 closest temperatures are both above or below the one of
    % interest, just use the single closest, otherwise use both
    if sign(Trcomp(mfdtidx(1),j))==sign(Trcomp(mfdtidx(2),j))
        t(:,j) = Tmfd_ndl.T_K(mfdtidx(1));
        wt(:,j) = [1,0];
    else
        t(:,j) = Tmfd_ndl.T_K(mfdtidx(1:2));
        wt(:,j) = 1-abst(1:2)/sum(abst(1:2));
    end

    if ~all(t(:,j)==tref)
        sprintf('j=%i, T=%.2gK, Tref=[%.2g,%.2g]K',j,Cpnum(i).t_phenomeno(j)*Tc_ndl,t(:,j))
        tref=t(:,j);
    end
    % Same for value of field
    % Find the rows in Tmfd that matches both t and h 
    rows = find(ismember(Tmfd_ndl.T_K,t(:,j)) & Tmfd_ndl.Hext_Oe==h);
    ndl_temps = [rows(rows<=56),rows(rows>57)];
    [ntemps,nsamples] = size(ndl_temps);
    Cph = zeros(size(Tmfd_ndl.binCenters(rows,:)));
    for ndl_idx=1:nsamples
        for temp_idx=1:ntemps
            for Hin=1:length(Tmfd_ndl.binCenters(ndl_temps(1,ndl_idx),:))
                Cph(ndl_idx,Hin) = Cph(ndl_idx,Hin) +...
                    Cp_TFIM(...
                    Cpnum(i).t_phenomeno(j),...
                    Tmfd_ndl.binCenters(ndl_temps(temp_idx,ndl_idx),Hin)...
                    ).*...
                    wt(temp_idx,j);
            end
        % Compute the corresponding value of heat capacity
        Cpnum(i).comsolpdf(j) = Cpnum(i).comsolpdf(j) +...
            sum(...
            Cph(ndl_idx,:).*...
            Tmfd_ndl.hc(ndl_temps(temp_idx,ndl_idx),:).*...
            Tmfd_ndl.binWidths(ndl_temps(temp_idx,ndl_idx),:)...
            )./...
            prod(size(ndl_temps));
        end
    end
%     Cph_ndl = mean(Cph,1);

end

%% For a given dataset, find closest values of temperature and field in COMSOL mfd
i = 6;
Hdata = unique(round(avgNdlData(i).H,-2));
[~,mfdhidx] = min(abs(Tmfd_ndl.Hext_Oe-Hdata));
h = Tmfd_ndl.Hext_Oe(mfdhidx);
tref = 0;
for j=1:length(Cpnum(i).t_phenomeno)
    % Find value of temperature in COMSOL mfd closest to that of actual data
    [~,mfdtidx] = min(abs(Tmfd_ndl.T_K/Tc_ndl-Cpnum(i).t_phenomeno(j)));
    % Improvement note: use sort instead of min, to be able to interpolate...
    t = Tmfd_ndl.T_K(mfdtidx);
    if t ~= tref
        sprintf('j=%i, T=%.2gK',j,t)
        tref=t;
    end
    % Same for value of field
    % Find the row in Tmfd that matches both t and h 
    row = find(Tmfd_ndl.T_K==t & Tmfd_ndl.Hext_Oe==h);
    ndl1_row = row(1);
    Cph = zeros(size(Tmfd_ndl.binCenters(ndl1_row,:)));
    for col=1:length(Tmfd_ndl.binCenters(ndl1_row,:))
        Cph(col) = Cp_TFIM(Cpnum(i).t_phenomeno(j),Tmfd_ndl.binCenters(ndl1_row,col));
    end
    % Compute the corresponding value of heat capacity 
%     Cpnum(i).comsolpdf(j) = trapz(Tmfd_ndl.binCenters(row,:),Cph.*Tmfd_ndl.hc(row,:));
    Cpnum(i).comsolpdf(j) = sum(...
        Tmfd_ndl.binWidths(ndl1_row,:).*Cph.*...
        Tmfd_ndl.hc(ndl1_row,:));
end

%% Plot Cp for COMSOL distribution of fields
figure
plot(avgNdlData(i).T,avgNdlData(i).Cpelr,'.','DisplayName','data')
hold on;
plot(Cpnum(i).t_single_h*Tc_ndl,Cpnum(i).single_h,'DisplayName','MF');
plot(Cpnum(i).t_phenomeno*Tc_ndl,Cpnum(i).comsolpdf,'DisplayName',sprintf('Hc=%.2dOe',h));
title(['Cp mean-field vs COMSOL pdf $H_{\mathrm{ext}}=$' sprintf('%.0fOe',fieldsNdl(i))]);
lgd = legend();% title(lgd,'TmVO4-Ndl-E');






%% Export figure
formatFigure;
% printPNG([todaystr '_TmVO4-2017-07-needle2_mfd@H=4750Oe']);
% printPDF(['2019-06-18_TmVO4-RF-E_fit_Schottky_' strrep(hrstr,'.','p') 'xHc']);






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







