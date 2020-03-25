%% Import neutrons diffraction data measured at 0.6K
% Import first set of data (format is different than the second set)
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_neutrons\2019-02_ORNL_Corelli\2019-02-14';
clear nData;
FilesOrtho = dir('p6K_*.txt');
for i=1:length(FilesOrtho)
    tbl1 = ImportNeutronsLineCut(FilesOrtho(i).name);
    nData(i).file = FilesOrtho(i).name;
    nData(i).hh0 = tbl1.hh0;
    nData(i).I = tbl1.I;
    nData(i).dI = tbl1.dI;
    [C,matches] = strsplit(FilesOrtho(i).name,{'p6K_','T\w*.txt'},...
        'DelimiterType','RegularExpression','CollapseDelimiters',true);
    % extract value of field from file name using regular expression
    str = replace(C{2},'p','.');% replace letter p (if any) with a decimal dot
    nData(i).field = str2num(str);
    nData(i).temp = 0.6;
end

%% Import second data set at 0.6K
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_neutrons\2019-02_ORNL_Corelli\2019-02-14\p6K\linecut2';
fieldinfo = importfieldinfo('field_info.txt');% fieldinfo = import manually data from "field_info.txt"
if exist('nData','var'); lnd = length(nData); else; lnd = 0; end
for j=1:length(fieldinfo.FileName)
    tbl2 = ImportNeutronsLineCut(strcat(fieldinfo.FileName(j,:),".txt"));
    i = j+lnd;
    nData(i).file = fieldinfo.FileName(j,:);
    nData(i).hh0 = tbl2.hh0;
    nData(i).I = tbl2.I;
    nData(i).dI = tbl2.dI;
    nData(i).field = fieldinfo.H_T(j);
    nData(i).temp = fieldinfo.T_K(j);
end

%% Sort structure according to ascendring value of magnetic field
nData = nestedSortStruct(nData,'field');%
% Note: nestedSortStruct is a Matlab add-on available at https://www.mathworks.com/matlabcentral/fileexchange/28573-nestedsortstruct

% %% Import neutrons diffraction data measured at 0.94K
% % fieldinfo0p94K = import manually data from "field_info.txt"
% fieldinfo0p94K.FileName = strcat("HH0_",num2str(fieldinfo0p94K.FileID));
% for i=1:length(fieldinfo0p94K.FileName)
%     tbl2 = ImportNeutronsLineCut(strcat(fieldinfo0p94K.FileName(i,:),".txt"));
%     nData(i).file = fieldinfo0p94K.FileName(i,:);
%     nData(i).hh0 = tbl2.hh0;
%     nData(i).I = tbl2.I;
%     nData(i).dI = tbl2.dI;
%     nData(i).field = fieldinfo0p94K.H_T(i);
%     nData(i).temp = fieldinfo0p94K.T_K(i);
% end

%% Magnetic field data
field = extractfield(nData,'field');
Hc_0 = 0.51;% value in Tesla units of the critical field at zero temperature
% in the absence of demagnetizing factor
% see data taken on needles of TmVO4-LS5200 in July 2017

%% Basic data treatment
nData(round(field,2)==0.86)=[];% Delete fields where the data does not make
% sense: H=0.86T and H=0.865T (which rounds to 0.86)
nData(round(field,2)==0).I = nData(round(field,2)==0).I*0.73/1.15;
% rescale data at zero field, as it has a higher intensity than the rest
field = extractfield(nData,'field');

%% Center of peak to be studied in the following
hcenter = -8.0;% center of unsplit peak in reciprocal space

%% Analysis parameters
% ufb = 0.99; % upper fit boundary = highest value of h-hc for which to include datapoints for the fit
istart = 1;
iend = length(field);
fullHrange = istart:1:iend;% full range of magnetic fields
aboveHc = iend;% range of magnetic fields above the critical field

%% Initialize fit using convolution of Ikeda-Carpenter function with pseudo-Voigt
clear fitPrms
I2 = 17e4;% fit parameters initial values
peakCenter = hcenter*(1-[-1:1]*1e-3);% position of peaks in reciprocal space
freeKeySet =        {'I',   'R',    'alpha','beta', 'gamma','sigma','x0'};%7 free parameters
for fitPeakIdx = 1:3% initialize fit parameter values
    freeValueSet =  {4e4,   0.1,    200,    0.1,    1e-3,   6.6e-3, peakCenter(fitPeakIdx)};%7 free parameters
    for keyIdx = 1:length(freeKeySet)
        fitPrms(fitPeakIdx).freePrms{keyIdx,1} = sprintf([freeKeySet{keyIdx} '%i'],fitPeakIdx);
        fitPrms(fitPeakIdx).freePrms{keyIdx,2} = freeValueSet{keyIdx};
    end
end
fitPrms(2).freePrms{1,2} = 8e4;%change initial intensity of central fit peak
% freePrms11 = {'I',I11;'R11',R1;'alpha11',a1;'beta',b1;'gamma',g1;'sigma',s1;'x011',peakCenter(1)};%7 free parameters
% freePrms1 = {'I1',I1;'R1',R1;'beta1',b1;'gamma1',g1;'x01',hc(1)};% 5 free parameters
% freePrms2 = {'I2',I2;'R2',R2;'beta2',b2;'gamma2',g2;'x02',hc(2)};% 5 free parameters
% freePrms2 = {'I2',I2;'x02',xc(2)};% free parameters
% Note: the order in which free parameters are defined matters because of
% how the array of initial fitting parameters 'initParams' is defined
fmt1 = 'ENS at T=%.2fK H=%.3fT %i params fit';% string format for plot title

%% Perform and plot fit 
for i=aboveHc
    nData1 = nData(i).hh0(nData(i).dI>0);
    datExcld = nData1<min(hcenter)-.7 | nData1>max(hcenter)+0.55 |...
        (nData1>min(hcenter)-.35 & nData1<min(hcenter)-0.15) |...
        (nData1>max(hcenter)+0.14 & nData1<max(hcenter)+0.2);%
%     (nData1>max(hc)+0.045 & nData1<max(hc)+0.065);% Exclude 
% data points that correspond to other peaks as well as those that are too far away
% datExcld should be defined on the x interval where obj.dY>0, otherwise
% the number of excluded points will be higher than the number of data points
    myfit = fitICpV(nData(i).hh0,nData(i).I,nData(i).dI,peakCenter); 
    myfit.dataExcl = datExcld;
%     ap1 = myfit.allParams{1}; ap1('gamma') = 0;% 'ap1' is shorter than 'myfit.allParams{1}'
%     ap1('alpha') = 140; ap1('sigma') = 6.6e-3;% It is not necessary to reassign the values that are assigned by default
%     ap1('R')=0.0; ap1('beta')=0; ap1('I') = I1; 
%     ap2 = myfit.allParams{2}; ap2('gamma') = 0;
    myfit.freeParams = {fitPrms.freePrms}
    myfit.indepFreePrms();% compute array of *independent* free parameters
    Nprms = length(myfit.indepFreeParams);% total number of independent free parameters
    label = sprintf(fmt1,nData(i).temp,field(i),Nprms);
    fitStr = ['fit'  int2str(length(peakCenter)) 'ICpV' int2str(Nprms)]; 
    gofStr = ['gof'  int2str(length(peakCenter)) 'ICpV' int2str(Nprms)];
    [nData(i).(fitStr), nData(i).(gofStr)] = myfit.compute_fit();
    if exist('H_c','var'); hfactor = H_c/cval(1);
    else hfactor = 1;
    end
    if 1%mod(i,20)==1% select data to plot
        myfit.plot_fit(nData(i).(fitStr));% title(label);
    xlim([hcenter-.1 hcenter+.1]);
    ann00 = annotation('textbox',[0.15 0.8 0.2 0.1],'interpreter','latex',...
        'String',{sprintf('T=%.2fK',nData(i).temp) sprintf('H=%.2fT',field(i)*hfactor)},...
        'FontSize',14,'FontName','Arial','LineStyle','-','EdgeColor','r',...
        'FitBoxToText','on','LineWidth',2,'BackgroundColor',[1 1 1],'Color','k');% add annotation
% hAnnotAxes = findall(gcf,'Tag','scribeOverlay');% retrieve annotation object, in case it is necessary to manipulate it
    end
    disp(label); disp(nData(i).(fitStr)); disp(nData(i).(gofStr));
end

%% Write fit parameters to a table 
np = Nprms;
fitStrnp = fitStr; gofStrnp = gofStr;
if ~exist('Stbl','var') || ~isfield(Stbl,(fitStrnp)); Stbl(1).(fitStrnp) = []; end
% if structure Stbl does not contain any field called (fitStrnp), create that field
flag = 0;% to check whether it is necessary to add a new row to the field or not
for i=1:numel(Stbl)
  if isempty(Stbl(i).(fitStrnp)); Ntbl = i; flag = 1; break; end%
  %if a row of Stbl.(fitStrnp) is empty, use it to store the following table
end
if ~flag; Ntbl = numel(Stbl)+1; end% if all the rows of Stbl.(fitStrnp) are non-empty, store table in a new row
Stbl(Ntbl).(fitStrnp) = table(field(fullHrange)','VariableNames',{'Field_Oe'});% first column of the table contains magnetic field values
cfn = coeffnames(nData(istart).(fitStrnp));% extract fit coefficient names, i.e. free parameters 
for nc=1:np% for each free parameter
    prm = extract_structure_field(nData(fullHrange),fitStrnp,cfn{nc});% get its value
    cft = cell2mat(arrayfun(@(c) confint(c.(fitStrnp)),nData(fullHrange).','Uniform',0));% get the 95% confidence intervals
    cftm = cft(1:2:end,nc); cftp = cft(2:2:end,nc);% get the lower and upper bounds of that interval, respectively
    prmErrm = prm-cftm;prmErrp = prm-cftp;% calculate corresponding negative and positive absolute errors
    relErrm = abs(prmErrm./prm);relErrp = abs(prmErrp./prm);% and relative errors
    Stbl(Ntbl).(fitStrnp).(cfn{nc}) = prm;% store fit parameter value in a new column
    Stbl(Ntbl).(fitStrnp).([cfn{nc} '_RelErr']) = relErrm;% and relative error in another new column
end
% if ~any(strcmp('R',Stbl(Ntbl).(fitStrnp).Properties.VariableNames))% if there is no field called 'R' in table Stbl.(fitStrnp)
%     Stbl(Ntbl).(fitStrnp).R = ap1('R')*ones(length(rng),1);% store value of parameter R in another column 
% % this is only useful when looking at the influence of the value of R on
% % the quality of fits
% end
rsquarenp = extract_structure_field(nData(fullHrange),gofStrnp,'rsquare');% extract r^2 value from goodness of fit
Stbl(Ntbl).(fitStrnp).Rsquare = rsquarenp;% store r^2 value in a new column
flag = 0;% reset flag for next run

%% Plot ratio of peak intensities as a function of field
if ismember('I1', Stbl(end).(fitStr).Properties.VariableNames) && ismember('I2', Stbl(end).(fitStr).Properties.VariableNames)
    figure; hold on;
    If1 = Stbl(end).(fitStr).I1; If1Err = Stbl(end).(fitStr).I1_RelErr.*If1;
    If2 = Stbl(end).(fitStr).I2; If2Err = Stbl(end).(fitStr).I2_RelErr.*If2;
    Iratio = If1./If2;
    IrErr = abs(If1Err.*If2-If2Err.*If1)./If2.^2;
    errorbar(Stbl(end).fit2ICpV4.Field_Oe,Iratio,IrErr,'.','MarkerSize',12);
    % errorbar(Stbl.fit2ICpV4.Field_Oe,If1,Stbl.fit2ICpV4.I1_RelErr,'.','MarkerSize',12);
    % errorbar(Stbl.fit2ICpV4.Field_Oe,If2,Stbl.fit2ICpV4.I2_RelErr,'.','MarkerSize',12);
    ylim([0 2*Iratio(1)]);
    title('TmVO$_4$ neutrons elastic 880 peak intensity ratio vs field at 0.6K');
    xlabel('Magnetic field (Oe)'); ylabel('I$_1$/I$_2$');
    nmaxavg = 10;
    annir = annotation('textbox',[0.15 0.8 0.2 0.1],'interpreter','latex',...
        'String',sprintf('Average of I$_1$/I$_2$ below %.2fT: %.2f',Stbl(end).fit2ICpV4.Field_Oe(nmaxavg),mean(Iratio(1:nmaxavg))),...
        'FontSize',14,'FontName','Arial','LineStyle','-','EdgeColor','r',...
        'LineWidth',2,'BackgroundColor',[1 1 1],'Color','k');% add annotation
    annir.FitBoxToText='on';% fit annotation box to text
else; warning('Need two intensities to plot ratio. Ignoring this part.');
end

%% Write table to file
fileChar = [sprintf('%1.eK',nData(1).temp) fitStr '.txt'];% '_R=' sprintf('%2.e',ap1('R'))
fileID = fopen(fileChar,'a');
% fprintf(fileID,'\nValues of free parameters after fit:\n');
% fprintf(fileID,'%s\n',nData(np).Stbl);
writetable(Stbl(Ntbl).(fitStrnp),fileChar);
fprintf(fileID,'\nInitial values of free parameters:\n');
S = string(vertcat(myfit.freeParams{:}));
fprintf(fileID,'%s\n',strcat(S(:,1)," = ",S(:,2)));
fprintf(fileID,'\nFixed values of all parameters (disregard free parameters):');
for il=1:length(peakCenter)
    fprintf(fileID,['\nPeak #' int2str(il) '\n']);
    K = keys(myfit.allParams{il}); V = values(myfit.allParams{il});
    fprintf(fileID,'%s\n',strcat(K," = ",string(V)));
end
fclose(fileID);

%% Identify single peak maximum and width
% xM = ones(length(rng),1);
% fwhm = ones(length(rng),1);
% for i=rng
%     I1 = nData(i).fitStr.I1;
% %     gamma = nData(i).fitStr.gamma;
%     x0 = nData(i).fitStr.x01;
%     ftot = @(x)-I1*voigtIkedaCarpenter_ord(x,[0,140,0,gamma,6.6e-3,0.05,x0]);%fnfit = -1*[fit function] so that the maximum becomes a minimum
%     xM(i) = fminbnd(ftot,hc-.2,hc+.2);% Identify position of the maximum on interval [hc-0.2 hc+0.2]
%     M = -ftot(xM(i));% compute value of the maximum
%     fd = @(x)abs(ftot(x)+M/2);
%     xhm1 = fminbnd(fd,xM(i)-0.2,xM(i));% identify x value for which function equals M/2 on interval [xM-0.2 xM]
%     xhm2 = fminbnd(fd,xM(i),xM(i)+0.2);% same on interval [xM xM+0.2]
%     fwhm(i) = xhm2 - xhm1;% FWHM of big peak
% end

%% Extract splitting between peaks as the distance between peak maxima
xM1 = ones(length(fullHrange),1); xM2 = ones(length(fullHrange),1); 
splitting = zeros(length(fullHrange),1);
hM1 = repmat(xM1,1); hM2 = repmat(xM2,1); 
Imax1 = repmat(xM1,1); Imax2 = repmat(xM2,1); 
for i=fullHrange
    label = sprintf(fmt1,nData(i).temp,field(i),Nprms);
%     I1 = nData(i).(fitStr).I1; 
    I2 = nData(i).(fitStr).I2;
%         gamma1 = nData(i).fitStr.gamma1; gamma2 = nData(i).fitStr.gamma2;
    x01 = nData(i).(fitStr).x01; x02 = nData(i).(fitStr).x02;
    f1 = @(x)-(I2*0.46*voigtIkedaCarpenter_ord(x,[0,140,0,0,0.05,6.6e-3,x01]));
    f2 = @(x)-(I2*voigtIkedaCarpenter_ord(x,[0,140,0,0,0.05,6.6e-3,x02]));
%fnfit = -1*[fit function] so that the maximum becomes a minimum
    xM1(i) = fminbnd(f1,peakCenter(1)-.1,peakCenter(1)+.1);% Identify position of the maximum of fit of peak 1 on interval [hc-0.2 hc+0.2]
    xM2(i) = fminbnd(f2,peakCenter(2)-.1,peakCenter(2)+.1);% Same for peak 2
    splitting(i) = -(xM2(i) - xM1(i))/hcenter;
    [~,idxM1] = min(abs((nData(i).hh0-xM1(i))));% extract index of max of peak 1 in dataset
    hM1(i) = nData(i).hh0(idxM1);% position of peak 1 in dataset
    Imax1(i) = nData(i).I(idxM1);% max of peak 1
    [~,idxM2] = min(abs((nData(i).hh0-xM2(i))));% same for peak 2
    hM2(i) = nData(i).hh0(idxM2);%
    Imax2(i) = nData(i).I(idxM2);%
end

% Estimation of error bars using error bars on peak position from fit
warning('Check the row number of the table in structure Stbl before running this section!')
j = 1;% check the row number of the table in structure Stbl before running!
spltRE = Stbl(j).(fitStr).x01_RelErr + Stbl(j).(fitStr).x02_RelErr;
wghts = min(spltRE)./spltRE;

%% Plot splitting
figure
% plot(field,splitting,'.')
errorbar(field,splitting,spltRE(fullHrange),'.','MarkerSize',12);
ylim([0 6e-3])

%% Fit splitting
[xData, yData, weights] = prepareCurveData( field, splitting, wghts );
HMaxFit = 0.74;% maximum value of field for which to include data for fit
% Set up fittype and options.
ft = fittype( 'delta0*sqrt(1-(H/Hc)^2)', 'independent', 'H', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Domain', [0 HMaxFit] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.8 6e-3];
opts.Weights = weights;
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

%% Print fit parameter values with error bars
cval = coeffvalues(fitresult);% extract fit parameter values; 
% cval(1) is the value of H_c, cval(2) that of the orthorhombic distortion
cft=confint(fitresult);% extract confidence intervals from fit
strSplit = [sprintf('$\\frac{(a-b)}{a_0}|_{H=0}$ = %.2f(%.0f)',...
    cval(2)*1e3,(cval(2)-cft(1,2))*1e5) '$\cdot 10^{-3}$'];% print out value of splitting at zero field, with error bars
strHMaxFit = sprintf('Max value of H/Hc for fit of orthorhombic distortion: %g',...
    HMaxFit/cval(1));

%% Calculate effective temperature
Tc0 = 2.15;% critical temperature at zero field (better to use the actual value of Tc from Cp data?)
max_splitting = 5.84e-3;% maximum value of splitting, from Segmuller et al. 1974
x = cval(2)/max_splitting;% reduced splitting calculated from fit
dx = (cval(2)-cft(1,2))/max_splitting;% error bar on reduced splitting calculated from fit
Teff = Tc0*x/atanh(x);% effective temperature calculated from data
% 2.2*x/atanh(x) is simply the result of inverting the self-consistent
% equation of pseudospin vs temperature
dTeff = Tc0*dx*abs(atanh(x)-x/(1-x^2))/(atanh(x)^2);% effective temperature calculated from data
% calculated by differentiating the expression of Teff wrt x
sTeff = sprintf('T = %.2f(%.0f)K',Teff,dTeff*1e2);% Effective temperature with error bars
sTdr = sprintf('$T_{DR}$=%.1fK',nData(1).temp);

%% Estimate critical field at the effective temperature in the absence of demagnetizing factor
% Function defining the critical field
t = Teff / Tc0;% reduced temperature
fnh = @(h) h - t*atanh(h);% when this function goes to zero, the value of h is the reduced critical field
figure
fplot(fnh,[0 1])

%% Compute critical field
h_c = fzero(fnh,[1e-3 1-1e-3]);
H_c = h_c*Hc_0;% value of the critical field at the effective temperature of
% the neutrons data, in the absence of demagnetizing factor
sprintf("Critical field at T/Tc0=%.2f: Hc(T)/Hc0 = %.2f",t,h_c)

%% Create dataset for 2D color plot
nData2 = nData;% create an independent copy of nData
for i=length(nData2):-1:1
    if size(nData2(i).I)~=2200% if none of the dimensions of nData2(i).I is 2200
        nData2(i)=[];% simply remove this dataset from the structure
    end
end

%% Treatment of magnetic field to correct for demag
demag_correction = H_c/cval(1);% cval(1) is the value of critical field 
% that comes from the fit of the splitting extracted from the data at 0.94K
% and 0.969 is the value of Hc/Hc0 at the effective temperature of the measurement
% see results of the below code when using demag_correction = 1
field2 = extractfield(nData2,'field')*H_c/cval(1);%
strHc = sprintf('$H_c$ = %.3f(%.0f)T',H_c,...
    (cval(1)-cft(1,1))*H_c/cval(1)*1e3);% print out value of critical field, with error bars

%% Compute conversion factor from position of maxima to orthorhombic lattice parameter 
at = 7.0426;% in-plane lattice parameter in the tetragonal phase
rstolp = sqrt(2)*at/hcenter;% conversion factor from reciprocal space units
% to units of the in-plane lattice parameter in the orthorhombic phase

%% Prepare matrix of peak maxima for exportation
Mmax = horzcat(field'/cval(1),xM1*rstolp,xM2*rstolp,splitting,spltRE);% matrix resulting from horizontal concatenation of 
% field normalized to the critical field, 
% orthorhombic lattice parameters a and b (calculated from positions of both peaks), 
% relative splitting of the peaks (equals orthorhombic distortion) and 
% relative error of the latter
% Note: horzcat corresponds to stacking columns side by side

%% Write matrix Mmax to a tab delimited file
filename = '2019-06-28_TmVO4_p6K_ENS_peak-max_vs_field.txt';
% Header
hdr1={'H/Hc','a_o','b_o','Orthorhombic distortion','Relative error on orthorh. distortion'};% First line header: quantities names
hdr2={'(no unit)','angstrom','angstrom','(no unit)','(no unit)'};% Second line header: units
hdr3={sTeff,strHMaxFit};% Third line header: comments
fmt1 = repmat('%s\t ', 1, length(hdr1)); fmt1(end:end+1) = '\n';% String formatting for header lines 1 and 2
fmt3 = repmat('%s\t ', 1, length(hdr3)); fmt3(end:end+1) = '\n';% String formatting for header line 3
fid = fopen(filename, 'wt');
fprintf(fid, fmt1, hdr1{:});% header 1
fprintf(fid, fmt1, hdr2{:});% header 2
fprintf(fid, fmt3, hdr3{:});% header 3
fclose(fid);
% Write matrix
% dlmwrite(filename,Mmax,'-append','Delimiter','\t')% Increase precision of
% exportation to full precision?

%% Prepare tight subplot
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m,n,p, [0 0], [0.1 0.1], [0.2 0.01]); 
if ~make_it_tight,  clear subplot;  end

%% Plot fit with data.
plotRng = field<HMaxFit;
xPlot = xData(plotRng)/cval(1);%xData(plotRng)*H_c/cval(1) to plot using calculated value of Hc(T)
yPlot = yData(plotRng);
spltPlot = spltRE(plotRng);
exclPlot = excludedPoints(plotRng);
Xfit = linspace(0,max(field)/cval(1),1000);%max(field)*H_c/cval(1) to plot using calculated value of Hc(T)
Yfit = cval(2).*sqrt(1-(Xfit).^2);% compute fit over a controlled number of points
% Xfit/(cval(1)*H_c/cval(1)) to plot using calculated value of Hc(T)
fig = figure; fig.Units = 'inches'; fig.Position(2:4)=[0.5 4 6];
ax2 = subplot(2,1,2);
pfit = plot(Xfit,Yfit*1e3,'r-'); hold on; grid on
pdat = plot(xPlot,yPlot*1e3,'.b','MarkerSize',18);
% pdat = errorbar(xPlot,yPlot*1e3,spltPlot,'.b','MarkerSize',10,'LineWidth',2);
% pexcl = plot(xPlot(exclPlot),yPlot(exclPlot),'xk','MarkerSize',6);
% legend([pdat,pexcl,pfit],'Splitting','Excluded','MF fit');
% legend([pdat,pfit],'Splitting','MF fit');
% title('TmVO$_4$ (8 8 0) peak splitting vs field');
xlabel('$H/H_c(T)$'); ylabel({'$\delta = 2 (a_o-b_o)/(a_o+b_o)$', '($\times 10^{-3}$)'});
xlim([0 1.1]); %ylim([0 6e-3]);
% ax = gca; 
ax2.TitleFontSizeMultiplier = 0.8; ax2.YTick = [0 2 4]; ax2.FontSize=16;
% ann21 = annotation('textbox',[0.15 0.3 0.2 0.1],'interpreter','latex',...
%     'String',{'Bragg peak at (8 8 0)' sTeff strSplit strHc},...
%     'FontSize',14,'LineStyle','-','EdgeColor','k','FitBoxToText','on',...
%     'LineWidth',2,'BackgroundColor',[1 1 1],'Color','k');% add annotation
ann22 = annotation('textbox',[ax2.Position(1)+.01 ax2.Position(2)+.01 0.2 0.1],'interpreter','latex',...
    'String',{'(b)'},'FontSize',ax2.FontSize,'LineStyle','-','EdgeColor','none',...
    'FitBoxToText','on','LineWidth',2,'BackgroundColor','none','Color','k');% add annotation
ann23 = annotation('textbox',[0.2 ax2.Position(2)+.05 0.2 0.1],'interpreter','latex',...
    'String',{'TmVO$_4$' sTeff},...
    'FontSize',14,'LineStyle','-','EdgeColor','k','FitBoxToText','on',...
    'LineWidth',1,'BackgroundColor',[1 1 1],'Color','k');% add annotation

%% Figure formatting
ann22.Position(1) = ax2.Position(1)+.01; ann22.Position(2) = ax2.Position(2)+.01;
ann23.Position(1) = ax2.Position(1)+.1/diff(ax2.XLim)*ax2.Position(3);% position textbox at x=0.2 in x axis units
ann23.Position(2) = ax2.Position(2)+.1;
ax2.LineWidth = 1;

%% Data formatting for curve fitting tool analysis
i=1;
hh1 = nData2(i).hh0;
I1dat = nData2(i).I;
dI = nData2(i).dI;

%% Color plot of intensity in H-(hh0) 2D-map
% Prepare plot
[X,Y] = meshgrid(field2,rstolp*hh1);
Ifull = cell2mat( arrayfun(@(c) c.I', nData2(1:length(nData2)).', 'Uniform', 0) );% intensity data combined in one big matrix
Ift = Ifull';

%% Select data to plot
ybounds = [9.85 10.05];
aocenter = sqrt(2)*at;
Ysel = Y>ybounds(1) & Y<ybounds(2)+.01;

%% Plot
% fsplt = figure;
ax1 = subplot(2,1,1); 
ax1.LineWidth = 1; ax1.Layer = 'top';% show ticks on top of plot
n = 300;
% contourf(Hg./5100,Tg,-d1Cpg,n,'EdgeColor','none');
sp = contourf(X(Ysel(:,1),:)/H_c,Y(Ysel(:,1),:),Ift(Ysel(:,1),:),n,'EdgeColor','None');
% sp = surf(X(Ysel(:,1),:),Y(Ysel(:,1),:),-Ifull(Ysel(:,1),:),'EdgeColor','None');
% I am plotting -Ifull instead of Ifull in order to be able to plot the
% position of peak maxima on top; see 'Matlab_debugging.pptx' for more info
ylim([ybounds(1) ybounds(2)]); %xlim([0 max(field2)]);
ax1.XTickLabel={}; ax1.YTick=[9.85 9.95 10.05]; ax1.FontSize=16;
colormap(jet);
cb = colorbar; cb.Location = 'north'; cb.AxisLocation = 'out';% cb.Direction = 'reverse';% reverse colorbar for negative intensities
cb.Position(1) = ax1.Position(1); cb.Position(3) = ax1.Position(3);% expand colorbar to make it match the width of the plot
cb.Position(2) = ax1.Position(2)+ax1.Position(4)+.01;% move colorbar just above figure
cb.Ticks = 0:2e6:6e6;% reduce number of colorbar ticks
cbtl = CbTixAbsolute(ax1,cb,0.0,0);% Reposition colorbar tick labels
for i=1:length(cbtl); cbtl(i).FontSize=ax1.FontSize; end% adjust fontsize of tick labels to match the rest of the figure
cbl = cb.Label; cbl.String = '(a.u.)';%
cbl.Interpreter = 'latex'; cb.TickLabelInterpreter = 'latex';
cbl.Position = [7e6 1 0]; cbl.Rotation = 0;% horizontal colorbar label
ylabel('$a_o,b_o$ (\AA)');% xlabel('$H/H_c(T)$');
if exist('sTeff','var'); Tstr = sTeff; elseif exist('sTdr','var'); Tstr = sTdr; end

% ann11 = annotation('textbox',[0.7 0.85 0.2 0.1],'interpreter','latex',...
%     'String',{Tstr},'FontSize',14,'FontName','Arial','LineStyle','-','EdgeColor','r',...
%     'FitBoxToText','on','LineWidth',1,'BackgroundColor',[1 1 1],'Color','k');% add annotation
ann12 = annotation('textbox',[ann22.Position(1) 0.51 0.2 0.1],'interpreter','latex',...
    'String',{'(a)'},'FontSize',ax1.FontSize,'LineStyle','-','EdgeColor','none',...
    'FitBoxToText','on','BackgroundColor','none','Color','w');% add annotation
% ann13 = annotation('textbox',[0.2 0.85 0.2 0.1],'interpreter','latex',...
%     'String',{'TmVO$_4$'},...
%     'FontSize',14,'LineStyle','-','EdgeColor','k','FitBoxToText','on',...
%     'LineWidth',1,'BackgroundColor',[1 1 1],'Color','k');% add annotation
caxis('auto');% auto rescale of color scale
grid off;
view(0,90);
hold on;
plot(field(plotRng)/cval(1),xM1(plotRng)*rstolp,'.k');% plot position of max of peak 1 vs H/Hc
plot(field(plotRng)/cval(1),xM2(plotRng)*rstolp,'.k');% same for peak 2
linkaxes([ax1,ax2],'x');

%% 
Iftplot = Ift(Ysel(:,1),:);
Xplot = X(Ysel(:,1),:)/H_c;
Yplot = Y(Ysel(:,1),:);

%% Write the intensity matrix to a tab delimited file
% dlmwrite('2019-06-28_TmVO4_p6K_ENS_I_vs_field_&_hh0.txt',Iftplot,'Delimiter','\t')
% Xplot and Yplot exported manually on 2019-06-28 by copying from Matlab 
% variables and pasting into Originlab worksheet














