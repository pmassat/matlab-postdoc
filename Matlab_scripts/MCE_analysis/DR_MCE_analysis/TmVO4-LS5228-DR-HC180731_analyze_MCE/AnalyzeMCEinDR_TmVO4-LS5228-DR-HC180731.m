%% Analysis of MCE data generated by the DR option
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2018-08_TmVO4-LS5228';
load('DR_thermometer_TvsHR_fit.mat');
% 'DR_thermometer_TvsHR_fit.mat' is an sfit variable called 'stf', exported from file 'DR_platform_therm_cal.mlx'
FullData=ImportMCEinDRmatrix('2018-08-01_MCE_TmVO4-LS5228-DR-HC180731.dat');%

%% 
% For some reason, Matlab does not want to import the MCE data in a table 
% directly as numerics, it wants to import it as text, then convert it to numerics, 
% which is not time efficient. The only way to import directly as numerics is 
% to use NumericMatrix as type in the import wizard. The drawback of this method 
% is that one has to indicate manually the number of each data column by counting 
% from the first column to be imported, which is 'Time Stamp' (see raw data file; 
% the 'Comment' column is not imported). However, this will hopefully not be a 
% problem on the long run, as all the data files generated by the DR option of 
% the PPMS are formatted the same way.
%% Parameters
Hmax = 10000;% high magnetic fields boundary for the plots
Tc180731 = 2.20;% from peak in dCp/dT in 'TmVO4_LS5228_201808_AnalyzeCpDR.mlx'

%% R-H MCE traces
Data = repmat(FullData,1);% make a copy of the data
Data(Data(:,20)<0,:)=[];% remove rows with negative bridge 2 resistance, which does not make any sense
t = Data(:,1);% Time Stamp in seconds is column #1
Tmce = Data(:,3);% in K; Temperature of the platform thermometer is column #3 in the data file,
Hmce = Data(:,4);% in Oe
excCurrent = Data(:,9);% Bridge 2 excitation current, in uA
Rb2 = Data(:,20);% Bridge 2 resistance, in Ohms
figure
plot(Hmce,Rb2,'.-')
% ylim([3000,9000])
xlim([0,Hmax])
xlabel('Field (Oe)')
ylabel('Platform resistance ($\Omega$)')
title('R-H MCE traces')

%% T-H MCE traces
% Convert Bridge 2 resistance into temperature by using the sfit generated from 
% the file 'DR_platform_therm_cal.mlx' and stored in 'DR_thermometer_TvsHR_fit.mat'
% 
% Then plot platform thermometer temperature and bridge 2 resistance temperature 
% together to check that the latter shows the MCE effect but not the former.

Tb2 = stf(Hmce,Rb2);% stf is the sfit variable created when loading 'DR_thermometer_TvsHR_fit.mat'
figure
plot(t,Tmce,'DisplayName','Thermometer')
hold on
plot(t,Tb2,'DisplayName','Bridge 2')
hold off
% ylim([0,2])
legend('show')
xlabel('Time stamp (s)')
ylabel('Temperature (K)')

%% Compute derivatives of relevant quantities
sigma = 3;% width of the smoothing gaussian, in number of data points
[Tb2smth, d1tb2s, d2tb2s] = gConvolve(Tb2,sigma);% smooth platform bridge 2 temperature data & derivative
[Hsmth, d1hs, d2hs] = gConvolve(Hmce,sigma);% smooth field data and derivative
[tmsth, d1ts, d2ts] = gConvolve(t,sigma);% smooth timestamp and derivative; 

%% Plot the actual T-H MCE traces
figure
plot(Hmce,Tb2,'.-')
% ylim([0,2])
xlim([0,Hmax])
xlabel('Field (Oe)')
ylabel('Temperature (K)')
title('T-H MCE traces')

%% Filter data according to bridge 2 excitation current
uec = unique(excCurrent)% values of excitation current without repetition
FilterEC = false(length(excCurrent),length(uec));% initialize logical array
for i = 1:length(uec)% for each value of excitation current
    figure
    FilterEC(:,i) = excCurrent==uec(i);% store the corresponding logical vector as a separate column
    plot(Hmce(FilterEC(:,i)),Tb2(FilterEC(:,i)),'.-',...
        'DisplayName',strcat('I=',num2str(round(uec(i),2,'significant')),'uA'))
    % plot data for each excitation current separately
    legend('show')
    xlim([0,Hmax])
    xlabel('Field (Oe)')
    ylabel('Temperature (K)')
end
%% Filter data according to magnetic field sweep rate
% swprt = diff(Hmce)./diff(t);% magnetic field sweep rate;
% Computation of sweeprate using diff is deprecated, as it removes a data point
sweeprate = d1hs./d1ts;
usr = unique(round(sweeprate,-1));% values of sweep rate rounded to nearest tenth without repetition
for ii = length(usr):-1:1% for each sweeprate value
    if length(sweeprate(abs(sweeprate-usr(ii))<5))<10% if there are less than 10 data points with this sweeprate
        usr(ii)=[];% remove it
    end
end
usr(usr==0)=[];%remove sweeprate = 0
lusr = length(usr);
FilterSR = false(length(sweeprate)+1,lusr);% initialize logical array
for j = 1:lusr% for each value of sweep rate
%     figure
    FilterSR(2:end,j) = abs(sweeprate-usr(j))<5;
    % gather together sweeprates that are less than 5 Oe/s away from each reference sweeprate
%     plot(H(FilterSR(:,j)),Tb2(FilterSR(:,j)),'.','DisplayName',strcat(num2str(usr(j)),'Oe/s'))
    % plot data for each magnetic field sweep rate separately
%     legend('show')
%     xlim([0,Hmax])
%     xlabel('Field (Oe)')
%     ylabel('Temperature (K)')
end

%% Split temperature and field arrays based on sweeprate
Tsrup = cell(lusr/2,1);% initialize cell arrays
Hsrup = cell(lusr/2,1);
Tsrdn = cell(lusr/2,1);
Hsrdn = cell(lusr/2,1);
for isr = 1:lusr/2% for each value of sweep rate
    Tsrup{isr} = Tmce(FilterSR(:,lusr-isr+1));% cell array of temperatures for upsweep data
    Hsrup{isr} = Hmce(FilterSR(:,lusr-isr+1));% cell array of fields for upsweep data
    Tsrdn{isr} = Tmce(FilterSR(:,isr));% cell array of temperatures for downsweep data
    Hsrdn{isr} = Hmce(FilterSR(:,isr));% cell array of fields for downsweep data
end

%% Plot upsweep and downsweep data with same sweeprate together
for isr = 3%1:lusr/2% for each absolute value of sweep rate
    figure; hold on
    plot(Hsrup{isr},Tb2(FilterSR(:,lusr-isr+1)),'.','DisplayName',sprintf('+%g Oe/s',usr(lusr-isr+1)));
	plot(Hsrdn{isr},Tb2(FilterSR(:,isr)),'.','DisplayName',sprintf('%g Oe/s',usr(isr)));
    hold off
    % plot data for each magnetic field sweep rate separately
    legend('show','Location','best')
    xlim([3200 6800]); ylim([1 1.8]);
    xticks(3500:1000:6500);
    xlabel('$H$ (Oe)')
    ylabel('$T$ (K)')
    grid on
end

%% Export data to text file
savepath = pwd;
todaystr = datestr(datetime('today'),29);
fname = [todaystr '_TmVO4-LS5228-DR-HC180731_MCE_'];
headerstr = sprintf('This data is restricted to sweeprates of +-%i Oe/s',abs(usr(isr)));
M(1).data = cat(2,Hsrup{isr},Tb2(FilterSR(:,lusr-isr+1)));
M(1).id = 'upsweep';
M(2).data = cat(2,Hsrdn{isr},Tb2(FilterSR(:,isr)));
M(2).id = 'downsweep';
for i=1:2
    expstr = fullfile(savepath, [fname M(i).id '.dat']);
    fid = fopen(expstr, 'wt');
    fprintf(fid, '%s\n', ['Exported from AnalyzeMCEinDR_TmVO4-LS5228-DR-HC180731 on ' todaystr]);  % header
    fprintf(fid, '%s\n', [headerstr]);  % header
    fprintf(fid, '%s\t%s\n','H (Oe)','T (K)');  % header
    fclose(fid);
    dlmwrite(expstr,M(i).data,'-append','delimiter','\t')
end

%% Figure exportation
% xlim([0 10000]); ylim([0.55 0.83]);
% printPNG('2019-05-20_TmVO4-LS5228-DR-HC180731_MCE_20Oeps_p6K-p7K-p8K');
% printPDF('2019-07-22_TmVO4-LS5228-DR-HC180731_MCE_20Oeps_1p1K-1p7K');

%% Select full usable dataset
% See labnotes for tables listing usable data
useusr = usr(abs(usr)<30); lu2 = length(useusr);
FullFilterSRup = false(size(sweeprate));
FullFilterSRdown = false(size(sweeprate));
for k=1%k=1 corresponds to 20Oe/s only
    FullFilterSRup = FullFilterSRup | abs(sweeprate-useusr(lu2+1-k))<5;% Data at sweeprate of +40 Oe/s are not usable
    FullFilterSRdown = FullFilterSRdown | abs(sweeprate-useusr(k))<5;% Data at sweeprate of -40 Oe/s are not usable
end
FullFilterEC = round(excCurrent,1)==0.1 | round(excCurrent,1)==0.4 |...
    round(excCurrent,1)==1.0;% Usable excitation currents are .1 uA, .4 uA and 1.0 uA
FullFilterT = round(Tb2,1)>=1.;% Usable temperature range is .6 K and above
FullFilterH = Hmce>0;% FullFilterH = Hmce>3300 & Hmce<6000;

FullFilterUp = FullFilterT & FullFilterH & FullFilterEC & FullFilterSRup;
FullFilterDown = FullFilterT & FullFilterH & FullFilterEC & FullFilterSRdown;

%% Plot full usable dataset
% Use this section to plot the MCE traces on top of the T-H dCp/dT colormap
% from '2017-07_TmVO4-RF-E_Analyze_Cp_under_field.m'
figure% comment this line out to combine with the 2D colormap
pup = plot(Hmce(FullFilterUp),Tb2(FullFilterUp),'.');
% Divide by Tc180731 to normalize temperature and by 5500 for magnetic field
hold on
pdown = plot(Hmce(FullFilterDown),Tb2(FullFilterDown),'.');
% xlabel('$H / H_c(T=0)$'); ylabel('$T / T_D(H=0)$');
xlabel('$H$ (Oe)'); ylabel('$T$ (K)');
% title('TmVO4-LS5228-DR-HC180731 MCE full usable dataset');
lgd = legend([pup,pdown],'Upsweep','Downsweep');

%% Print the phase diagram combining dCp/dT + MCE traces
% printPNG('2019-05-21_TmVO4-LS5228-DR-HC180731_MCE_+_2017-07-20_TmVO4_dCp-dT')

%% Rescale temperature of MCE data by factor (1+(1-T)/10)
%% Note about the interpretation of data
% The MCE observed here is different from what is reported in Kohama et al. 

%% Plot the first derivative of temperature change with respect to magnetic field
mceder = d1tb2s./d1hs*10^3.*exp(2-Tmce);% derivative of the temperature change
% the exponential term allows the low temperature data to be multiplied by a higher factor than the high T data
for isr = 3%1:lusr/2% for each value of sweep rate
    figure
    plot(Hsrup{isr},Tsrup{isr}+mceder(FilterSR(:,lusr-isr+1)),'.','DisplayName',strcat(num2str(usr(lusr-isr+1)),'Oe/s'))
    hold on
    plot(Hsrdn{isr},Tsrdn{isr}+mceder(FilterSR(:,isr)),'.','DisplayName',strcat(num2str(usr(isr)),'Oe/s'))
    hold off
    % plot data for each magnetic field sweep rate separately
    legend('show')
    xlim([0,Hmax])
    ylim([0 max(Tmce)])
    xlabel('Field (Oe)')
    ylabel('T+d$\Delta$T')
    title(['$\sigma = $' sprintf('%i',sigma)])
end

%% Same using differentiation
% In order to not get confused with shifting of datasets due to reduction
% of number of data points when using diff, I prefer the above convolution method
mcediff = diff(Tb2smth)*10^2;
for isr = 1:lusr/2% for each value of sweep rate
    figure
    plot(Hsrup{isr},Tsrup{isr}+mcediff(FilterSR(:,lusr-isr+1)),'.','DisplayName',strcat(num2str(usr(lusr-isr+1)),'Oe/s'))
    hold on
    plot(Hsrdn{isr},Tsrdn{isr}+mcediff(FilterSR(:,isr)),'.','DisplayName',strcat(num2str(usr(isr)),'Oe/s'))
    hold off
    % plot data for each magnetic field sweep rate separately
    legend('show')
    xlim([0,Hmax])
    ylim([0 max(Tmce)])
    xlabel('Field (Oe)')
    ylabel('T+d$\Delta$T')
end

%% Plot the second derivative of temperature change with respect to magnetic field
% sigma = 3;% width of the smoothing gaussian, in number of data points
% [~, d1tb2s, d2tb2s] = gConvolve(Tb2,sigma);% smooth platform bridge 2 temperature data & derivative
% [~, d1hs, d2hs] = gConvolve(Hmce,sigma);% smooth field data and derivative
mceder2 = d2tb2s./d1hs.^2*10^5.*exp(2-Tmce);% derivative of the temperature change
% the exponential term allows the low temperature data to be multiplied by a higher factor than the high T data
for isr = 3%1:lusr/2% for each value of sweep rate
    figure
    plot(Hsrup{isr},Tsrup{isr}+mceder2(FilterSR(:,lusr-isr+1)),'.','DisplayName',sprintf('+%g Oe/s',usr(lusr-isr+1)))
    hold on
    plot(Hsrdn{isr},Tsrdn{isr}+mceder2(FilterSR(:,isr)),'.','DisplayName',sprintf('%g Oe/s',usr(isr)))
    hold off
    % plot data for each magnetic field sweep rate separately
    legend('show')
%     xlim([0,Hmax]); ylim([0 max(Tmce)])
    xlim([3200,6500]); ylim([1.0 1.8]);
    xticks(3500:1000:6500);
    xlabel('$H$ (Oe)'); ylabel('$T$ (K) $+\frac{\partial^2\Delta T}{\partial^2 H}$ (K/Oe$^{2})$');
%     title(['$\sigma = $' sprintf('%i',sigma)])
end

%% Split mceder2 array based on sweeprate
d2mceus = cell(lusr/2,1);% initialize cell arrays
d2mceds = cell(lusr/2,1);
for isr = 1:lusr/2% for each value of sweep rate
    d2mceus{isr} = mceder2(FilterSR(:,lusr-isr+1));% cell array of 2nd derivative of MCE data for field swept up
    d2mceds{isr} = mceder2(FilterSR(:,isr));% cell array of 2nd derivative of MCE data for field swept down
end

%% Identify the maximum of the second derivative at each temperature
clear utus utds
isr = 3;% index of sweep rate; 3 corresponds to +-20 Oe/s
utus = unique(round(Tsrup{isr},1));% unique values of temperatures for the selected sweeprate 
Hfilterus = Hsrup{isr} > 3300 & Hsrup{isr} < 6500;
for itmce = 1:length(utus)
    Tfilter = round(Tsrup{isr},1)==utus(itmce);
    maxd2mce(itmce).T = utus(itmce);
    maxd2mce(itmce).upsweep(2) = max(d2mceus{isr}(Tfilter & Hfilterus));
    maxd2mce(itmce).upsweep(1) = Hsrup{isr}(Tfilter & Hfilterus & d2mceus{isr}==maxd2mce(itmce).upsweep(2));
end

%% Same for downsweeps
utds = unique(round(Tsrdn{isr},1));%
if utds ~= utus
    warning("The values of temperature measured for up- and downsweeps are not the same.")
end
Hfilterds = Hsrdn{isr} > 3300 & Hsrdn{isr} < 6500;
for itmce = 1:length(utds)
    Tfilter = round(Tsrdn{isr},1)==utds(itmce);
    maxd2mce(itmce).T = utus(itmce);
    maxd2mce(itmce).downsweep(2) = max(d2mceds{isr}(Tfilter & Hfilterds));
    maxd2mce(itmce).downsweep(1) = Hsrdn{isr}(Tfilter & Hfilterds & d2mceds{isr}==maxd2mce(itmce).downsweep(2));
end


%% 
% Next steps:
% 
% # Combine excitation current and sweeprate filters to plot figures ==> identify 
% relevant traces; perhaps remove sweeprate of +-40Oe/s, where traces do not seem 
% to contain useful information
% # Plot derivative of relevant traces
% # Identify max of derivative => transition temperature (with error bar)
% # Plot phase diagram
%% Below this point: code imported from AnalyzeSR830MCEPierre >> to be adpated
% ref_sweep = 0.6e-3;%Reference sweep rate of magnetic field (in Tesla/s) used 
%to discriminate between high and low sweep rates
FilterUp = diff(Data(1).H)>0;
FilterDown = diff(Data(1).H)<0;
Hup = Data(1).H(FilterUp);
ResUp = Data(1).ResPlat(FilterUp);
Hdown = Data(1).H(FilterDown);
ResDown = Data(1).ResPlat(FilterDown);
figure
plot(Hup,ResUp,'.-')
hold on
plot(Hdown,ResDown,'.-')
legend('Up','Down')
%%
sgm=10;%number of points for smoothing
newH=gConvolve(SR830MCE(3).H,sgm);%smooth the field data
newT=gConvolve(1./BetaModel(SR830MCE(3).H,SR830MCE(3).ResPlat),sgm);%smooth the temperature data
NewD=diff(newT);%compute the derivative of the smoothed temperature
%%
figure
plot(newH(30:(end-30)),NewD(30:(end-29))*50+SR830MCE(3).TPuck(30:(end-30)),'-')%plot dT/dH versus H
%ylim([0.3 1.5])
xlabel('Field (Oe)')
ylabel('dT/dH + T_{puck}')
title('MCE effect')
%%
figure;
yyaxis left;
plot(diff(SR830MCE(3).H),'DisplayName','Raw');
hold on;
plot(diff(newH),'-g','DisplayName','Smoothed');
xlabel('Measurement time t (s)')
ylabel('dH/dt (T/s)')
title('Sweep rate and puck temperature')
legend('show')

yyaxis right;
plot(SR830MCE(3).TPuck,'DisplayName','Puck temperature')
ylabel('T_{puck} (K)')
%%
ref_sweep = 0.6e-3;%Reference sweep rate of magnetic field (in Tesla/s) used 
%to discriminate between high and low sweep rates
FilterFastUp = diff(newH)>ref_sweep;
FilterSlowUp = diff(newH)>0 &  diff(newH)<ref_sweep;
FilterFastDown = diff(newH)<-ref_sweep;
FilterSlowDown = diff(newH)<0 &  diff(newH)>-ref_sweep;
%%
figure
SelectOnlyFastUp.newH=newH(FilterFastUp);
SelectOnlyFastUp.newH=SelectOnlyFastUp.newH(30:(end-30));
SelectOnlyFastUp.NewD=NewD(FilterFastUp);
SelectOnlyFastUp.NewD=SelectOnlyFastUp.NewD(30:(end-30));
SelectOnlyFastUp.TPuck=SR830MCE(3).TPuck(FilterFastUp);
SelectOnlyFastUp.TPuck=SelectOnlyFastUp.TPuck(30:(end-30));
plot(SelectOnlyFastUp.newH,SelectOnlyFastUp.NewD*50+SelectOnlyFastUp.TPuck,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
maintitle = dataname + date + newline;
title(maintitle + "Sweeping up in field at 8 Oe/s (0.4K) and 10 Oe/s (1.0K, 1.1K, 1.2K, 1.3K)")
%%
figure
SelectOnlySlowUp.newH=newH(FilterSlowUp);
SelectOnlySlowUp.NewD=NewD(FilterSlowUp);
SelectOnlySlowUp.TPuck=SR830MCE(3).TPuck(FilterSlowUp);
%pts_ctff = 50;%Number of points to cut off when plotting the data in order 
%to avoid the Gaussian aberrations at the edges of the plots
%SelectOnlySlowUp.newH=SelectOnlySlowUp.newH(pts_ctff:(end-pts_ctff));
%SelectOnlySlowUp.NewD=SelectOnlySlowUp.NewD(pts_ctff:(end-pts_ctff));
%SelectOnlySlowUp.TPuck=SelectOnlySlowUp.TPuck(pts_ctff:(end-pts_ctff));
plot(SelectOnlySlowUp.newH,SelectOnlySlowUp.NewD*100+SelectOnlySlowUp.TPuck,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH*100')
title(maintitle + "Sweeping up in field at 4 Oe/s")
%%
figure
SelectOnlyFastDown.newH=newH(FilterFastDown);
SelectOnlyFastDown.NewD=NewD(FilterFastDown);
SelectOnlyFastDown.TPuck=SR830MCE(3).TPuck(FilterFastDown);
plot(SelectOnlyFastDown.newH,SelectOnlyFastDown.NewD*50+SelectOnlyFastDown.TPuck,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title(maintitle + "Sweeping down in field at 8 Oe/s (0.4K) and 10 Oe/s (1.0K, 1.1K, 1.2K, 1.3K)")
%%
figure
SelectOnlySlowDown.newH=newH(FilterSlowDown);
SelectOnlySlowDown.NewD=NewD(FilterSlowDown);
SelectOnlySlowDown.TPuck=SR830MCE(3).TPuck(FilterSlowDown);
plot(SelectOnlySlowDown.newH,SelectOnlySlowDown.NewD*100+SelectOnlySlowDown.TPuck,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title(maintitle + "Sweeping down in field at 4 Oe/s")
%%
figure
plot(SelectOnlyFastDown.newH,SelectOnlyFastDown.NewD*50+SelectOnlyFastDown.TPuck,'.','DisplayName','10 Oe/s down')
hold on
plot(SelectOnlySlowDown.newH,SelectOnlySlowDown.NewD*100+SelectOnlySlowDown.TPuck,'.','DisplayName','4 Oe/s down')
plot(SelectOnlySlowUp.newH,SelectOnlySlowUp.NewD*100+SelectOnlySlowUp.TPuck,'.','DisplayName','4 Oe/s up')

plot(SelectOnlyFastUp.newH,SelectOnlyFastUp.NewD*50+SelectOnlyFastUp.TPuck,'.','DisplayName','10 Oe/s up')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title(maintitle + "Comparing sweep rates")
legend('show')