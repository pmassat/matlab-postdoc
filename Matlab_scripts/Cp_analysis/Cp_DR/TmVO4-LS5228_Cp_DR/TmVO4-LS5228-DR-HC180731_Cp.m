% Code was derived from 'YTmVO4_primary_analysis.m' on 2019-08-06
M = 283.87326;% Molar mass of each sample, in g/mol
m = 1e-3*0.38;% mass in g of the sample on which each dataset was measured

%% Analyze heat capacity from DR
% This routine is intended at anaLzing Cp data acquired with our DR used in 
% the Dynacool PPMS of the Lee lab
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2018-08_TmVO4-LS5228\'

%% Import data YTmVO4
Data5228HC18 = ImportCpDR('2018-07-31_TmVO4-LS5228-DR-HC180731.dat');
% Data0 = importCpSharedPPMS_('20180322_TmVO4-LS5228-MP3-Plt-HC1803_Cp.dat');

%% Store data into variables
% H=[DATA.FieldOersted; DATA(2).FieldOersted];
% T=[DATA.SampleTempKelvin; DATA(2).SampleTempKelvin];
% Cp=[DATA.SampHCJmoleK; DATA(2).SampHCJmoleK];
H5228HC18=[Data5228HC18.FieldOersted];
T5228HC18=[Data5228HC18.SampleTempKelvin];
Cp5228HC18=[Data5228HC18.SampHCJK];
CpErr5228HC18=[Data5228HC18.SampHCErrJK];

whichPoints = isfinite(H5228HC18) & isfinite(T5228HC18) & isfinite(Cp5228HC18);
H5228HC18 = H5228HC18(whichPoints);
T5228HC18 = T5228HC18(whichPoints);
Cp5228HC18 = Cp5228HC18(whichPoints)*1e-6*M/m;% molar heat capacity, converted from uJ/K to J/mol/K
CpErr5228HC18 = CpErr5228HC18(whichPoints)*1e-6*M/m;

[uh5228,~,X] = unique(round(H5228HC18(~isnan(H5228HC18))));
% split = accumarray(X,1:length(H(~isnan(H))),[],@(r){Data(r,:)});
luh5228 = length(uh5228);

%% Store data into cell array (outdated)
% isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
% % function to test if a column exists in a table
% for i=1:luh5228% for doped samples measured in DR
%     if isTableCol(split{i},'SampleTempKelvin')
%         split{i}.Properties.VariableNames{'SampleTempKelvin'} = 'T';% rename the temperature column
%     end
%     if isTableCol(split{i},'FieldOersted')
%         split{i}.Properties.VariableNames{'FieldOersted'} = 'H';% rename the magnetic field column
%     end
%     if isTableCol(split{i},'SampHCJK')% for samples measured in DR
%         split{i}.Properties.VariableNames{'SampHCJK'} = 'Cp';% rename the heat capacity column
%     elseif isTableCol(split{i},'SampHCmJmoleK')% for samples measured in shared PPMS
%         split{i}.Properties.VariableNames{'SampHCmJmoleK'} = 'Cp';
%     elseif isTableCol(split{i},'SampHCJmoleK')% for samples measured in Fisher He4 PPMS
%         split{i}.Properties.VariableNames{'SampHCJmoleK'} = 'Cp';
%     end
%     if isTableCol(split{i},'SampHCErrJmoleK')% for samples measured in DR
%         split{i}.Properties.VariableNames{'SampHCErrJmoleK'} = 'CpErr';
%     elseif isTableCol(split{i},'SampHCErrJK')% for samples measured in DR
%         split{i}.Properties.VariableNames{'SampHCErrJK'} = 'CpErr';        
%     end
% end
% 
% %% Remove NaN rows
% for i=1:luh5228
%     split{i}(any(isnan(split{i}.T), 2), :) = [];% Remove rows where T is NaN
% end
% 
% %% Keep only data under zero magnetic field
% for i=1%for all datasets
%     split{i} = split{i}(round(split{i}.H,-1)==0,:);% keep only data at zero field
% end
% %% Compute molar heat capacity
% for i=1:luh5228
%     split{i}.Cpmol = split{i}.Cp *1e-6*M/m;% molar heat capacity, in J/mol/K
%     split{i}.CpmolErr = split{i}.CpErr *1e-6*M/m;% molar heat capacity, in J/mol/K
%     % starting from a heat capacity measured in microJoules per Kelvin, as is measured in the DR
%     % Cpmol is calculated per mole of Tm3+ ions, hence the (1-dpg) factor in the denominator
% end
% 
% %% Sort each dataset by increasing value of temperature
% srtd = repmat(split,1);
% for i=1:luh5228
%     srtd{i} = sortrows(split{i},{'T'});
% %     [split{i}.T,wo] = sort(split{i}.T);
% %     split{i}.Cp = split{i}.Cp(wo);
% end
% srtd = srtd';
% 

%% Parameters for plotting heat capacity
xlblTemp = '$T$ (K)';
ylblCp = 'C$_p$ (J$\cdot$mol$^{-1}\cdot$K$^{-1}$)';
ttlCp = 'Heat capacity of TmVO$_4$';

%% Separate data according to value of magnetic field
clear sepH5228HC18

for i = 1:luh5228
    wp = abs(H5228HC18-uh5228(i))<50;
    sepH5228HC18(i).H = H5228HC18(wp);
    sepH5228HC18(i).T = T5228HC18(wp);
    sepH5228HC18(i).Cp = Cp5228HC18(wp);
    sepH5228HC18(i).CpErr = CpErr5228HC18(wp);

    [sepH5228HC18(i).T,wo] = sort(sepH5228HC18(i).T);
    sepH5228HC18(i).H = sepH5228HC18(i).H(wo);
    sepH5228HC18(i).Cp = sepH5228HC18(i).Cp(wo);
    sepH5228HC18(i).CpErr = sepH5228HC18(i).CpErr(wo);
end

%% Compute average of data points taken
clear avg5228HC18
R = 8.314;% gas constant in J/mol/K
Tp = 26;% Temperature scale of phonons contribution to Cp in TmVO4, in K; see 'TmVO4_Cp_phonons.m'
for i = luh5228:-1:1
    avg5228HC18(i) = averageCpwithH2(6e-3, sepH5228HC18(i).T,...
        sepH5228HC18(i).Cp, sepH5228HC18(i).CpErr, sepH5228HC18(i).H);
end
for i = luh5228:-1:1
    avg5228HC18(i).CpFull = avg5228HC18(i).Cp;
    avg5228HC18(i).Cp = avg5228HC18(i).CpFull - R*(avg5228HC18(i).T/Tp).^3;% electronic contribution to Cp, after subtracting phonons contribution
    avg5228HC18(i).Cpr = avg5228HC18(i).Cp/R;
    avg5228HC18(i).CprErr = avg5228HC18(i).CpFullErr/R;
end

%% Plot averaged data 
figure
plot(avg5228HC18(i).T,avg5228HC18(i).CpFull,'.','DisplayName','$C_p^{\mathrm{full}}$')
hold on
plot(avg5228HC18(i).T,R*(avg5228HC18(i).T/Tp).^3,'.','DisplayName','$C_p^{\mathrm{phonons}}$')
plot(avg5228HC18(i).T,avg5228HC18(i).Cp,'.','DisplayName','$C_p^{\mathrm{full}}-C_p^{\mathrm{phonons}}$')
legend('show')
title('TmVO$_4$ heat capacity')
xlabel(xlblTemp); ylabel('$C_p$ (J/K/mol)');

%% Compute Cp for Gaussian distribution of fields
clear Cpnum
Tc5228 = 2.193;% (2.192, 2.194) Value of transition temperature obtained from fit with Cp_LFIM function in Curve Fitting Tool including data poitns up to T = 1.191 K
e5228 =  6.5e-4;% (0.0005145, 0.0007925) Value of longitudinal field obtained the same way as Tc5228
Hc = 5000;
rescaling = .9;
rngNum = 1:luh5228;
for i=rngNum
Cpnum(i).h = uh5228(i)/Hc*rescaling;
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
[~,~,Cpnum(i).single_h_w_e] = FSCp_TLFIM(Cpnum(i).t_single_h,Cpnum(i).h,e5228);
Cpnum(i).single_h_w_e = Cpnum(i).single_h_w_e';
end
toc

%% Plot averaged data at each field separately
figure; hold on
rngAvg = 1:luh5228;
clr = lines(luh5228);
eb = cell(size(rngAvg));
for i=rngAvg
%     fp = fplot(@(t) Cp_TFIM(t/Tc_5228HC18,fields5228HC18(i)/Hc),[0 3],'LineWidth',2);
    fp = plot(Cpnum(i).t_single_h(2:end-1)*Tc5228,Cpnum(i).single_h_w_e,'DisplayName','MF');
    clr(i,:) = get(fp,'Color');
end
for i=rngAvg
    eb{rngAvg==i} = errorbar(avg5228HC18(i).T,avg5228HC18(i).Cpr,avg5228HC18(i).CprErr,...
        '.','MarkerSize',18,'DisplayName',sprintf('%d',uh5228(i)),...
        'Color',clr(i,:),'LineWidth',2);
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
% cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2018-08_TmVO4-LS5228\TmVO4_TFIM_article1_fig1'
% formatFigure;
% printPNG([todaystr mfd_ID '_Cp_vs_T_H'])
% printPDF([todaystr '_TmVO4-LS5228-DR-HC1807_Cp+single-h-fits_e=6p5e-4']);
printPNG([todaystr '_Cp_vs_T_H=0_fits_MF_MF+strain']);



%% Prepare plot of theoretical heat capacity curve 
% Tc = 2.19;%  Value of transition temperature obtained from fit with Cp_LFIM_NF function in Curve Fitting Tool
Tc5228 = 2.185;%  Value of transition temperature obtained from fit with Cp_LFIM_CW function in Curve Fitting Tool
s =  5e-4;%  Value of longitudinal field obtained from fit with Cp_LFIM_NF and Cp_LFIM_CW and Cp_LFIM_MEC_cst_stress function in Curve Fitting Tool
T_CW = 2.09;% Curie-Weiss temperature obtained from fit with function Cp_LFIM_CW 
A_CW = 0.0136;% amplitude of CW divergence obtained from fit with function Cp_LFIM_CW 
T0 = 1.6;% Divergence temperature obtained from fit with function Cp_LFIM_MEC_cst_stress
A = 5e-3;% Amplitude of "high"-temperature divergence obtained from fit with function Cp_LFIM_MEC_cst_stress
tvec = linspace(1e-3,4,300);% temperature vector for computation of Cp

%% Compute pure mean-field heat capacity
Cptheo0 = Cp_TFIM(tvec/Tc5228,0);

%% Compute theoretical heat capacity in LFIM
% Cptheo = Cp_LFIM(tvec/Tc,s);
% Cptheo = Cp_LFIM_CW(tvec/Tc,s,A,T0/Tc);
Cptheo = Cp_LFIM_MEC_cst_stress(tvec/Tc5228,s,A,T0/Tc5228);

%% Subtract Cptheo from experimental data to analyze the residual Cp
Cptheodat = Cp_LFIM(avg5228HC18(i).T/Tc5228,s);% theoretical Cp computed at same temperatures as data
avg5228HC18(i).Cpres = avg5228HC18(i).Cp - Cptheodat*R;% Residual Cp after subtraction of the fit

%% Prepare fit of Cp vs Temperature 
% Use these variables in Curve fitting tool
i=1;
fitT = avg5228HC18(i).T;
fitCp = avg5228HC18(i).Cpr;
fitCpErr = avg5228HC18(i).CprErr;
fitwghts = 1./fitCpErr;
% fitCpres = avg5228HC18(i).Cpres/R;

%% Prepare tight subplot for combining plots in single figure
% For figure 1 of paper 1
sf = 1;%5./6;% scaling factor of figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.0 0.0], [0.1 0.02], [0.13 0.02]*1/sf);
figure; formatFigure([6 7]*sf);
ax = gca; ax.delete

%% Plot averaged data + fit including external stress
figure; 
% coefficients of fit of residual Cp, obtained in Curve Fitting Tool with Adjusted R-square of 0.9927
c1 =     0.02607;  %(0.01489, 0.03725)
c2 =    0.004759;  %(0.003322, 0.006196)
% ax2 = subplot(5,1,[3:5]); 
p0 = plot(tvec,Cptheo0,'Color',[.5 .5 .5],'DisplayName','Mean-field');
hold on
p1 = plot(tvec,Cptheo,'-r','DisplayName',['MF + strain']);
% fp = plot(avg5228HC18(i).T,Cptheodat,'-r','DisplayName',sprintf('MF: e=%.2g',e));
pdat = errorbar(avg5228HC18(i).T,avg5228HC18(i).Cp/R,avg5228HC18(i).CpFullErr/R,'.b',...
    'MarkerSize',18,'LineWidth',1,'DisplayName','Experiment');
xlabel(xlblTemp); ylabel('$C_p^{4f}/R$');
ylim([0 1.55])
% annttl = annotation('textbox',[0.2 0.45 0.17 0.1],'interpreter','latex',...
%     'String',{'TmVO$_{4}$'}, 'LineStyle','-','EdgeColor','black',...
%     'BackgroundColor','none','Color','k','LineWidth',1,'FitBoxToText','on');% add annotation
% annttl.FontSize = ax.XAxis.Label.FontSize;
grid on
% title(ttlCp);
lgd = legend([pdat, p0, p1], 'Location', 'northwest');

%% Plot averaged data + CW fit 
figure; 
Tfit = 2.22; T66 = 1.8;
% coefficients of fit of residual Cp, obtained when fitting in Curve Fitting Tool
T0 =       2.091;%  (2.079, 2.103); fit a/(T-T0),  Adjusted R-square of 0.9992
a =     0.02932;%  (0.02726, 0.03138); fit a/(T-T0),  Adjusted R-square of 0.9992
a66 =     0.03713;%  (0.03533, 0.03892); fit a66/(T-T66),  Adjusted R-square of 0.9926
% fp = fplot(@(t)a/(t-T0),[0 4],'-r','LineWidth',1.5,...
%     'DisplayName',[sprintf('%.1g/(T-%.2g)',a,T0)]);
p1 = fplot(@(t)a66/(t-T66)^2,[0 4],'-r','LineWidth',1.5,...
    'DisplayName',[sprintf('%.2g/(T-%.2g)',a66,T66) '$^2$']);
hold on
errorbar(fitT(fitT<Tfit),fitCpres(fitT<Tfit),fitCpErr(fitT<Tfit),'xk',...
    'MarkerSize',9,'LineWidth',1,'DisplayName','Excluded')
errorbar(fitT(fitT>Tfit),fitCpres(fitT>Tfit),fitCpErr(fitT>Tfit),'.b',...
    'MarkerSize',18,'LineWidth',1,'DisplayName','Data')
% errorbar(fitT,fitCpres,fitCpErr,'.b','MarkerSize',18,'LineWidth',1,'DisplayName','Data')
xlabel(xlblTemp); ylabel('$\delta C_p/R$');
xlim([0 4]); ylim([-.05 .25]);
annttl = annotation('textbox',[0.2 0.8 0.1 0.1],'interpreter','latex',...
    'String',{'TmVO$_{4}$'}, 'LineStyle','-','EdgeColor','none',...
    'BackgroundColor','none','Color','k');% add annotation
% annttl.FontSize = ax.XAxis.Label.FontSize;
grid on
% title(ttlCp);
lgd = legend('show');






%% % % Plot figure % % % 
%% Compute splitting of GS doublet vs temperature
% For figure 1 of paper 1
sz = zeros(size(tvec));
for j=1:length(tvec)
sz(j) = 2.95/2*OP_TFIM(tvec(j)/Tc5228,0,s);%2.95 cm-1 is the total splitting of the GS doublet, see e.g. Melcher1976, fig.2(a)
end

%% Plot splitting vs T
% Top of figure 1 of paper 1
ax1 = subplot(5,1,[1:2]); 
plot(tvec,sz)
hold on
plot(tvec,-sz,'Color',lines(1))
% ylim([-2 3]);
% ylabel('$E$ (cm$^{-1}$)')
% arr = annotation('doublearrow',[.2 .2],[0.76 .925]);% add annotation
% anngap = annotation('textbox',[0.91 0.79 0.2 0.1],'interpreter','latex',...
%     'String',{'$E_{g}$'}, 'LineStyle','-','EdgeColor','none',...
%     'BackgroundColor','none','Color','k');% add annotation

%% Figure formatting when combining plots
ax1.YLabel.Position(1) = ax2.YLabel.Position(1);
anngap.FontSize = ax1.XAxis.Label.FontSize;
ax1.XTickLabel={};
ax1.YTick = [];
ax2.XLabel.Position(2) = -.15;

%% Print figure to PNG file
% cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_summary_figures\TmVO4_GS-splitting+multipoles+Cp'
% formatFigure;
printSVG([todaystr '_TmVO4-LS5228-DR-HC180731_Cp+fit'])
% printPNG('2019-11-19_TmVO4_MFIM_Cp')
% printPNG('2019-11-19_TmVO4_Cp_data+MF')
% printPNG('2019-11-19_TmVO4_Cp_data+stress')
% printPNG('2019-11-20_TmVO4_CpRes_data+C66fit')

%% Export averaged Cp data to text file
% savepath = 'C:\Users\Pierre\Desktop\Postdoc\YTmVO4\YTmVO4_HeatCapacity\YTmVO4_Cp_anaLsis\2018-10-17_YTmVO4_averaged_Cp\';
% for i=L-3:L
%     dpg(i);
%     expstr(i) = fullfile(savepath, ['AnaLzeCpDRdoping_2018-11-30_YTmVO4_x=',num2str(dpg(i)),'.txt']);
%     fid = fopen(expstr(i), 'wt');
%     fprintf(fid, '%s\n', ['Exported from AnaLzeCpDRdoping_2018-10-10_YTmVO4_VTmAsO4.mlx on ',date]);  % header
%     fclose(fid);
%     dlmwrite(expstr(i),cat(2,avg5228HC18(i).T',avg5228HC18(i).Cp',avg5228HC18(i).stdCp'),'-append','delimiter','\t')
% end
% 





%% % % % % % % % Failed attempts at fitting tail of Cp data at T>Tc % % % % % % % % % % 

%% Fit at low temperatures
maxTfit = 2.23;
[fitresult, gof] = fitCpLFIM(fitT,fitCp,fitwghts,maxTfit);

%% Plot residual Cp and power law fit
% Power law a*T^b fit parameters for residual Cp above T_D
% a = 22.58;%  (-6.487, 51.65); b = -6.062;%  (-7.608, -4.516);% When fitting all data points above T_D
a = 0.9996;%  (0.203, 1.796)
b = -3.035;%  (-3.795, -2.275);% When fitting only the highest three data points
figure
errorbar(avg5228HC18(i).T,avg5228HC18(i).Cpres/R,avg5228HC18(i).CpFullErr/R,'.b',...
    'MarkerSize',18,'LineWidth',1,'DisplayName','Residual $C_p$')
hold on
% fpr = plot(tvec(2:end-1),Cpmagres','-r','DisplayName',['Mag. dip.' newline '${\Delta e}/T_c=$ ' sprintf('%.2g',d0md)]);
fpwr = fplot(@(t) a*t^b, [0 4],'-g','LineWidth',2,...
    'DisplayName',[sprintf('$T^{%.2g}$ fit',b)]);
xlabel(xlblTemp); ylabel('$\Delta C_p/R$');
grid on
title({'Residual $C_p$ data +' 'theoretical correction from magnetic dipoles'});% title(ttlCp);
lgd = legend('show');
ylim([0 .2])

%% Cp LFIM with arbitrary value of longitudinal field
e2 = 1e-2;
Cpt2 = Cp_LFIM(tvec/Tc5228,e2);

%% Compute Cp correction due to magnetic dipole interactions
C_0 = 0.088;% 1/T^2 coefficient, see paper
d0md = 0.01;
szcp = sz(2:end-1);
% tadr = linspace(3e-3,2,500);
deltaCpdr = (1-szcp'*2/2.95).*C_0/Tc5228.*Cp_magDipRandFields(tvec'/Tc5228,d0md);

%% Cp in random strains
d0dr = .3;
Tcrs = 2.3;
Cpma = Cp_full_random_strains(d0dr,6*s,tvec/Tcrs);

%% Plot averaged data + fit including random strains
figure; 
% ax2 = subplot(5,1,[3:5]); 
p1 = plot(tvec,Cptheo,'-r','DisplayName',sprintf('MF: e=%.2g',s));
% fp = plot(avg5228HC18(i).T,Cptheodat,'-r','DisplayName',sprintf('MF: e=%.2g',e));
hold on
fpr = plot(tvec,Cpma,'-g','DisplayName',['Rand. strains:' newline '${\Delta e}/T_c$=' sprintf('%.2g',d0dr)]);
errorbar(avg5228HC18(i).T,avg5228HC18(i).Cp/R,avg5228HC18(i).CpFullErr/R,'.b',...
    'MarkerSize',18,'LineWidth',1,'DisplayName','Experiment')
xlabel(xlblTemp); ylabel('$C_p/R$');
ylim([0 1.55])
% annttl = annotation('textbox',[0.2 0.75 0.1 0.1],'interpreter','latex',...
%     'String',{'TmVO$_{4}$'}, 'LineStyle','-','EdgeColor','none',...
%     'BackgroundColor','none','Color','k');% add annotation
% annttl.FontSize = ax.XAxis.Label.FontSize;
grid on
title(sprintf('$T_c=$ %.2g K',Tcrs));
lgd = legend('show');

%% Plot averaged data + theoretical correction from magnetic dipole interactions
figure; 
p1 = plot(tvec,Cptheo,'-r','DisplayName',sprintf('MF: $e=$ %.2g',s));
hold on
fpr = plot(tvec(2:end-1),deltaCpdr','-g','DisplayName',['Mag. dip.:' newline '${\Delta e}/T_c$=' sprintf('%.2g',d0md)]);
fpt = plot(tvec(2:end-1),Cptheo(2:end-1)+deltaCpdr','-m','DisplayName','MF + mag. dip.');
errorbar(avg5228HC18(i).T,avg5228HC18(i).Cp/R,avg5228HC18(i).CpFullErr/R,'.b',...
    'MarkerSize',18,'LineWidth',1,'DisplayName','Experiment')
xlabel(xlblTemp); ylabel('$C_p/R$');
ylim([0 1.55])
grid on
title(['TmVO$_4$ heat capacity: mean-field fit' newline '+ correction from magnetic dipoles']);% title(ttlCp);
lgd = legend('show');

%% Mag dip correction to Cp for comparison with Cpres
Cpmagres = C_0/Tc5228.*Cp_magDipRandFields(tvec'/Tc5228,d0md);

