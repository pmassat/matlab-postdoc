% Code was derived from 'YTmVO4_primary_analysis.m' on 2019-08-06

%% Analyze heat capacity from DR
% This routine is intended at anaLzing Cp data acquired with our DR used in 
% the Dynacool PPMS of the Lee lab
cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2018-08_TmVO4-LS5228\TmVO4_TFIM_article1_fig1'

%% To do: 
% * Error bars:
% *     For dT data, use standard deviation

%% Import data YTmVO4
Data = ImportCpDR('2018-07-31_TmVO4-LS5228-DR-HC180731.dat');
% Data0 = importCpSharedPPMS_('20180322_TmVO4-LS5228-MP3-Plt-HC1803_Cp.dat');
%% Concatenate them in a cell array
split = {Data};
L = length(split);

%% Rename variables
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
% function to test if a column exists in a table
for i=1:L% for doped samples measured in DR
    if isTableCol(split{i},'SampleTempKelvin')
        split{i}.Properties.VariableNames{'SampleTempKelvin'} = 'T';% rename the temperature column
    end
    if isTableCol(split{i},'FieldOersted')
        split{i}.Properties.VariableNames{'FieldOersted'} = 'H';% rename the magnetic field column
    end
    if isTableCol(split{i},'SampHCJK')% for samples measured in DR
        split{i}.Properties.VariableNames{'SampHCJK'} = 'Cp';% rename the heat capacity column
    elseif isTableCol(split{i},'SampHCmJmoleK')% for samples measured in shared PPMS
        split{i}.Properties.VariableNames{'SampHCmJmoleK'} = 'Cp';
    elseif isTableCol(split{i},'SampHCJmoleK')% for samples measured in Fisher He4 PPMS
        split{i}.Properties.VariableNames{'SampHCJmoleK'} = 'Cp';
    end
    if isTableCol(split{i},'SampHCErrJmoleK')% for samples measured in DR
        split{i}.Properties.VariableNames{'SampHCErrJmoleK'} = 'CpErr';
    elseif isTableCol(split{i},'SampHCErrJK')% for samples measured in DR
        split{i}.Properties.VariableNames{'SampHCErrJK'} = 'CpErr';        
    end
end

%% Remove NaN rows
for i=1:L
    split{i}(any(isnan(split{i}.T), 2), :) = [];% Remove rows where T is NaN
end

%% Keep only data under zero magnetic field
for i=1:L%for all datasets
    split{i} = split{i}(round(split{i}.H,-1)==0,:);% keep only data at zero field
end

%% Parameters for plotting heat capacity
xlblTemp = 'Temperature (K)';
ylblCp = 'C$_p$ (J$\cdot$mol$^{-1}\cdot$K$^{-1}$)';
ttlCpY = 'Heat capacity of Tm$_{1-x}$Y$_x$VO$_4$';

%% Compute molar heat capacity
M = [283.87326% Molar mass of each sample, in g/mol
275.3902538
271.869006
269.4681552
266.3470492
266.2670208
264.6664536
258.6643266];
m = 1e-3*[0.38% mass in g of the sample on which each dataset was measured
0.13
0.24
0.045
0.25
0.59
0.76
0.14];
dpg = [0% Y content for each dataset, after interpolation of uprobe 
% results; see 'YTmVO4_phase_diagram.m'
0.106
0.164
0.196
0.219
0.237
0.264
0.315];
% split{1}.Cpmol = split{1}.Cp *1e-3;% molar heat capacity, in J/mol/K,
% starting from a heat capacity in mJ/mol/K as is measured in the shared PPMS
for i=1:L
    split{i}.Cpmol = split{i}.Cp *1e-6*M(i)/(m(i)*(1-dpg(i)));% molar heat capacity, in J/mol/K
    split{i}.CpmolErr = split{i}.CpErr *1e-6*M(i)/(m(i)*(1-dpg(i)));% molar heat capacity, in J/mol/K
    % starting from a heat capacity measured in microJoules per Kelvin, as is measured in the DR
    % Cpmol is calculated per mole of Tm3+ ions, hence the (1-dpg) factor in the denominator
end

%% Sort each dataset by increasing value of temperature
srtd = repmat(split,1);
for i=1:L
    srtd{i} = sortrows(split{i},{'T'});
%     [split{i}.T,wo] = sort(split{i}.T);
%     split{i}.Cp = split{i}.Cp(wo);
end
srtd = srtd';

%% Compute average of data points taken
for i = 1:L
    avgData(i) = averageCp(6e-3,srtd{i}.T,srtd{i}.Cpmol,srtd{i}.CpmolErr);
end

%% Compute splitting of GS doublet vs temperature
tvec = linspace(1e-3,4,2000);
sz = zeros(size(tvec));
for j=1:length(tvec)
sz(j) = 2.95/2*OP_TFIM(tvec(j)/Tc,0,e);%2.95 cm-1 is the total splitting of the GS doublet, see e.g. Melcher1976, fig.2(a)
end

%% Prepare tight subplot
subplot = @(m,n,p) subtightplot (m, n, p, [0.0 0.0], [0.1 0.02], [0.13 0.02]);
figure; formatFigure([6 7]);
ax = gca; ax.delete

%% Plot splitting vs T
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

%% Prepare plot of theoretical curve 
Tc = 2.193;%  (2.192, 2.194) value of transition temperature obtained from fit using Curve Fitting Tool
e =  0.000643;%  (0.000583, 0.000703) value of longitudinal field obtained from fit using Curve Fitting Tool
Cptheo = Cp_LFIM(tvec/Tc,e);

%% Plot averaged data
R = 8.314;
% figure; 
ax2 = subplot(5,1,[3:5]); 
fp = plot(tvec,Cptheo,'-r');
hold on
errorbar(avgData(i).T,avgData(i).Cp/R,avgData(i).CpFullErr/R,'.b',...
    'MarkerSize',18,'LineWidth',1)
xlabel(xlblTemp); ylabel('$C_p/R$');
ylim([0 1.55])
% annttl = annotation('textbox',[0.2 0.75 0.1 0.1],'interpreter','latex',...
%     'String',{'TmVO$_{4}$'}, 'LineStyle','-','EdgeColor','none',...
%     'BackgroundColor','none','Color','k');% add annotation
% annttl.FontSize = ax.XAxis.Label.FontSize;
grid on
% title(ttlCpY);
% lgd = legend('show');

%% Figure formatting when combining plots
ax1.YLabel.Position(1) = ax2.YLabel.Position(1);
anngap.FontSize = ax1.XAxis.Label.FontSize;
ax1.XTickLabel={};
ax1.YTick = [];
ax2.XLabel.Position(2) = -.15;

%% Prepare fit of Cp vs Temperature 
% Use these variables in Curve fitting tool
fitT = avgData(i).T;
fitCp = avgData(i).Cp/R;
fitCpErr = avgData(i).CpFullErr/R;
fitwghts = 1./fitCpErr;

%%
maxTfit = 2.23;
[fitresult, gof] = fitCpLFIM(fitT,fitCp,fitwghts,maxTfit);

%% Print figure to PNG file
% formatFigure
% printPNG('2019-05-28_YTmVO4_Cp_vs_T_0-p11-p15-p22')
% printSVG('2019-08-06_TmVO4-LS5228-DR-HC180731_Cp+fit')

%% Export averaged Cp data to text file
% savepath = 'C:\Users\Pierre\Desktop\Postdoc\YTmVO4\YTmVO4_HeatCapacity\YTmVO4_Cp_anaLsis\2018-10-17_YTmVO4_averaged_Cp\';
% for i=L-3:L
%     dpg(i);
%     expstr(i) = fullfile(savepath, ['AnaLzeCpDRdoping_2018-11-30_YTmVO4_x=',num2str(dpg(i)),'.txt']);
%     fid = fopen(expstr(i), 'wt');
%     fprintf(fid, '%s\n', ['Exported from AnaLzeCpDRdoping_2018-10-10_YTmVO4_VTmAsO4.mlx on ',date]);  % header
%     fclose(fid);
%     dlmwrite(expstr(i),cat(2,avgData(i).T',avgData(i).Cp',avgData(i).stdCp'),'-append','delimiter','\t')
% end
% 