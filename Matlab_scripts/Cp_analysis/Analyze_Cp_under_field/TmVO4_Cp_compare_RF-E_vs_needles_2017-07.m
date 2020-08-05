%% Requirements
% Before running this file, need to run:
% - AnalyzeTmVO4Cp2017-07-28.m to compute structure Smfd_ndl
% - 2017-07_TmVO4-RF-E_Analyze_Cp_under_field to compute structure Smfd_RF

%% Create structure containing both datasets
Scomp(1).Tmfd = Tmfd_RF;
Scomp(2).Tmfd = Tmfd_ndl;

%% Select distributions of fields at a given value of T and Hext for both samples
% phyisical parameters for data selection 
T_plt = [1, 1.1];% temperatures of the data to be plotted; 
% first element is for sample RF-E, second is for the needle-shaped samples
h_plt = [.98, .78];% Hext/Hc ratio of the data to be plotted, rounded to 2nd decimal

% logical arrays to select the relevant data
for comp_idx=2:-1:1
sel{comp_idx} = Scomp(comp_idx).Tmfd.T_K==T_plt(comp_idx) &...
    round(Scomp(comp_idx).Tmfd.Hext_Oe/Hc,2)==h_plt(comp_idx);%
end

%% plot distribution of magnetic fields
figure;
hold on
for comp_idx=2:-1:1
    x = Scomp(comp_idx).Tmfd.binCenters(sel{comp_idx},:);
    y = Scomp(comp_idx).Tmfd.hc(sel{comp_idx},:);
    sample = Scomp(comp_idx).Tmfd.label(sel{comp_idx});
    T = Scomp(comp_idx).Tmfd.T_K(sel{comp_idx});
    Hext = Scomp(comp_idx).Tmfd.Hext_Oe(sel{comp_idx});
    [nrows,~] = size(x);
    for row=1:nrows
    lgd_str = [sprintf('%s, ', sample(row,:)),...
        sprintf('%.2g K, ',T(row,:)), sprintf('%.2d Oe',Hext(row,:))];
    p = plot(x(row,:), y(row,:), '.-', 'DisplayName', lgd_str);
    end
end
lgd_mfd = legend('show'); lgd_mfd.Title.String = 'Sample, $T$, $H_{\mathrm{ext}}$';
title(['Distributions of magnetic fields in TmVO$_4$'])
xlabel('$H_{\mathrm{in}}/H_{\mathrm{c}}$')
ylabel('Normalized PDF')

%% Compare Cp data at similar fields 
Scomp(1).data = struct2table(avgRFData);
Scomp(2).data = struct2table(avgNdlData);
Scomp(1).sample = 'RF-E';
Scomp(2).sample = 'Needles';

%% Find index of Cp datasets corresponding to above field values
data_idx = zeros(size(h_plt));
for comp_idx=1:2
    data_idx(comp_idx) = find(round(Scomp(1).data.uh/Hc,2)==h_plt(comp_idx));
end

Tc = [Tc0, Tc_ndl];

%% plot two Cp datasets at equivalent value of fields
figure;
hold on
for comp_idx=1:2
    T = Scomp(comp_idx).data.T{data_idx(comp_idx)};
    y = Scomp(comp_idx).data.Cpelr{data_idx(comp_idx)};
    sample = Scomp(comp_idx).sample;
    Hext = Scomp(comp_idx).data.uh(data_idx(comp_idx));
    lgd_str = [sprintf('%s, ', sample), sprintf('%.2d Oe',Hext)];
    p = plot(T/Tc(comp_idx), y, '.-', 'DisplayName', lgd_str);
end
lgd = legend('show'); lgd.Title.String = 'Sample, $T$, $H_{\mathrm{ext}}$';
title(['Heat capacity of TmVO$_4$'])
xlabel('$T/T_c$')
ylabel('$C_p/R$')

%% Change to directory where figures should be stored
cd 'C:\Users\Pierre\Desktop\Postdoc\Software\Matlab\Matlab_scripts\Cp_analysis\Analyze_Cp_under_field\TmVO4_Cp_compare_RF-E_vs_needles_2017-07'

%% Export figure
% formatFigure;
printPNG([todaystr '_TmVO4_compare_RF_vs_needles_mfd']);
% printPDF(['2019-06-18_TmVO4-RF-E_fit_Schottky_' strrep(hrstr,'.','p') 'xHc']);













