cd 'C:/Users/Pierre/Desktop/Postdoc/Software/Matlab/Matlab_scripts/TFIM/Magnetization_TFIM/';

%% 
fig = figure; 
hold on; 
temps = [.1, 1.3, 1.8, 2.0, 2.14];% measurement temperatures, in Kelvin units, in figure 2 of Cooke et al. 1972
hmax = 3;
h  = linspace(0,hmax,hmax*100+1);
for jt = length(temps):-1:1
    T = temps(jt);
    t = T/2.15; 
    pm{jt} = plot(h, magnetization_TFIM(t,h),...
        'DisplayName',sprintf('%.2f', t));
end

Hc = 5000;% critical field, in Oersted
Hext = 4000:4000:8000;%
cmap = [0 .5 1; 0 0 1];

for jh=1:length(Hext)
    xhext = Hext(jh)/Hc;
    xline(xhext, 'LineWidth', 2, 'LineStyle', '--', 'Color', cmap(jh,:),...
        'HandleVisibility', 'off');
    text(xhext+.1,1.05,['  ' sprintf('%.2i Oe',Hext(jh))],...
        'Color', cmap(jh,:),'HorizontalAlignment','center');
    xhinmin = 0.63*xhext;
    xhinmax = 0.9*xhext;
%     xline(xhin, 'LineWidth', 2, 'LineStyle', '--', 'Color', cmap(2*jh,:), 'HandleVisibility', 'off');
%     area([xhinmin xhinmax], [1 1], 'FaceColor', cmap(2*jh,:),'EdgeColor', 'None',...
%         'FaceAlpha', 0.25, 'HandleVisibility', 'off');
%     text(mean([xhinmin xhinmax]), .25, '$H_{\mathrm{in}}$', 'HorizontalAlignment', 'center');
end

txtHext = text(.2,1.05,'$H_{\mathrm{ext}}=$');
% txttc = text(2.0,.7,'$T_c=2.15$ K');
% txtHunit = text(xhext*1.25,1.05,'Oe');
lgd = legend([pm{:}],'Location','east'); 
lgd.Title.String = '$T/T_{c,0}$';
% title('Magnetization of TmVO$_4$ vs magnetic field')
xlabel('$H/H_{c,0}$')
ylabel('$M/M_{\mathrm{sat}}$')
ylim([0 1.1])

% Plot MFD's at 4kOe and 8kOe for sample TmVO4-RF-E
% Run the first 100 lines of file
% 2017-07_TmVO4-RF-E_Analyze_Cp_under_field.m first, i.e. those lines that
% import and compute the relevant variables for the MFD's
rf_rng = [21:4*6:48];
yyaxis right

for jmfd = 1:length(rf_rng)
    i = rf_rng(jmfd);
    Trf = Smfd_RF(i).T_K;%
    Hext_rf = Smfd_RF(i).Hext_Oe%
    % p = plot(Smfd_RF(i).binCenters, Smfd_RF(i).hc, '.-', 'DisplayName', sprintf('%.2g, %.2g',T/Tc0rf,Hext/Hc));
    p = plot(Smfd_RF(i).binCenters, Smfd_RF(i).hc, '.-', 'Color', cmap(jmfd,:));
end
ylabel('Probability density')

lgd.String = lgd.String(1:5);
ax = gca; 
ax.YColor = [0 0 1];% aka 'blue'
ar = annotation('arrow', [.5 .7], [.25 .25], 'Color', 'blue');


%% Export figure
% cd 'C:\Users\Pierre\Desktop\Postdoc\Software\Matlab\Matlab_scripts\TFIM\Magnetization_TFIM'
% formatFigure 
printPNG([todaystr '_magnetization_TFIM_5values'])




