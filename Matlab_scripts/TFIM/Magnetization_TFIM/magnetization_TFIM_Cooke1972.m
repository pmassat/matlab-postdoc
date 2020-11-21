figure; 
hold on; 
% temps = [0.5, 1.34, 1.94];% measurement temperatures, in Kelvin units, in figure 2 of Cooke et al. 1972
temps = 2.15*[1e-2, .5, .9];% measurement temperatures, in Kelvin units, in figure 2 of Cooke et al. 1972
hmax = 3;
h  = linspace(0,hmax,hmax*100+1);
for T = temps
    t = T/2.15; 
    plot(h, magnetization_TFIM(t,h),...
        'DisplayName',sprintf('%.2g', t))
%         'DisplayName',sprintf('%.2g, %.2g', t, critical_field(t)))
end
lgd = legend('show', 'Location', 'southeast'); 
% lgd.Title.String = '$T/T_c$, $H_c/H_{c,0}$';
lgd.Title.String = '$t$';
% title('Magnetization of TmVO$_4$ vs magnetic field')
title('Transverse magnetization in the TFIM')
xlabel('$h$')
ylabel('$m$')
ylim([0 1.1])

%% Export figure
% cd 'C:\Users\Pierre\Desktop\Postdoc\Software\Matlab\Matlab_scripts\TFIM\Magnetization_TFIM'
formatFigure
% printPNG([todaystr '_magnetization_TFIM_Cooke1972'])