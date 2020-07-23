fig = figure; 
hold on; 
temps = [.1, 1.0, 1.75, 2.0, 2.14];% measurement temperatures, in Kelvin units, in figure 2 of Cooke et al. 1972
hmax = 3;
h  = linspace(0,hmax,hmax*100+1);
for T = temps
    t = T/2.15; 
    plot(h, magnetization_TFIM(t,h),...
        'DisplayName',sprintf(' %.3g, %.2f, %.2g', T, t, critical_field(t)))
end
Hc = 5100;% critical field, in Oersted
Hext = 4000:4000:8000;%
cmap = prism(6);
for jh=1:length(Hext)
    xhext = Hext(jh)/Hc;
    xline(xhext, 'LineWidth', 2, 'LineStyle', '--', 'Color', cmap(2*jh,:),...
        'HandleVisibility', 'off');
    text(xhext+.1,1.05,['  ' sprintf('%.2i Oe',Hext(jh))],...
        'Color', cmap(2*jh,:),'HorizontalAlignment','center');
    xhinmin = 0.63*xhext;
    xhinmax = 0.9*xhext;
%     xline(xhin, 'LineWidth', 2, 'LineStyle', '--', 'Color', cmap(2*jh,:), 'HandleVisibility', 'off');
    area([xhinmin xhinmax], [1 1], 'FaceColor', cmap(2*jh,:),'EdgeColor', 'None',...
        'FaceAlpha', 0.25, 'HandleVisibility', 'off');
    text(mean([xhinmin xhinmax]), .25, '$H_{\mathrm{in}}$', 'HorizontalAlignment', 'center');
end
txtHext = text(.2,1.05,'$H_{\mathrm{ext}}=$');
txttc = text(2.0,.7,'$T_c=2.15$ K');
% txtHunit = text(xhext*1.25,1.05,'Oe');
lgd = legend('show','Location','southeast'); 
lgd.Title.String = '$T$(K), $T/T_c$, $H_c/H_{c,0}$';
title('Magnetization of TmVO$_4$ vs magnetic field')
xlabel('$H/H_{c,0}$')
ylabel('$M/M_{\mathrm{sat}}$')
ylim([0 1.1])

%% Export figure
% cd 'C:\Users\Pierre\Desktop\Postdoc\Software\Matlab\Matlab_scripts\TFIM\Magnetization_TFIM'
% formatFigure 
printPNG([todaystr 'magnetization_TFIM_5values'])




