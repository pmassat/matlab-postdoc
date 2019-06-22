%% Plot the order parameter of the TFIM 
% e = 0.01;
T = linspace(0,3,1000);
Y = repmat(T,1);
f = @(h) h/atanh(h);
figure; hold on;
for h=0%[0 0.6 0.9]
    for e = 0.01%[0, 0.01]
        for j=1:length(T)
        Y(j) = OP_TFIM(T(j),h,e);
        end
    plot(T,Y,'LineWidth',2,'DisplayName',sprintf('$h=%.2f$, $e=%1.1d$',h,e))
    end
    line([f(h),f(h)],ylim,'LineStyle','--','DisplayName',sprintf('$T_c(h=%.2f)$',h));legend('show')
end
legend('show')
xlabel('$T/T_{c}(h=0)$'); ylabel('Order parameter');
% xlim([0 1.2]); 
ylim([-.1 1.1]);
title(sprintf('OP vs T in the TFIM'));
% fplot(@(t)1./CpTFIMdenominator(t,e),[1e-3 2]);
% xlabel('$T/T_c$'); ylabel('$C_p$ denominator');
% fplot(@(t)Cp_TFIM_offset_strain(t,e,h),[1e-3 2]);
% xlabel('$T/T_c$'); ylabel('$C_p$');
% ann01 = annotation('textbox',[0.6 0.75 0.2 0.1],'interpreter','latex',...
%     'String',{['$e=$ ' sprintf('%.2f',e)] ['$h=$ ' sprintf('%.2f',h)]},...
%     'LineStyle','-','EdgeColor','k',...
%     'FitBoxToText','on','LineWidth',1,'BackgroundColor','w','Color','k');% add annotation

%% Export figure
formatFigure;
% printPDF('2019-06-19_OP_TFIM')
