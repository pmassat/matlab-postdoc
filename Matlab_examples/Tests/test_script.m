% x = [0,.1,.2,.3];
% y = [1,.8,.5,0];
% figure
% area(x,y,'FaceColor',[0 .5 0]);% dark green [0 .5 0]; light green [0 1 0]
% xlim([])
figure;
hold on
e = 0.01;
for h=[0 0.6 0.8 0.9 0.99]
fplot(@(t)OP_TFIM(t,h,e),[0 1.5],'LineWidth',2,'DisplayName',sprintf('$h=%.2f$',h))
end
legend('show')
xlabel('$T/T_{c,0}$'); ylabel('Order parameter');
xlim([0 1.2]); ylim([-.1 1.1]);
title(sprintf('OP vs T in the TFIM w/ $e=%.3f$',e));
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
