%% Change to relevant directory
cd C:\Users\Pierre\Desktop\Postdoc\Software\Matlab\Matlab_data_analysis\TFIM\OP_TFIM

%% Plot the order parameter of the TFIM 
% The order parameter in TFIM with non zero longitudinal field must show a rounding 
% because the longitudinal field has the same symmetry as the OP (by
% definition of a longitudinal field)
% e = 0.01;
T = linspace(0,3,1000);
den = repmat(T,1);
Y = repmat(T,1);
f = @(h) h/atanh(h);
fig1 = figure; ax1 = subplot(1,1,1); hold(ax1,'on');
fig2 = figure; ax2 = subplot(1,1,1); hold(ax2,'on');
for h=.0%[0 0.6 0.9]
    counter = 0;% reset counter used to determine which color to use for the plot
    disp(counter)
    for e = [0, 0.01]
        for j=1:length(T)
        r = OP_TFIM(T(j),h,e)/T(j);
        den(j) = sech(r)^2./T(j);
        Y(j) = 1./(1 - den(j));
        end
        if counter>0
            c1 = get(p1,'Color');
            plot(ax1,T,den,'LineWidth',2,'Color',c1,...
                'DisplayName',sprintf('$h=%.2g$, $e=%.2g$',h,e));
            plot(ax2,T,Y,'LineWidth',2,'Color',c1,...
                'DisplayName',sprintf('$h=%.2g$, $e=%.2g$',h,e));
        else
            p1 = plot(ax1,T,den,'LineStyle',':','LineWidth',2,...
                'DisplayName',sprintf('$h=%.2g$, $e=%.2g$',h,e));
            p2 = plot(ax2,T,Y,'LineStyle',':','LineWidth',2,...
                'DisplayName',sprintf('$h=%.2g$, $e=%.2g$',h,e));
        end
        counter = counter+1;
    end
%         line([f(h),f(h)],ylim,'LineStyle','--','LineWidth',1,'Color','k',...
%         'DisplayName',sprintf('$T_c(h=%.2f)$',h));
end
l = line(ax1,T,ones(1,length(T)),'Color','k','Linewidth',1);
legend(ax1,'show'); legend(ax2,'show')
xlabel(ax1,'$T/T_{c}(h=0)$'); ylabel(ax1,'Denominator');
xlabel(ax2,'$T/T_{c}(h=0)$'); ylabel(ax2,'Y');
% xlim([0 1.2]); 
% ylim([-.1 1.1]);
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
printPDF('2019-07-02_OP_TFIM_vs_t_h_e')
