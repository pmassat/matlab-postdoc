function TLFIM_plot(x,y,h,e)
%% Plot
figure;
ax = axes('LineStyleOrder',{'-',':','--','-.'});
hold(ax,'on');
for jh = 1:length(h)
    for je = 1:length(e)
        ax.ColorOrderIndex = jh;% Same line color for same value of h
        ax.LineStyleOrderIndex = je;% Same line style for same value of e
        plot(x,y(:,jh,je),'LineWidth',2,...
            'DisplayName',sprintf('%.2g, %.2g',h(jh),e(je)));
    end
end
legend('show');
lgd = ax.Legend;
title(lgd, '$h,e$');
