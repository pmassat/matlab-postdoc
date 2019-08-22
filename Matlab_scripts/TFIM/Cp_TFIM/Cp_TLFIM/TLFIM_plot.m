function TLFIM_plot(x,y,h,e)
%% Plot
figure;
ax = axes('LineStyleOrder',{'-',':','--','-.'});
hold(ax,'on');
for jh = 1:length(h)
    for je = 1%1:length(e)
        ax.ColorOrderIndex = jh;% Same line color for same value of h
        ax.LineStyleOrderIndex = je;% Same line style for same value of e
%         plot(x,y(:,jh,je),'LineWidth',2,'DisplayName',sprintf('%.2g,%.2g',h(jh),e(je)));% when looping over values of e for each value of h
        plot(x,y(:,jh),'LineWidth',2,'DisplayName',sprintf('%.2g, %.2g',h(jh),e(jh)));% when using only value of e for each value of h
    end
end
legend('show');
lgd = ax.Legend;
title(lgd, '$h,e$');
