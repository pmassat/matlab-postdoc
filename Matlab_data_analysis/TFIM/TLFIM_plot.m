function TLFIM_plot(x,y,h,e)
%% Plot
figure; hold on;
for jh = 1:length(h)
    counter = 0;% reset counter used to determine which color to use for the plot
    disp(counter)
    for je = 1:length(e)
        if counter>0
            c1 = get(p1,'Color');
            plot(x,y(:,jh,je),'LineWidth',2,'Color',c1,...
                'DisplayName',sprintf('$h=%.2g$, $e=%.2g$',h(jh),e(je)));
        else
            p1 = plot(x,y(:,jh,je),'LineStyle',':','LineWidth',2,...
                'DisplayName',sprintf('$h=%.2g$, $e=%.2g$',h(jh),e(je)));
        end
        counter = counter+1;
    end
end
legend('show')
