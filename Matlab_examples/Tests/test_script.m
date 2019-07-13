% x = [0,.1,.2,.3];
% y = [1,.8,.5,0];
% figure
% area(x,y,'FaceColor',[0 .5 0]);% dark green [0 .5 0]; light green [0 1 0]
% xlim([])
fig1 = figure; ax1 = subplot(1,1,1); hold(ax1,'on');
fig2 = figure; ax2 = subplot(1,1,1);
for h=[0 0.6 0.8 0.9 0.99]
    plot(ax1,0:10,h*rand(1,11));
    line(ax2,0:5,6*h*ones(1,6));
end
legend(ax1,'show'); legend(ax2,'show')
xlabel(ax1,'$T/T_{c}(h=0)$'); ylabel(ax1,'Order parameter');

% legend('show')
% xlabel('$T/T_{c,0}$'); ylabel('Order parameter');
% xlim([0 1.2]); ylim([-.1 1.1]);
% title(sprintf('OP vs T in the TFIM w/ $e=%.3f$',e));

%% Export figure
formatFigure;
% printPDF('2019-06-19_OP_TFIM')
