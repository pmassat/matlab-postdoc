make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.0 0.0], [0.15 0.1], [0.15 0.01]);
if ~make_it_tight,  clear subplot;  end

%% Upper and Lower Subplots with Titles
figure
income = [15.52,15.51,15.50,15.56]*1e6;
outgo = [2.5,4.0,3.35,4.9];
ax1 = subplot(2,1,1); plot(income); ylabel('test');
title('Income');
ax2 = subplot(2,1,2); plot(outgo); ylabel('$\frac{C}{T}$')
title('Outgo'); xlabel('test');

%% Upper subplot with colorbar, lower subplot with 1D curve
figure; 
ax1 = subplot(2,1,1); contourf(1:4,1:4,meshgrid(income,outgo)); 
ylabel('color map'); cb = colorbar; cb.Location = 'north'; cb.AxisLocation = 'out';
ax1.XTickLabel = {};
ax2 = subplot(2,1,2); plot(outgo); ylabel('$\frac{crap}{Test}$')
 
%% Subplots in Quadrants
figure
subplot(2,2,1)
text(.5,.5,{'subplot(2,2,1)';'or subplot 221'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,2)
text(.5,.5,{'subplot(2,2,2)';'or subplot 222'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,3)
text(.5,.5,{'subplot(2,2,3)';'or subplot 223'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,4)
text(.5,.5,{'subplot(2,2,4)';'or subplot 224'},...
    'FontSize',14,'HorizontalAlignment','center')

%% Asymmetrical Subplots
figure
subplot(2,2,[1 3])
text(.5,.5,'subplot(2,2,[1 3])',...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,2)
text(.5,.5,'subplot(2,2,2)',...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,4)
text(.5,.5,'subplot(2,2,4)',...
    'FontSize',14,'HorizontalAlignment','center')

%%  
figure
subplot(2,2,1:2)
text(.5,.5,'subplot(2,2,1:2)',...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,3)
text(.5,.5,'subplot(2,2,3)',...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,4)
text(.5,.5,'subplot(2,2,4)',...
    'FontSize',14,'HorizontalAlignment','center')

%% Plotting Axes Over Subplots
figure
y = zeros(4,15);
for k = 1:4
    y(k,:) = rand(1,15);
    subplot(2, 2, k)
    plot(y(k,:));
end
hax = axes('Position', [.35, .35, .3, .3]);
bar(hax,y,'EdgeColor','none')
set(hax,'XTick',[])
