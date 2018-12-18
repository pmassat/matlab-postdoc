%% Graphics properties
set(groot,'defaultTextFontsize',18);
set(groot,'defaultAxesFontsize',18);
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultAxesBox','on');
set(groot,'DefaultLineMarkerSize',18);
% set(0,'defaultAxesXLimSpec', 'tight')
set(0,'defaultAxesXMinorTick', 'on')

% Ensure that the figure toolbar appears when opening a figure window
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
