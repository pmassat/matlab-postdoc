%% Define physical constants
R = 8.314;% gas constant, in J/mol/K

%% Graphics properties
set(groot,'defaultTextFontsize',18);
set(groot,'defaultAxesFontsize',18);
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultAxesBox','on');
set(groot,'DefaultLineMarkerSize',18);
% set(0,'defaultAxesXLimSpec', 'tight')
set(0,'defaultAxesXMinorTick', 'on');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Ensure that the figure toolbar appears when opening a figure window
% This setting is only for MatlabR2018b
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
