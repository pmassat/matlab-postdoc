%% Plot the phase diagram of (Tm,Y)VO4 vs Y content
cd 'C:\Users\Pierre\Desktop\Postdoc\Software\Matlab\Matlab_data_analysis\Cp_analysis\AnalyzeCpDoping\YTmVO4_Cp_Matlab';
%% Actual doping vs nominal doping
% Actual doping measured in microprobe

x_nom = [0.050; 0.102; 0.201; 0.301];% nominal x=Y/(Y+Tm) atomic ratio
x_up = [0.048; 0.106; 0.219; 0.315];% x value measured in microprobe
dx_up = [0.0048; 0.0049; 0.0109; 0.0119];% standard deviation of x measured in microprobe
figure
errorbar(x_nom,x_up,dx_up,'.-','Markersize',18)

%% Interpolate actual x_up values to estimate for doping levels not yet measured in microprobe
x_nm = [0.153; 0.181; 0.220; 0.248];% nominal doping levels of samples not measured in microprobe
% for which samples have not yet been measured in microprobe
x_dop = sort([x_nom; x_nm]);
x_intpl = interp1(x_nom,x_up,x_dop);% interpolate x_up at x values including x_nm
[xval,xind]=intersect(x_dop,x_nm);% extract values and indices of x_nm (interpolated x) in x_dop
x_intpl(xind)% show interpolated values of x at x_nm
dx_intpl = interp1(x_nom,dx_up,x_dop);% interpolate dx_up at x values including x_nm

%% Corrections of errorbars on certain doping levels
x_full = [0;x_intpl];
dxp = [0;dx_intpl];
dxm = [0;dx_intpl];
dxp(round(x_full,2)==0.20) = 0.05;
% nominal doping value is likely to differ significantly from the actual doping

%% Plot data + interpolation
figure
plot(x_nom,x_up,'o','DisplayName','From microprobe');% plot interpolated x and circle the values of x_up
hold on
errorbar(x_dop,x_intpl,dx_intpl,':.','Markersize',18,'LineWidth',1,'DisplayName','Interpolated')
xlim([0 0.35])
xlabel('Nominal Y content')
ylabel('Actual Y content')
title('Y content measured in $\mu$probe interpolated')
legend('show','Location','northwest')
strValues = strtrim(cellstr(num2str([x_dop(:)  x_intpl(:)],'(%.2f,%.2f)')));
text(x_dop,x_intpl,strValues,'VerticalAlignment','top','HorizontalAlignment','left','fontsize',16);

%% Transition temperatures
% The transition temperature for each substitution level is determined from 
% the maximum of dCp/dT for the corresponding sample
Tc = [2.20, 1.85, 1.56, 1.11,0,0.70,0,0,0]';%from Heat capacity data
dTc = [0.02, 0.03, 0.03, 0.02,0.4, 0.02, 0.4,0.4,0.4]';

%% Concatenate data into table
X = horzcat(x_full,dxm,dxp,Tc,dTc);
pdtable = table(x_full,dxm,dxp,Tc,dTc);

pdtable(5:6,:) = [];% remove x=0.18 and 0.22 for which data points are uncertain
%% Plot
figure()
hold on
% Y-TmVO4 transition temperatures
p0 = errorbar(pdtable.x_full,pdtable.Tc./Tc(1),pdtable.dTc./Tc(1),...
    pdtable.dTc./Tc(1),pdtable.dxm,pdtable.dxp,'.','Markersize',24);
xfit = [1, 1, 1, 0, 0];% weigh doping values in order to fit only up to x=0.1
x0 = 1;% Value of doping at which a linear fit of Tc would reach 0
% value determined after using curve fitting tool
xl = 0:0.01:x0;% range of linear plot
yl = 1/x0*(x0-xl);
p1 = plot(xl,yl,'--');

annttl = annotation('textbox',[0.15 0.175 0.2 0.1],'interpreter','latex',...
    'String',{'Tm$_{1-x}$Y$_x$VO$_4$'},'LineStyle','-','EdgeColor','k',...
    'FitBoxToText','on','LineWidth',1,'BackgroundColor','w','Color','k');% add annotation

% Cumulative plot parameters
xlim([0 0.35]); ylim([0 1.1])
ylabel('$T_D/T_D(x=0)$')
xlabel('Y content $x$')
lgd = legend([p0,p1],{'Experimental $T_D(x)$','$T_D(x=0)\cdot (1-x)$'});

%% Export figure 
cd 'C:\Users\Pierre\Desktop\Postdoc\YTmVO4\YTmVO4_phase_diagram';
printPNG('2019-05-30_YTmVO4_phase_diagram');

%% Export table
expname = '2018-11-12_YTmVO4_phase_diagram.txt';
% writetable(pdtable,expname);