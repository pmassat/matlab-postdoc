%% Analyze heat capacity from DR
% This routine is intended at fitting Cp data on Y-substituted TmVO4 acquired 
% with our DR used in the Dynacool PPMS of the Lee lab
%% Average data
% Import data and compute average of data points taken at same setpoint temperature

YTmVO4_primary_analysis;% add measurement error bars!
% run script that imports and computes the average of Cp data of all
% desired compositions
%% Plot averaged Cp data for YTmVO4 x<xc
Lord = L-3;% index of the last dataset with Cp jump
plotAvgCp(avgData,dpg,1,Lord)
title(ttlCpY)
%% Plot averaged data for x>xc
%%
plotAvgCp(avgData,dpg,Lord-1,L)
title(ttlCpY)
%% Fit data
%% Fit parameters
%%
Tmfmax = [2.1,1.4,0.9,0.2,0.68];% Maximum temperature of mean-field jump for x=0 to x=0.219
% This is a visual estimate of the temperature below which Cp shows a mean-field jump. It only applies for x<xc.
Tc = [2.2,1.6,1.11,0.24,0.69];% Transition temperature; try higher value of Tc(2) for fit like Gehring 1976
Tschmin = [2.3,1.7,1.2,0.25,0.7];% temperature above which the data are essentially Schottky-like
x0 = [1e-7 0.5];% range of the pseudospin
Hc = 0.51;% critical field in Tesla, in the absence of demag factor
Tmax = [2.15; 1.51; 1.05];% temperature of the maximum of the Cp jump
%% Plot Cp data vs reduced temperature
figure; %
ax = gca;
co = ax.ColorOrder;
for i=1:3
    avgData(i).t = (avgData(i).T-Tmax(i))/Tmax(i);
    avgData(i).tp = avgData(i).t(avgData(i).t>0);
    avgData(i).tm = avgData(i).t(avgData(i).t<0);
    semilogx(avgData(i).tp,avgData(i).Cp(avgData(i).t>0),'.',...
        'Color',co(i,:),'DisplayName',sprintf('x=%.2f t$>$0',dpg(i)));
    hold on;
    semilogx(-avgData(i).tm,avgData(i).Cp(avgData(i).t<0),'x',...
        'MarkerSize',8,'Color',co(i,:),'DisplayName',sprintf('x=%.2f t$<$0',dpg(i)));
end
title([ttlCpY]);
xlabel('$|t|=\left|\frac{T-T_D}{x\cdot T_D(x=1)}\right|$');
ylabel(ylblCp);
legend('show','Location','best')
% xlim([0 1])
%% Fit mean-field jump for pure TmVO4
for i=1
    [avgData(i).fitmf, avgData(i).ffgof] = fitCpTFIM(avgData(i).T',...
        avgData(i).Cp',1./avgData(i).CpFullErr,Tc(i),0,0.95);
end

%% Fit data for x>0
%% Reproduce figure 3a of Gehring1976a 
% Equation defining the transition temperature as a function of parameter delta
td = @(d) random_strains_phase_boundary(d,0);
%% Plot phase boundary
figure;
fplot(@(d)td(d),[1e-3 1.25],'LineWidth',2)
title("T$_D$ vs spread in strain distribution");
xlabel('$\Delta_0$/($x$(Tm)$\cdot T_D$($\Delta_0$=0))'); 
ylabel('$T_D$/($x$(Tm)$\cdot T_D$($\Delta_0$=0))');

%% Compute average gap delta0 from transition temperature
d0 = ones(size(Tc));
delta0 = repmat(d0,1);
for i=2
    x = 1-dpg(i);
    tc = Tc(i)/(x*Tc(1));
    d0(i) = random_strains_energy_scale(tc,0);% Note: only works for i>1 (not for i=1)
    % do not include x in d0(i) as we want the value of d0 that corresponds to the actual value of Tc(i)
    delta0(i) = d0(i)*x*Tc(1);% delta0(i) should be of the same OM as Tc(i)
end

%% Compute the pseudospin as a function of temperature
% See Gehring1976a equation 4

%% Compute numerical arrays of sz and dsz for faster plotting with plot than fplot
% Initialize
ta1 = linspace(3e-3,tc-1e-5,1000);
sz1 = repmat(ta1,1);

%% Compute sz
for k=1:length(sz1)
    sz1(k) = order_parameter_random_strains(d0(i),ta1(k));
end
%% Plot sz
figure
plot(ta1,sz1);
title(sprintf('Order parameter vs temperature at $x$=%.2f',1-dpg(i)));
xlabel('$t=\frac{T_D}{T_D(x=1)}$');
ylabel('$\frac{\left<S^{z}\right>}{\left<S^{z}\right>_{x=1,T=0}}$');
ylim([0 1]);

%% Arrays with size that allow for computation of dsz 
ta = ta1(ta1<=tc-1e-3);
sz = sz1((ta1<=tc-1e-3));
dsz = repmat(ta,1);

%% Compute dsz 
for k=1:length(dsz)
    try
    dsz(k) = order_parameter_derivative_random_strains(d0(i),ta(k),sz(k));
    catch
        disp(k)
    end
end
%% Plot dsz
figure
plot(ta,dsz);
title(sprintf('Derivative of the order parameter vs temperature at $x$=%.2f',1-dpg(i)));
xlabel('$t=\frac{T_D}{T_D(x=1)}$');
ylabel('$\frac{d\left<S^{z}\right>}{dt}$');

%% Compute molar heat capacity Cpm
Cpma =  Cpm_random_strains(d0(i),ta,sz,dsz,x);
%% Plot Cpm
figure; hold on
errorbar(avgData(i).T/(x*Tc(1)),avgData(i).Cp/R,avgData(i).CpFullErr/R,'.','MarkerSize',18,'DisplayName',['x = ',num2str(dpg(i))])
plot(ta,Cpma);
R = 8.314;
% plot(ta,Cpintegrand(0));
% plot(ta,Cpma*1.45);
title([ttlCpY sprintf(' at $x$=%.2f',1-dpg(i))]);
xlabel('$t=\frac{T_D}{x\cdot T_D(x=1)}$');
ylabel(ylblCp);

%% arrays for fit
X = avgData(i).T/Tc(1);
X = X(X<Tc(i)/Tc(1));
Y = avgData(i).Cp(X<Tc(i)/Tc(1))/R*x;
szx = zeros(size(X));
dszx = repmat(szx,1);
%%
for kx=1:length(X)
    szx(kx) = order_parameter_random_strains(d0(i),X(kx));
end
%%
for kx=1:length(X)
    dszx(kx) = order_parameter_derivative_random_strains(d0(i),X(kx),szx(kx));
end

% 
% 
% 
% 
% 
%