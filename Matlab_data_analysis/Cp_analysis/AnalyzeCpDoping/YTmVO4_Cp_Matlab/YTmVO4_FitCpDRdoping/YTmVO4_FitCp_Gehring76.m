%% Analyze heat capacity from DR
% This routine is intended at fitting Cp data on Y-substituted TmVO4 acquired 
% with our DR used in the Dynacool PPMS of the Lee lab
%% Average data
% Import data and compute average of data points taken at same setpoint temperature

YTmVO4_primary_analysis;% add measurement error bars!
% run script that imports and computes the average of Cp data of all
% desired compositions
%% Plot averaged data
%% Plot averaged Cp data for YTmVO4 x<xc
%%
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
Tc = [2.2,1.57,1.11,0.24,0.69];% Transition temperature
Tschmin = [2.3,1.7,1.2,0.25,0.7];% temperature above which the data are essentially Schottky-like
x0 = [1e-7 0.5];% range of the pseudospin
Hc = 0.51;% critical field in Tesla, in the absence of demag factor

%% Fit mean-field jump for pure TmVO4
for i=1
    [avgData(i).fitmf, avgData(i).ffgof] = fitCpTFIM(avgData(i).T',...
        avgData(i).Cp',1./avgData(i).CpFullErr,Tc(i),0,0.95);
end

%% Fit data for x>0
%% Compute average gap delta0 from transition temperature
% See Gehring1976a equation 5.

pb = @(t,d) t-1/sqrt(pi)*integral(@(u)exp(-u^2)/(cosh(d*u/t)^2),-inf,inf,'ArrayValued',true);
% pb stands for phase boundary; 
% here t is a reduced temperature t = T_D/(x.lambda), where lambda = T_D(x=1)
% and d is a reduced spread in the strain distribution: d = Delta_0/(x.lambda)
% x stands for the concentration of the JT-active (Tm) ion, which therefore equals
% 1-dpg(i), since dpg(i) is the concentration of the JT-inactive (Y) ion.

%% 
% x*lambda = T_D(Delta_0=0) = T_D(x=0) = T_D(pureTmVO4) = Tc(1) with 
% the notations of this code
figure
fplot(@(t)pb(t,0.1) ,[-1 1])
title('$y = 0$ when $\tau = \tau_D$');% 
xlabel('$\tau$ = T/($x$(Tm)$\cdot\lambda$) = T/($x$(Tm)$\cdot$T$_D$(x=1))'); 
ylabel('$y(\tau) = \tau - \frac{1}{\sqrt \pi} \int e^{-u^2}/\cosh^2(\frac{\Delta_0}{x\cdot\lambda}\cdot \frac{u}{\tau}) \cdot$d$u$');
legend('$\Delta_0$/($x\cdot\lambda$)=0.1');

%% Reproduce figure 3a of Gehring1976a 
% Equation defining the transition temperature as a function of parameter delta
td = @(d) fzero(@(t)pb(t,d),[-2 2]);
%% Plot
figure;
fplot(@(d)td(d),[1e-3 1.25])
title("T$_D$ vs spread in strain distribution");
xlabel('$\Delta_0$/($x$(Tm)$\cdot T_D$($\Delta_0$=0))'); 
ylabel('$T_D$/($x$(Tm)$\cdot T_D$($\Delta_0$=0))');
% ylim([0 1]);

%%
delta0 = ones(size(Tc));
figure; hold on;
for i=2
%     delta = @(x) (1-dpg(i))*lambda/sqrt(pi)*integral(@(u)exp(-u^2)/(cosh(x*u/Tc(i))^2),-inf,inf,'ArrayValued',true)-Tc(i)
    delta = @(d)pb(Tc(i)/Tc(1),d);
% The width \Delta_0 of the gaussian strain distribution is the zero of this function
%     delta0(i) = fzero(@(x)delta(x),[0 20])
    d0 = fzero(delta,[0 1.25]);% for all 0<Tc(i)/Tc(1)<1, 0<delta<1.13
    % Note: only works for i>1 (not for i=1)
    fplot(delta,[0 1.25])
    line(xlim,[0 0],'color','black','linestyle','--')
    line([d0 d0],ylim,'color','red','linestyle','--')
    title('$y=0$ when $\delta$ equals the strain distribution $\delta_0$');
    xlabel('$\delta = \Delta/(x$(Tm)$\cdot\lambda)$');
    ylabel('$y(\delta) = \tau - \frac{1}{\sqrt \pi} \int e^{-u^2}/\cosh^2(\delta \cdot \frac{u}{\tau}) \cdot$d$u$');
    delta0(i) = d0*Tc(i);% delta0(i) should be of the same OM as Tc(i)
end

%% Compute the pseudospin as a function of temperature
% See Gehring1976a equation 4

syms s t% symbolic math used to derive the expression of pseudospin sz 
%%
i=2;
    %% compute
%     sigma = @(x,T) 1/sqrt(pi)*integral(@(u)exp(-u^2)*tanh((delta0(i)*u+lambda*x)./T),-inf,inf,'ArrayValued',true)-x
%     sigma = @(x,T) 1/sqrt(pi)*vpaintegral(exp(-u^2)*tanh((delta0(i)*u+lambda*x)./T),u,[-inf,inf])-x
    sigma = @(t,s) s-1/sqrt(pi)*int(exp(-u^2)*tanh(1/t*(s+delta0(i)*u/Tc(1))),u,[-inf,inf]);
% s is the variable that defines the order parameter sz that goes from 0 to 1
% t = T/Tc(1) is a reduced temperature, which also goes from 0 to 1

% Sigma is the self-consistent equation defining pseudospin x as a function 
% of temperature T. The value of the pseudospin Sz is simply the solution of sigma=0, 
% hence the use of fzero below.
    sgm = matlabFunction(sigma(t,x));% converting into non-symbolic function for computation of values
    sz = @(t) fzero(@(x)sgm(t,x),[1e-7 1]);
    
    %% plot
    figure
    tplot = 0.8;
    fplot(@(x)sgm(tplot,x),[0 1])
    line(xlim,[0 0],'color','black','linestyle','--')
    line([sz(tplot) sz(tplot)],ylim,'color','red','linestyle','--')
    title('$y=0$ when $x = \sigma_z$');
    xlabel('$x$');
    ylabel(['$f_{\sigma}(x) = x-\frac{1}{\sqrt\pi} \int e^{-u^2}\cdot \tanh \left(\frac{1}{t}\cdot'...
        '\left(x+\frac{\delta_0(i)\cdot u}{Tc(1)}\right) \right)$d$u$']);
    legend(sprintf('T/T$_c$(x=1)=%.1f',tplot));

%% 
% Test function sz

sz(1e-2)
sz(0.5)
sz(0.99)

%% Compute the derivative of sigma, to determine the max temperature of the ordering
% sigma goes to zero when increasing T from 0, until dsigma/dx(T,x=0) becomes 
% positive, then the equation sigma(T,x>0)=0 does not have any solution anymore. 
% Therefore, the max temperature of the ordering is reached when dsigma/dx(T,x=0)=0.

dsigma = diff(sigma(t,x),x);
g = matlabFunction(dsigma);
%% Plot derivative of f_sigma wrt to x at a given temperature
figure
fplot(@(x)g(tplot,x),[0 1])
line(xlim,[0 0],'color','black','linestyle','--')
% line([sz(Tplot) sz(Tplot)],ylim,'color','red','linestyle','--')
title('Derivative of $f_{\sigma}$');
xlabel('x');
ylabel(sprintf('$\frac{\partial f_{\sigma}}{\partial x}$(t=%.1f,x)',tplot));
g(tplot,0)

%% Plot df_sigma/dx(x=0) vs temperature
figure 
fplot(@(t)g(t,0),[0 1])
g(1e-3,0)
%% Determine max value of reduced temperature for which the order parameter is defined
maxT = fzero(@(t)g(t,0),[1e-3 1]);

%% Plot the order parameter vs reduced temperature
figure;
fplot(sz,[1e-2 maxT-1e-3]);
ylim([0 1])

%% Compute derivative of order parameter wrt temperature
dsz = diff(@(t)sz(t),t);% does not work; compute dsz self-consistently











% 
% 
% 
% 
% 
%