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
        avgData(i).Cp',1./avgData(i).stdCp,Tc(i),0,0.95);
end

%% Fit data for x>0
%% Compute average gap delta0 from transition temperature
% See Gehring1976a equation 5.

pb = @(t,d) 1/sqrt(pi)*integral(@(u)exp(-u^2)/(cosh(d*u/t)^2),-inf,inf,'ArrayValued',true)-t;
% pb stands for phase boundary; 
% here t is a reduced temperature t = T/(x.lambda), where lambda is defined in the paper
% and d is a reduced spread in the strain distribution: d = Delta_0/(x.lambda)
% It is not clear to me however what x stands for. If x were the doping
% level, this expression would not make sense for pure TmVO4.

%% 
% x*lambda = T_D(Delta_0=0) = T_D(x=0) = T_D(pureTmVO4) = Tc(1) with 
% the notations of this code
figure
fplot(@(t)pb(t,0.1) ,[1e-3 5])
title('Function defining T$_Q$ for $\Delta_0$/(k$_B\cdot$dpg(i)$\cdot\lambda$)=0.1');
xlabel('$\tau$ = T/(dpg(i)$\cdot\lambda$)'); ylabel('pb(T,d=0.1)');

%% Reproduce figure 3a of Gehring1976a
td = @(d) fzero(@(t)pb(t,d),[1e-3 5]);
figure;
fplot(@(d)td(d),[1e-3 1.25])
title("Transition temperature as a function of spread in strain distribution");
xlabel('$\Delta_0$/(x$\cdot\lambda$)'); ylabel('$\tau_D$/$\tau_D$($\Delta_0$=0)');
ylim([0 1]);

%%
delta0 = ones(size(Tc));
for i=2
%     delta = @(x) (1-dpg(i))*lambda/sqrt(pi)*integral(@(u)exp(-u^2)/(cosh(x*u/Tc(i))^2),-inf,inf,'ArrayValued',true)-Tc(i)
    delta = @(d)pb(Tc(i)/Tc(1),d);
% The width \Delta_0 of the gaussian strain distribution is the zero of this function
%     delta0(i) = fzero(@(x)delta(x),[0 20])
    d0 = fzero(delta,[0 20]);% note: only works for i>1 (not for i=1)
    fplot(delta,[0 10])
    line(xlim,[0 0],'color','black','linestyle','--')
    line([d0 d0],ylim,'color','red','linestyle','--')
    title("Function defining \delta_0 = \Delta_0/(k_B\cdotdpg(i)\cdot\lambda)");
    xlabel("x");ylabel("f(x)-k_B\cdotT_D")
    delta0(i) = d0*Tc(i);
end
%% Compute the pseudospin as a function of temperature
% See Gehring1976a equation 4

syms x T u% symbolic math used to derive the expression of pseudospin sigma 
for i=2
%     sigma = @(x,T) 1/sqrt(pi)*integral(@(u)exp(-u^2)*tanh((delta0(i)*u+lambda*x)./T),-inf,inf,'ArrayValued',true)-x
%     sigma = @(x,T) 1/sqrt(pi)*vpaintegral(exp(-u^2)*tanh((delta0(i)*u+lambda*x)./T),u,[-inf,inf])-x
    sigma = @(T,x) 1/sqrt(pi)*int(exp(-u^2)*tanh((delta0(i)*u+x*Tc(1)/(1-dpg(i)))./T),u,[-inf,inf])-x
%% 
% Sigma is the self-consistent equation defining pseudospin x as a function 
% of temperature T. The value of the pseudospin Sz is simply the solution of sigma=0, 
% hence the use of fzero below.

    sgm = matlabFunction(sigma(T,x))% converting into non-symbolic function for computation of values
    sz = @(T) fzero(@(x)sgm(T,x),[1e-3 1])
    Tplot = 1;
    fplot(@(x)sgm(Tplot,x),[0 1])
    line(xlim,[0 0],'color','black','linestyle','--')
    line([sz(Tplot) sz(Tplot)],ylim,'color','red','linestyle','--')
    title("Function defining pseudospin Sz at T=1K");
    xlabel("x");ylabel("f(x)-x");
end
%% 
% Test function sz

sz(1e-2)
sz(1.6)
sz(2.1)

%% Compute the derivative of sigma, to determine the max temperature of the ordering
% sigma goes to zero when increasing T from 0, until dsigma/dx(T,x=0) becomes 
% negative, then the equation sigma(T,x>0)=0 does not have any solution anymore. 
% Therefore, the max temperature of the ordering is reached when dsigma/dx(T,x=0)=0.

dsigma = diff(sigma(T,x),x)
g = matlabFunction(dsigma)
%%
fplot(@(x)g(1,x),[0 1])
line(xlim,[0 0],'color','black','linestyle','--')
% line([sz(Tplot) sz(Tplot)],ylim,'color','red','linestyle','--')
title(["dsigma/dx(T,x) at T=",num2str(Tplot),"K"]);
xlabel("x");ylabel("d\sigma/dx(T=1K,x)");
g(1,0.39)
%%
fplot(@(T)g(T,0),[0 4])
g(1e-3,0)
maxT = fzero(@(T)g(T,0),[1e-3 4])
%%
figure
fplot(sz,[1e-2 2])
%% 
% Need to check: delta0(i) should be of the same OM as Tc(i)
% 
% 
% 
% 
% 
%