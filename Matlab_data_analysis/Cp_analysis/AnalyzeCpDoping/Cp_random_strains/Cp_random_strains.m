function y = Cp_random_strains(t)
% t is the reduced temperature T/Tc(x=1)
y = zeros(size(t));
% In order to finish this function:
% - compute dsz as an independent function <=> independent script
% - comment out plots 
% Once done, try using it to fit data of YTmVO4 with curve fitting tool
%% Compute derivative of order parameter wrt temperature
ds = @(s1,t) s1-1/sqrt(pi)*integral(@(u)(s1./t-(sz(t)+u.*delta0(i)/Tc(1))./(t.^2))...
    .*exp(-u.^2)./cosh((sz(t)+u.*delta0(i)/Tc(1))./t).^2,-Inf,Inf);
dsz = @(t) fzero(@(s1)ds(s1,t),[-30 0]);

%% Plot ds
figure; hold on;
tplot = 0.8;
fplot(@(s)ds(s,tplot))
line(xlim,[0 0],'color','black','linestyle','--')

%% Plot dsz
figure
tplot = 0.8;
fplot(@(t)dsz(t),[1e-2 maxT-1e-3])
title(sprintf('Derivative of the order parameter vs temperature at $x$=%.2f',1-dpg(i)));
xlabel('$t=\frac{T_D}{T_D(x=1)}$');
ylabel('$\frac{\left<S^{z}\right>}{\left<S^{z}\right>_{x=1,T=0}}$');

%% Compute numerical arrays of sz and dsz for faster plotting with plot than fplot
% Initialize
ta = linspace(1e-2,maxT-1e-3,1000);
sza = repmat(ta,1);
dsza = repmat(ta,1);
%% Compute for sz
for k=1:length(sza)
    sza(k) = sz(ta(k));
end
%% Same with dsz 
for k=1:length(sza)
    dsza(k) = dsz(ta(k));
end

%% Compute molar heat capacity Cpm
x = 1-dpg(i);
d0r = delta0(i)/Tc(1);
Er = @(t,u) (x.*sz(t)+u.*d0r)./t;
Era = @(u) (x.*sza+u.*d0r)./ta;
%% Compute integral of molar heat capacity
%% Compute integrand
sintegrand = @(u)-(x.*sza/2 + u*d0r).*tanh(Era(u));
%% Molar entropy
sm = x/sqrt(pi).*integral(@(u)exp(-u^2).*sintegrand(u),...
    -Inf,Inf,'ArrayValued',true);
%% Functional form of Cp
integrand = @(t,u) -(x/2.*dsz(t).*tanh(Er(t,u))+...
    (x.*sz(t)/2 + u*d0r).*(x/t.*dsz(t)-Er(t,u)./t)./cosh(Er(t,u)).^2);
Cpm = @(t) x/sqrt(pi)*integral(@(u)exp(-u.^2).*integrand(t,u),-Inf,Inf);
%% Array computation (faster than functional form)
Cpintegrand = @(u) -(x/2.*dsza.*tanh(Era(u))+...
    (x.*sza/2 + u*d0r).*(x.*dsza-Era(u))./ta./cosh(Era(u)).^2);
Cpma = x/sqrt(pi)*integral(@(u)exp(-u.^2).*Cpintegrand(u),-Inf,Inf,'ArrayValued',true);
%% Plot Cpm
figure; hold on
% fplot(@(t)Cpm(t),[1e-2 maxT-1e-3]);
plot(ta,Cpma);
R = 8.314;
errorbar(avgData(i).T/Tc(1),avgData(i).Cp/R,avgData(i).CpFullErr/R,'.','MarkerSize',18,'DisplayName',['x = ',num2str(dpg(i))])
plot(ta,Cpma*1.675);

end



% 
% 
% 
% 
% 
%

