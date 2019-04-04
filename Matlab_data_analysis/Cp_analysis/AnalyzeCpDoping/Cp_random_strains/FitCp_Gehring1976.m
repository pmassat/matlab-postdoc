%% Import data
% from 'C:\Users\Pierre\Desktop\Postdoc\Bibliographic_resources\Articles_CJTE\RVO4\TmVO4\TmVO4_heat_capacity\Gehring1976a_parts\Gehring1976a_fig2a.txt';
data = Gehring1976afig2a;
x = 0.9;
Tc0 = 2.15;
% tc = 1.58/Tc0;
% d0 = random_strains_energy_scale(tc);
d0 = 1.14/2.15;
tc = random_strains_phase_boundary(d0);
% Try using tc = x;
%% Compute numerical arrays of sz and dsz for faster plotting with plot than fplot
% Initialize
ta1 = linspace(3e-3,tc-1e-5,1000);
sz1 = repmat(ta1,1);

%% Compute sz
for k=1:length(sz1)
    sz1(k) = order_parameter_random_strains(d0,ta1(k));
end
%% Plot sz
figure
plot(ta1,sz1);
title(sprintf('Order parameter vs temperature at $x$=%.2f',x));
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
    dsz(k) = order_parameter_derivative_random_strains(d0,ta(k),sz(k));
    catch
        disp(k)
    end
end
%% Plot dsz
figure
plot(ta,dsz);
title(sprintf('Derivative of the order parameter vs temperature at $x$=%.2f',x));
xlabel('$t=\frac{T_D}{T_D(x=1)}$');
ylabel('$\frac{d\left<S^{z}\right>}{dt}$');

%% Compute molar heat capacity Cpm
Cpma =  Cpm_random_strains(d0,ta,sz,dsz,x);
%% Plot Cpm
Tplot = ta*Tc0;
figure; hold on
% plot(ta,Cpintegrand(0));
plot(data.T,data.Cp,'.','DisplayName','Data')
plot(Tplot,Cpma,'DisplayName','Eqn (6)');
title(sprintf('Heat capacity of Tm$_{%.1f}$Lu$_{%.1f}$VO$_4$',x,1-x));
xlabel('$T$ (K)');
ylabel('C$_p$ (J$\cdot$mol$^{-1}\cdot$K$^{-1}$)');
legend('show','Location','northwest');
plot(Tplot,Cpma*1.4);

%% Compute Cp from magnetic dipole interactions
% Gehring et al. 1976, equation 7
syms t
M = 1/(t*cosh(1/t)^2)+tanh(1/t);
fm = matlabFunction(M);
d1M = diff(M,t);
fd1m = matlabFunction(d1M);
d2M = diff(M,t,2);
fd2M = matlabFunction(d2M);
% g = @(t,u) 1/(d0*u).^2.*fd2M(t./(u*d0));% doesn't work, most likely for numerical reasons
g2 = @(t,u) - 6*d0*(u*d0)^4/(t^5*cosh(u*d0/t)^4) -...
    10*d0*(u*d0)^3*tanh(u*d0/t)/(t^4*cosh(u*d0/t)^2) +...
    4*d0*(u*d0)^4/(t^5*cosh(u*d0/t)^2) +...
    4*d0*(u*d0)^2/(t^3*cosh(u*d0/t)^2);
% g2 is obtained as follows: use the output of simplify(d2M) and rearrange
% it into 4 distinct terms, then replace all occurrences of 't' by
% t/(u*d0), then dividing the whole expression by u as the full
% integrand of equantion 7 of Gerhing et al 1976 is 1/u*M
C_0 = 0.088;% 1/T^2 coefficient, see paper
%%
Cpmd = @(t) t*x/(4*sqrt(pi))*C_0/Tc0*integral(@(u)exp(-u.^2).*g2(t,u),-Inf,Inf,'ArrayValued',true);

%%
figure
fplot(Cpmd,[1e-3 2])

%%
Mu = @(u) 1/(cosh(u)^2)+tanh(u)/u;
figure
fplot(Mu, [-20 20])












