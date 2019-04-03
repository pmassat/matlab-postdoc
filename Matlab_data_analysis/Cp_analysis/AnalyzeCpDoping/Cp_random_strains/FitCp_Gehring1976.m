%% Import data
% from 'C:\Users\Pierre\Desktop\Postdoc\Bibliographic_resources\Articles_CJTE\RVO4\TmVO4\TmVO4_heat_capacity\Gehring1976a_parts\Gehring1976a_fig2a.txt';
data = Gehring1976afig2a;
x = 0.9;
Tc0 = 2.15;
tc = 1.58/Tc0;
d0 = random_strains_energy_scale(tc);
%% Compute numerical arrays of sz and dsz for faster plotting with plot than fplot
% Initialize
ta1 = linspace(1.5e-3,tc-1e-5,1000);
sz1 = repmat(ta1,1);

%% Compute sz
for k=1:length(sz1)
    sz1(k) = order_parameter_random_strains(d0,ta1(k));
end
%% Plot sz
figure
plot(ta1,sz1);
title(sprintf('Order parameter vs temperature at $x$=%.2f',1-dpg));
xlabel('$t=\frac{T_D}{T_D(x=1)}$');
ylabel('$\frac{\left<S^{z}\right>}{\left<S^{z}\right>_{x=1,T=0}}$');
ylim([0 1]);

%% Arrays with size that allow for computation of dsz 
ta = ta1(ta1<=tc-1e-3);
sz = sz1((ta1<=tc-1e-3));
dsz = repmat(ta,1);

%% Compute dsz 
for k=1:length(dsz)
    dsz(k) = order_parameter_derivative_random_strains(d0,ta(k),sz(k));
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
figure; hold on
% plot(ta,Cpintegrand(0));
plot(data.T/2.15,data.Cp,'.')
plot(ta,Cpma);
% plot(ta,Cpma*1.45);

%% Compute Cp from magnetic dipole interactions
% Gehring et al. 1976, equation 7
syms t
M = 1/(t*cosh(t)^2)+tanh(1/t);
d1M = diff(M,t);
d2M = diff(M,t,2);
fd2M = matlabFunction(d2M);
g = @(t,u) 1/(d0*u)*fd2M(t/(u*d0));
C_0 = 0.088;% 1/T^2 coefficient, see paper
%%
Cpmd = @(t) t*x/(4*sqrt(pi))*C_0/Tc0^2*integral(@(u)exp(-u.^2).*g(t,u),-Inf,Inf,'ArrayValued',true);

%%
figure
fplot(Cpmd,[0 2])














