%% Import data
% from 'C:\Users\Pierre\Desktop\Postdoc\Bibliographic_resources\Articles_CJTE\RVO4\TmVO4\TmVO4_heat_capacity\Gehring1976a_parts\Gehring1976a_fig2a.txt';
data = Gehring1976afig2a;
x = 0.9;
Tc0 = 2.15;
% tc = 1.58/Tc0;
% d0 = random_strains_energy_scale(tc);
da = 0.0;
d0paper = 1.14/(x*Tc0);
tcpaper = random_strains_phase_boundary(d0paper,da);

tc_onset = 1.65/(x*Tc0);
d0_onset = random_strains_energy_scale(tc_onset,da);
ta_onset = linspace(3e-3,tc_onset-1e-3,1000);

d0_test = 0.01;
tc_test = random_strains_phase_boundary(d0_test,da);
ta_test = linspace(3e-3,tc_test-1e-3,1000);

%% Compute numerical arrays of sz and dsz for faster plotting with plot than fplot
% Initialize
ta1 = linspace(3e-3,tcpaper-1e-3,1000);
ta_paper = ta1(ta1<=tcpaper-1e-3);
sz1 = repmat(ta1,1);

%% Compute sz
for k=1:length(sz1)
    try
    sz1(k) = order_parameter_random_strains(d0paper,ta1(k));
    catch
        disp(k)
    end
end
%% Plot sz
figure
plot(ta1*x,sz1*x);
title(sprintf('Order parameter vs temperature at $x$=%.2f',x));
xlabel('$t=\frac{T_D}{T_D(x=1)}$');
ylabel('$\frac{\left<S^{z}\right>}{\left<S^{z}\right>_{x=1,T=0}}$');
xlim([0 1]); ylim([0 1]);

%% Arrays with size that allow for computation of dsz 
% sz = sz1((ta1<=tc-1e-3));
% dsz = repmat(ta,1);
% 
% %% Compute dsz 
% for k=1:length(dsz)
%     try
%     dsz(k) = order_parameter_derivative_random_strains(d0,ta(k),sz(k));
%     catch
%         disp(k)
%     end
% end
% %% Plot dsz
% figure
% plot(ta,dsz);
% title(sprintf('Derivative of the order parameter vs temperature at $x$=%.2f',x));
% xlabel('$t=\frac{T_D}{T_D(x=1)}$');
% ylabel('$\frac{d\left<S^{z}\right>}{dt}$');

%% Compute molar heat capacity Cpm
% Cpma =  Cpm_random_strains(d0,ta,sz,dsz,x);
Cpma = Cp_sz_dsz_random_strains(d0paper,ta_paper,x);

%%
Cpma2 = Cp_sz_dsz_random_strains(d0paper,ta_paper,x);
%%
Cpm_onset = Cp_sz_dsz_random_strains(d0_onset,ta_onset,x);

%%
Cpm_test = Cp_sz_dsz_random_strains(d0_test,ta_test,x);

%% Plot Cpm
Tplot = ta_paper*x*Tc0;
figure; hold on
% plot(ta,Cpintegrand(0));
plot(ta_paper*x*Tc0,Cpma*1.4,'DisplayName','Paper fit');
plot(data.T,data.Cp,'.','DisplayName','Data')
plot(ta_paper*x*Tc0,Cpma,'DisplayName','Eqn (6)');
title(sprintf('Heat capacity of Tm$_{%.1f}$Lu$_{%.1f}$VO$_4$',x,1-x));
xlabel('$T$ (K)');
ylabel('C$_p$ (J$\cdot$mol$^{-1}\cdot$K$^{-1}$)');
legend('show','Location','northwest');
xlim([0 3]); ylim([0 1.4]);
plot(ta_test*x*Tc0,Cpm_test*1.2,'DisplayName','Test');
plot(ta_onset*x*Tc0,Cpm_onset*1.4,'DisplayName','T$_{c,onset}$');
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
g2 = @(t,u) - 6*d0paper*(u*d0paper)^4/(t^5*cosh(u*d0paper/t)^4) -...
    10*d0paper*(u*d0paper)^3*tanh(u*d0paper/t)/(t^4*cosh(u*d0paper/t)^2) +...
    4*d0paper*(u*d0paper)^4/(t^5*cosh(u*d0paper/t)^2) +...
    4*d0paper*(u*d0paper)^2/(t^3*cosh(u*d0paper/t)^2);
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












