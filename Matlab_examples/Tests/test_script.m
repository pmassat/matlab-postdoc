% % generate data
% x = 0:.1:10;
% y = x.*x + randn(size(x));
% w = linspace(0.1, 10.1,length(x));
% x = x(:);
% y = y(:);
% w = w(:);
% %plot data
% figure
% errorbar(x,y,1./w,'.');
% %fit
% ft = fittype('poly2');
% cf = fit(x,y,ft,'Weight',w);
% % Plot fit
% hold on
% plot(cf,'fit',0.95);
%%
delta0 = 0.5;
ea = 0.01;
x = 0.9;
t = linspace(1e-3,1,500);
sza = zeros(size(t));
% function [sz,maxT] = order_parameter_random_strains(delta0,ea,temp)
% delta0 is a *reduced* energy scale of random strains, delta0 = Delta0/(x*Tc(1))
% temp = T/Tc(1) is a reduced temperature
% sgm = @(s,temp)s-x/sqrt(pi).*integral(@(u)tanh((s+u.*delta0)./temp).*...
%     exp(-(u+ea).^2),-Inf,Inf,'ArrayValued',true);
%%
tc = random_strains_phase_boundary(delta0,ea);

%% Compute the pseudospin as a function of temperature
% See Gehring1976a equation 4
 
% % s is the variable that defines the order parameter sz that goes from 0 to 1
% % t = T/Tc(1) is a reduced temperature, which also goes from 0 to 1
% %     sgm = matlabFunction(sigma);% converting into non-symbolic function
% %     for computation of values (not necessary when using direct definition below)
for i=1:length(t)
    temp = t(i);
% Sigma and sgm are self-consistent equations defining pseudospin as a function 
% of reduced temperature t. The value of the pseudospin Sz is simply the 
% solution of sigma=0, hence the use of fzero below.
%     sgm = @(s)s-1/sqrt(pi).*integral(@(u)tanh((s+u.*delta0)./temp).*...
%         exp(-(u+ea).^2),-Inf,Inf,'ArrayValued',true);
    try
%     sz = fzero(@(s)sgm(s,temp),1);
    sza(i) = order_parameter_random_strains_offset(delta0,ea,t(i));
    catch
        disp(i); break
    end
end

%% Important note
% the alphabetical order of function variables matters to Matlab! Always
% check that the order of variables when calling a function matches the 
% function the definition of the function in the workspace.

%% Plot the order parameter vs reduced temperature
% fplot is slow; it is faster to compute an array of sz and use plot
figure; 
plot(t,sza);
title(sprintf('Order parameter vs temperature'));
xlabel('$t=\frac{T_D}{T_D(x=1)}$');
ylabel('$\frac{\left<S^{z}\right>}{\left<S^{z}\right>_{x=1,T=0}}$');
ylim([0 1]);

%% 
tcomp = t(250);
[sz,sgm] = order_parameter_random_strains_offset(delta0,ea,tcomp);
%% plot sgm
    figure
    tplot = tcomp;
    fplot(@(s)sgm(s),[0 1])
    line(xlim,[0 0],'color','black','linestyle','--')
%     line([sz(tplot) sz(tplot)],ylim,'color','red','linestyle','--')
    title('$f_{\sigma}=0$ when $x = \sigma_z$');
    xlabel('$x$');
    ylabel(['$f_{\sigma}(x) = x-\frac{1}{\sqrt\pi} \int e^{-u^2}\cdot \tanh \left(\frac{1}{t}\cdot'...
        '\left(x+\frac{\delta_0(i)\cdot u}{Tc(1)}\right) \right)$d$u$']);
    legend(sprintf('t=%.1f',tplot));
