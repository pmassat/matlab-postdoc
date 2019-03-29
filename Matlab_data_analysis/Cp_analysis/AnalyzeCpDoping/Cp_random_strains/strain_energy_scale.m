function d0 = strain_energy_scale(tc)
% tc is a reduced critical temperature, i.e. Tc(x)/Tc(x=0)

% % % Note: copy figure plotting from YTmVO4_FitCp_Gehring76.m % % %

pb = @(t,d) t-1/sqrt(pi)*integral(@(u)exp(-u^2)/(cosh(d*u/t)^2),-inf,inf,'ArrayValued',true);
% pb stands for phase boundary; 
% here t is a reduced temperature t = T_D/(x.lambda), where lambda = T_D(x=1)
% and d is a reduced spread in the strain distribution: d = Delta_0/(x.lambda)
% x stands for the concentration of the JT-active (Tm) ion, which therefore equals
% 1-dpg(i), since dpg(i) is the concentration of the JT-inactive (Y) ion.

delta = @(d)pb(tc,d);
% The width \Delta_0 of the gaussian strain distribution is the zero of this function
%     delta0(i) = fzero(@(x)delta(x),[0 20])
d0 = fzero(delta,[0 1.25]);% for all 0<Tc(i)/Tc(1)<1, 0<delta<1.13
% Note: only works for x<1; for x=1, d0=0 and therefore fzero will output an error
end