function d0 = random_strains_energy_scale(tc,ea,varargin)
% tc is a reduced critical temperature, i.e. Tc(x)/Tc(x=1), with x being
% the concentration of the JT-active ion (Tm)

% % % Note: copy figure plotting from YTmVO4_FitCp_Gehring76.m % % %

delta = @(d) random_strains_phase_boundary_equation(d,ea,tc);
% The width \Delta_0 of the gaussian strain distribution is the zero of this function
%     delta0(i) = fzero(@(x)delta(x),[0 20])
d0 = fzero(delta,[0 1.25]);% for all 0<tc<1, 0<delta<1.13
% Note: only works for x<1; for x=1, d0=0 and therefore fzero will output an error

if nargin>2
figure; hold on;
fplot(delta,[0 1.25]);
line(xlim,[0 0],'color','black','linestyle','--');
line([d0 d0],ylim,'color','red','linestyle','--');
title('$y=0$ when $\delta$ equals the strain distribution $\delta_0$');
xlabel('$\delta = \Delta/(x$(Tm)$\cdot\lambda)$');
ylabel('$y(\delta) = \tau - \frac{1}{\sqrt \pi} \int e^{-u^2}/\cosh^2(\delta \cdot \frac{u}{\tau}) \cdot$d$u$');
end
end