function y = random_strains_phase_boundary_equation(delta,ea,t)
% See Gehring1976a equation 5.
% t is a reduced temperature t = T_D/(x.lambda), where lambda = T_D(x=1)
% and d is a reduced spread in the strain distribution: d = Delta_0/(x.lambda)
% x stands for the concentration of the JT-active (Tm) ion, which therefore equals
% 1-dpg(i), since dpg(i) is the concentration of the JT-inactive (Y) ion.
y = t-1/sqrt(pi)*integral(@(u)exp(-(u+ea)^2)/(cosh(delta*u/t)^2),-inf,inf,'ArrayValued',true);

% % % Plot equation vs reduced temperature for a given spread in the strain
% % % distribution, just to check that it performs well
% % % See Gehring1976a: lambda = T_D(Delta_0=0) = T_D(x=1) = T_D(pureTmVO4) 
% pb = @(t,d) random_strains_phase_boundary_equation(d,t);
% figure
% fplot(@(t)pb(t,0.1) ,[-1 1])
% title('$y = 0$ when $\tau = \tau_D$');% 
% xlabel('$\tau$ = T/($x$(Tm)$\cdot\lambda$) = T/($x$(Tm)$\cdot$T$_D$(x=1))'); 
% ylabel('$y(\tau) = \tau - \frac{1}{\sqrt \pi} \int e^{-u^2}/\cosh^2(\frac{\Delta_0}{x\cdot\lambda}\cdot \frac{u}{\tau}) \cdot$d$u$');
% legend('$\Delta_0$/($x\cdot\lambda$)=0.1');
end