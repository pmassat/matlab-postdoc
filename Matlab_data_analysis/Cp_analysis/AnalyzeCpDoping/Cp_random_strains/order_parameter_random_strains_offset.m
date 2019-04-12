function [sz,sgm] = order_parameter_random_strains_offset(delta0,ea,temp)
% delta0 is a *reduced* energy scale of random strains, delta0 = Delta0/(x*Tc(1))
% temp = T/Tc(1) is a reduced temperature

%% Compute the pseudospin as a function of temperature
% See Gehring1976a equation 4

% % s is the variable that defines the order parameter sz that goes from 0 to 1
% % t = T/Tc(1) is a reduced temperature, which also goes from 0 to 1
% %     sgm = matlabFunction(sigma);% converting into non-symbolic function
% %     for computation of values (not necessary when using direct definition below)

% Sigma and sgm are self-consistent equations defining pseudospin as a function 
% of reduced temperature t. The value of the pseudospin Sz is simply the 
% solution of sigma=0, hence the use of fzero below.
    sgm = @(s)s-1/sqrt(pi).*integral(@(u)tanh((s+u.*delta0)./(temp)).*exp(-(u-ea).^2),-Inf,Inf,'ArrayValued',true);
%     sz = fzero(@(s)sgm(s),[1e-7 1]);% this doesn't work for negative
%     values of the offset strain ea
    sz = fzero(@(s)sgm(s),1);
%     tc = random_strains_phase_boundary(delta0);

%% Important note about the definition of sz
% Whether sgm should be defined as s-x/sqrt(pi)... or s-1/sqrt(pi)... as
% above depends on the definition of the reduced quantities delta0 and
% temp: if these quantities are defined over Tc0 (which is the transition
% temperature of the parent compound), then one should include the x in the
% definition of sgm; however if delta0 = Delta0/(x*Tc0) and temp = temp/(x*Tc0), 
% i.e. if they are defined per JT-active ion, then sgm should not include x

%% Important note about Matlab
% the alphabetical order of function variables matters to Matlab! Always
% check that the order of variables when calling a function matches the 
% function the definition of the function in the workspace.
 
%% Plot the order parameter vs reduced temperature
% % fplot is slow; it is faster to compute an array of sz and use plot
% figure;
% fplot(sz,[1e-2 tc-1e-3]);
% title(sprintf('Order parameter vs temperature at $x$=%.2f',1-dpg(i)));
% xlabel('$t=\frac{T_D}{T_D(x=1)}$');
% ylabel('$\frac{\left<S^{z}\right>}{\left<S^{z}\right>_{x=1,T=0}}$');
% ylim([0 1]);

%% plot sgm
%     figure
%     tplot = 0.8*tc;
%     fplot(@(s)sgm(s,tplot),[0 1])
%     line(xlim,[0 0],'color','black','linestyle','--')
%     line([sz(tplot) sz(tplot)],ylim,'color','red','linestyle','--')
%     title('$f_{\sigma}=0$ when $x = \sigma_z$');
%     xlabel('$x$');
%     ylabel(['$f_{\sigma}(x) = x-\frac{1}{\sqrt\pi} \int e^{-u^2}\cdot \tanh \left(\frac{1}{t}\cdot'...
%         '\left(x+\frac{\delta_0(i)\cdot u}{Tc(1)}\right) \right)$d$u$']);
%     legend(sprintf('t=%.1f',tplot));

%% Compute the derivative of sigma, to determine the max temperature of the ordering
% % sigma goes to zero when increasing T from 0, until dsigma/dx(T,x=0) becomes 
% % positive, then the equation sigma(T,x>0)=0 does not have any solution anymore. 
% % Therefore, the max temperature of the ordering is reached when dsigma/dx(T,x=0)=0.
% syms s t u% symbolic math used to derive the expression of pseudospin sz 
% fsigma = s-1/sqrt(pi)*int(exp(-u^2)*tanh((s+delta0*u)./t),u,[-inf,inf]);
% dsigma = diff(fsigma,s);
% dfsigma = matlabFunction(dsigma);

%% Determine max value of reduced temperature for which the order parameter is defined
% % Note: because there is a univoque correspondence between d0 and Tc(i)/Tc,
% % maxT should always be equal to Tc(i)/Tc(1). Hence this is only a verification.
% maxT = fzero(@(t)dfsigma(0,t),[1e-3 1]);

%% Plot derivative of f_sigma wrt to x vs x at a given temperature
% figure
% fplot(@(s)dfsigma(s,tplot),[0 1])
% line(xlim,[0 0],'color','black','linestyle','--')
% % line([sz(Tplot) sz(Tplot)],ylim,'color','red','linestyle','--')
% title('$s$ derivative of $f_{\sigma}(s,t)$ vs $s$');
% xlabel('s');
% ylabel(['$\frac{\partial f_{\sigma}}{\partial s}$' sprintf('(t=%.1f,$s$)',tplot)]);
% dfsigma(0,tplot)

%% Plot df_sigma/dx(x=0) vs temperature
% figure 
% fplot(@(t)dfsigma(0,t),[0 1])
% line(xlim,[0 0],'color','black','linestyle','--');
% title('$s$ derivative of $f_{\sigma}(s,t)$ at $s=0$ vs $t$');
% xlabel('t');
% ylabel(['$\frac{\partial f_{\sigma}}{\partial s}(s=0,t)$']);
% dfsigma(0,1e-3)% output value at t=0 to check that it is negative

end
