function dsz = order_parameter_derivative_random_strains_offset(delta0,ea,temp,sz)
dsgm = @(ds) ds-1/sqrt(pi)*integral(@(u)(ds./temp-(sz+u.*delta0)./(temp.^2))...
    .*exp(-(u-ea).^2)./cosh((sz+u.*delta0)./temp).^2,-Inf,Inf);
dsz = fzero(@(ds)dsgm(ds),[-30 0]);% lowest empirical value is ~-27

%% Plot ds
% % fplot is slow; it is faster to compute an array and use plot
% figure; hold on;
% fplot(@(s)ds(s,tplot))
% line(xlim,[0 0],'color','black','linestyle','--')

%% Plot dsz
% % fplot is slow; it is faster to compute an array and use plot
% figure
% fplot(@(t)dsz(t),[1e-2 tc-1e-3])
% title(sprintf('Derivative of the order parameter vs temperature at $x$=%.2f',1-dpg(i)));
% xlabel('$t=\frac{T_D}{T_D(x=1)}$');
% ylabel('$\frac{d\left<S^{z}\right>}{dt}$');

end