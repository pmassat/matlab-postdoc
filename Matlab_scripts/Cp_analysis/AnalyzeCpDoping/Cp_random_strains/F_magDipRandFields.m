function F_magDipRandFields(T)
% High temperature correction to the free energy F from magnetic dipole
% interactions in the mean-field Ising model including random fields
% See equation (7) of paper by Gehring et al. (1976)
f = @(t,u) 1/t * sech(u*d0/t)^2 + 1/(u*d0) * tanh(u*d0/t);
y = zeros(size(T));
for j=1:length(y)
    y(j) = 1/(4*sqrt(pi))*integral(@(u)exp(-u.^2).*f(T(j),u),-Inf,Inf,'ArrayValued',true);
% If x(Tm) (i.e. concentration of Tm ion) is *NOT* at the denominator of y, then
% y is the correction per Tm ion; otherwise, it is the correction per
% rare-earth ion
end
end