function y = Cp_magDipRandFields(t,d0)
% t is the reduced temperature T/Tc
% High temperature correction to the free energy F from magnetic dipole
% interactions in the mean-field Ising model including random fields
% See equation (7) of paper by Gehring et al. (1976)
T = t';% transpose temperature array into column vector
F = zeros(length(T),length(d0));
for j1=1:length(F)
    for k=1:length(d0)
    f = @(t,u) 1/t * sech(u*d0(k)/t)^2 + 1/(u*d0(k)) * tanh(u*d0(k)/t);
    F(j1,k) = -1/(4*sqrt(pi))*integral(@(u)exp(-u.^2).*f(T(j1),u),-Inf,Inf,'ArrayValued',true);
    end
% If x(Tm) (i.e. concentration of Tm ion) is *NOT* at the denominator of y, then
% y is the correction per Tm ion; otherwise, it is the correction per
% rare-earth ion
end

dT = diff(T);
diff1f = diff(F,1,1);% first order difference along 1st dimension (i.e. between rows) of F
[~,m] = size(F);% in case F is a multi-column array, e.g. calculated for various values of a parameter
S = zeros(size(diff1f));
for j2 = 1:m% loop over all columns of F
    S(:,j2) = -diff1f(:,j2)./dT;% dF/dT = entropy
end

diff1s = diff(S,1,1);% first order difference between rows (1st dimension) of S
y = zeros(size(diff1s));
Tcp = T(2:end-1);
for j3 = 1:m% loop over all columns of S, which has the same number of columns as F
    y(Tcp<1,j3) = 0;
    y(Tcp>=1,j3) = Tcp(Tcp>=1).*diff1s(Tcp>=1,j3)./dT(Tcp>=1);% T*dS/dT = Cp
end

end
