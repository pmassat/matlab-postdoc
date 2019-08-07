function y = Cpmf(x,Tc)
% x is the reduced pseudospin 2*Sz/J0, hence Cpmf is defined on the range [0 1]
%     beta = (2./x).*atanh(2.*x);% beta = 4./temperature 
%     f = 1/4*(sech(beta.*x./2)).^2;% factor appearing in the expression of 
%     f = (sech(x*Tc/T))^2;
%     y = (x.*beta).^2.*f./(1-beta.*f);% Cp in the ordered phase
    T = Tc*x/atanh(x);% temperature as a function of 
    r = x*Tc/T;% ratio of reduced field to reduced temperature
    y = r^2*sech(r)^2 / (1 - Tc/T*sech(r)^2);
end