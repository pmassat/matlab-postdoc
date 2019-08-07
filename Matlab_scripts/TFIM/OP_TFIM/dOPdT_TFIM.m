function dg = dOPdT_TFIM(t,h,e)
% Temperature derivative of the order parameter in the TFIM
% Does not work for h>0...
% Correction needed: the current code computes dgamma/dT, instead of dOP/dT
    x = OP_TFIM(t,h,e);
    gamma = sqrt((x+e)^2+h^2);
    r = gamma./t;
    dg = t^2*gamma./(1/t-cosh(r).^2);
end
