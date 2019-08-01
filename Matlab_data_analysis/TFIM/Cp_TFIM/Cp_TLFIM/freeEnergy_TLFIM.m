function F = freeEnergy_TLFIM(t,h,e)
x = OP_TFIM(t,h,e);
gamma = sqrt((x+e)^2+h^2);
% bf = 0.5*beta*x^2 - log(2*cosh(beta*gamma));
F = 0.5*x^2 - t.*log(2*cosh(gamma./t));
end
