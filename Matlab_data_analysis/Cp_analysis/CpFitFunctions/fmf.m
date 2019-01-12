function y = fmf(x,Amf)% Mean-field expression of heat capacity as a function of pseudospin x
    y = Amf*(atanh(2*x)).^2.*(1-4*x.^2)./(1-(1-4*x.^2).*atanh(2*x)./(2*x));
end
