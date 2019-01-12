function y = fnrmps(An,mu,sgm,Ttrans,x)% Normal PDF of Schottky anomaly as a function of pseudospin x
    y = integral(@(D)(fschps(x,Ttrans,An,D).*normpdf(D,mu,sgm)),0,inf,'ArrayValued',true);
end
