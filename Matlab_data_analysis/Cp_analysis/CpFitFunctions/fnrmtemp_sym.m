function y = fnrmtemp_sym(An,mu,sgm,T)% Normal PDF of Schottky anomaly as a function of temperature T
    y = integral(@(D)fSchTemp(T,An,D).*normpdf(D,mu,sgm),-Inf,inf,'ArrayValued',true);
end
