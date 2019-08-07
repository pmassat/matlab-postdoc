function y = fnrmtemp_sym(mu,sgm,T)% Normal PDF of Schottky anomaly as a function of temperature T
% It is symmetrized wrt the energy gap, i.e. integrated from -inf to inf
    y = integral(@(D)fSchTemp(T,1,2*D).*normpdf(D,mu,sgm),-Inf,inf,'ArrayValued',true);
% Here we use 2*D for consistency with function 'Cp_full_random_strains' 
% where D is defined as half of the splitting of the doublet 
end
