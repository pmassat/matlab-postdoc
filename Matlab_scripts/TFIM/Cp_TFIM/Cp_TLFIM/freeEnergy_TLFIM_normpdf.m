function y = freeEnergy_TLFIM_normpdf(t,hm,e,sgm)% Free energy in the case when there is a normal PDF of fields in the sample
    y = integral(@(x)freeEnergy_TLFIM(t,x,e).*normpdf(x,hm,sgm),0,10,'ArrayValued',true);
% Integral is performed only up to x=10 because the computation results in
% NaN when computing up to values of order 100 or more. However, for hm<~1,
% which is the case that we are interested in,
% the result of the integral is unchanged when computing up to x>=5; hence
% xmax=10 is conservative enough.
end
