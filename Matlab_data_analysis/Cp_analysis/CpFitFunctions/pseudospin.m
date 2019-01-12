function y = pseudospin(T,Ttrans,Tmf)
%%% calculate pseudospin as a function of temperature T from T=0 to T=Tmf
%%% for a given transition temperature Ttrans
    Tord = T(T<Tmf);
    y = 0.2*ones(1,length(Tord));% initialize pseudospin array
    g = @(t) fzero(@(x) Ttrans*2*x./atanh(2*x)-t,[1e-7 0.5]);
    for jj=1:length(y)
        y(jj) = g(Tord(jj));
    end
end
% The expression of the anonymous function inside fzero is that of a reduced 
% temperature defined as T_caluclated - T_measured, such that it takes the value 
% zero when T_calculated = T_measured. This expression results from inverting 
% the self-consistent equation 2*gamma = th(beta*gamma/2), where gamma is the 
% total pseudospin and beta is the inverse temperature beta = 1/(k_B*T). Therefore, 
% when the function takes the value zero, the correponding abscissa is the value 
% of the total pseudospin at the measured temperature.
% [1e-7 0.5] is the range of values taken by the pseudospin.
% 
% Note: the equation 2*gamma = th(beta*gamma/2) is for a reduced beta = beta/beta_c 
% such that k_B*T_c/J_0 = 1/4. In order to get the actual value of the critical 
% temperature, one needs to multiply k_B*T by 4*T_c. The factor 2 in the numerator 
% of the function results from the multiplication by 4 and the factor Tc is the 
% actual value of the transition temperature (estimated visually).

