function y = funiftemp(Au,di,df,T)% Uniform PDF of Schottky anomaly as a function of temperature T
    y = integral(@(D) fSchTemp(T,Au,D),di,df,'ArrayValued',true);
end
