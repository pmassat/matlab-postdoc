function y = funifps(Au,di,df,Ttrans,x)% Uniform PDF of Schottky anomaly as a function of temperature T
    y = integral(@(D) fschps(x,Ttrans,Au,D),di,df,'ArrayValued',true);
end
