function y = fschps(x,Ttrans,A,D)% Schottky anomaly as a function of pseudospin x
    y = A*(D./(2*Ttrans*2*x./atanh(2*x))).^2.*sech(D./(2*Ttrans*2*x./atanh(2*x))).^2 ;
end
