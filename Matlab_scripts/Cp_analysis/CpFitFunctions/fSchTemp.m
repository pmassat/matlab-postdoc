function y = fSchTemp(T,A,D)% Schottky anomaly as a function of temperature T
    y = A*(D./(2*T)).^2.*sech(D./(2*T)).^2;% D is the total splitting of the doublet
end
