function y = Schottky(A,D,T)% Schottky anomaly as a function of temperature T
    y = A*(D./(2*T)).^2.*sech(D./(2*T)).^2;% D is the total splitting of the doublet
end
% Note: this function was created as a duplicate of fSchTemp on 2020-07-09