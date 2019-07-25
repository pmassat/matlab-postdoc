 function y = Cp_TLFIM(t,h,e)
% t is the reduced temperature T/Tc
S = zeros(size(t));
for i=1:length(t)
    S = Entropy_TLFIM(T(i),h,e);
end
dT = diff(T);
y = t(2:end-1).*diff(S)./dT(2:end);
end