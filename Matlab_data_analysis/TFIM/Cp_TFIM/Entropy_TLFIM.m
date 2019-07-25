function S = Entropy_TLFIM(t,h,e)
% t is the reduced temperature T/Tc
F = zeros(size(t));
for i=1:length(t)
    F = freeEnergy_TLFIM(T(i),h,e);% free energy
end
diff1f = diff(F);
dT = diff(T);
S = -diff1f./dT;% entropy
end