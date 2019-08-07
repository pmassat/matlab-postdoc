function S = Entropy_TLFIM(t,F)
% t is the reduced temperature T/Tc
% F = zeros(size(t));
% for i=1:length(t)
%     F = freeEnergy_TLFIM(T(i),h,e);% free energy
% end
diff1f = diff(F,1,1);
dT = diff(t);
[~,n,p] = size(F);
S = zeros(size(diff1f));
for jh = 1:n
    for je = 1:p
        S(:,jh,je) = -diff1f(:,jh,je)./dT';% entropy
    end
end
end