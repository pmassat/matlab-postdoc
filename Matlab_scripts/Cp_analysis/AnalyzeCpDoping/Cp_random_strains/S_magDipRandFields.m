function S = S_magDipRandFields(t,F)
% t is the reduced temperature T/Tc
dT = diff(t);
diff1f = diff(F,1,1);% first order difference along 1st dimension (i.e. between rows) of F

[~,n] = size(F);% in case F is a multi-column array, e.g. calculated for various values of a parameter
S = zeros(size(diff1f));
for jh = 1:n
    S(:,jh) = -diff1f(:,jh)./dT';% entropy
end

end