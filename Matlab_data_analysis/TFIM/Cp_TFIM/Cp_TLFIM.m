 function y = Cp_TLFIM(t,S)
% t is the reduced temperature T/Tc
% S = zeros(size(t));
% for i=1:length(t)
%     S = Entropy_TLFIM(T(i),h,e);
% end
T = t';% transpose temperature array into column vector
diff1s = diff(S,1,1);% first order difference between rows of S
dT = diff(T);
[~,n,p] = size(S);
y = zeros(size(diff1s));
for jh = 1:n
    for je = 1:p
        y(:,jh,je) = T(2:end-1).*diff1s(:,jh,je)./dT(2:end);
    end
end
 end
