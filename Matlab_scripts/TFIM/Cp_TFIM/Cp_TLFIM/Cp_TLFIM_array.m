function y = Cp_TLFIM_array(t,S)
% t is the reduced temperature T/Tc
T = t';% transpose temperature array into column vector
dT = diff(T);
diff1s = diff(S,1,1);% first order difference between rows of S

% Note: this subsection would work even if p=1
% [~,n,p] = size(S);
% y = zeros(size(diff1s));
% for jh = 1:n
%     for je = 1:p
%         y(:,jh,je) = T(2:end-1).*diff1s(:,jh,je)./dT(2:end);
%     end
% end

[~,n] = size(S);
y = zeros(size(diff1s));
for jh = 1:n
    y(:,jh) = T(2:end-1).*diff1s(:,jh)./dT(2:end);
end

end
