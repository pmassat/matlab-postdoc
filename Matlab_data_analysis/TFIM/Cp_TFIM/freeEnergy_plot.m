h = 1e-3;% transverse field
e = 0;% longitudinal field
T = linspace(1e-3,1-1e-3,998);
F = repmat(T,1);
for i=1:length(T)
    F(i) = betaFreeEnergy_TFIM(T(i),h,e).*T(i);% free energy
end
%%
figure;
plot(T,F)
title('Free energy')
%%
beta = 1./T;%
diff1beta = diff(beta);
diff1f = diff(F);
d1betabf = diff1f./diff1beta;
dT = diff(T);
% d1bf = - beta(2:end).^2 .* d1betabf;% df/dT = -1/beta^2 * df/dbeta
S = -diff1f./dT;% entropy
%%
figure;
plot(T(2:end),S)
title('Entropy')
%%
diff2f = diff(S);
Cp = T(3:end).*diff2f./dT(2:end);
%% 
figure;
plot(T(3:end),Cp)
title('Cp')
