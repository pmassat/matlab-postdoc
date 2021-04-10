%% This function is for finding the expectation value of an operator (X), such as O22
% with states given by the Hamiltonian Hcef at a temperature T
function [exvalue]=expect(X,Hcef,T)

%kb=1;
%exvalue=trace(expm(-1/(kb*T)*Hcef)*X);
[CEF,E]=eig(Hcef);
E=diag(E);
[OrderE,P]=sort(E);
OrderE = OrderE-ones(size(OrderE))*min(OrderE);
L=length(P);
%OrderE=OrderE ; %-OrderE(1);
%groundstate=OrderE(1);
Eigen=zeros(L);
for i=1:L
   Eigen(:,i)=CEF(:,P(i));
end

Z=sum(exp(-1/(T)*OrderE*(1))); %(11605) is 1 eV/ 1K in k_B units

Pe=exp(-1/(T)*OrderE*(1))/Z; %so 1 eV/1 K is 11605

exvalue=sum(Pe.*diag(Eigen'*X*Eigen));
%whos;

%Sum=0;
% for i=1:length(Pe)
%     Sum=Sum+Pe(i)*Eigen(:,i)'*X*Eigen(:,i);
% end
% exvalue=Sum;


end

