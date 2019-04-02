function y = Cpm_random_strains(d0,t,sz,dsz,x)
% % % Compute integral of molar heat capacity
% See Gehring et al. 1976
% t is an array of reduced temperature T/Tc(x=1)
% d0 is the width of the gaussian strain distribution
% sz is an array containing values of the order parameter, with same size as t
% dsz is an array of dsz/dt, with same size as t
% x is the concentration of the JT-active ion(Tm); x=1 for pure TmVO4
y = zeros(size(t));

% Er = cell(size(t));
% Cpintegrand = cell(size(t));
for k=1:length(t)
    %% Compute expression in integrand of Cpm
    Er = @(u) (x.*sz(k)+u.*d0)./t(k);

    Cpintegrand = @(u) -(x/2.*dsz(k).*tanh(Er(u))+...
        (x.*sz(k)/2 + u*d0).*(x.*dsz(k)-Er(u))./t(k)./cosh(Er(u)).^2);

    y(k) = x/sqrt(pi)*integral(@(u)exp(-u.^2).*Cpintegrand(u),-Inf,Inf,'ArrayValued',true);
    % When Cpm is divided by x, it compares data *per Tm ion*
end

end



% 
% 
% 
% 
% 
%

