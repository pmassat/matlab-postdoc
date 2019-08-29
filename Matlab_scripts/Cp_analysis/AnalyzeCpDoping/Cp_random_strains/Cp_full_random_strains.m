function [Cp,sz,dsz] = Cp_full_random_strains(d0,da,t)
% % % Compute integral of molar heat capacity
% See Gehring et al. 1976
% t is an array of reduced temperature T/Tc(x=1)
% d0 is the width of the gaussian strain distribution
% sz is an array containing values of the order parameter, with same size as t
% dsz is an array of dsz/dt, with same size as t
% x is the concentration of the JT-active ion(Tm); x=1 for pure TmVO4
sz = zeros(size(t));
dsz = zeros(size(t));
Cp = zeros(size(t));
tc = random_strains_phase_boundary(d0,da);
dt = 1e-2;

for k=1:length(t)
    if t(k)<tc-dt || da>0
    sz(k) = order_parameter_random_strains_offset(d0,da,t(k));
    try
    dsz(k) = order_parameter_derivative_random_strains_offset(d0,da,t(k),sz(k));
    catch
        disp(k);% the following commands do not solve the errors: dt = tc-t(k); k = k-1;% continue;
    end
    %% Compute expression in integrand of Cpm
    Er = @(u) (sz(k)+u.*d0)./t(k);

    Cpintegrand = @(u) -(dsz(k)./2.*tanh(Er(u))+...
        (sz(k)/2 + u*d0).*(dsz(k)-Er(u))./t(k)./cosh(Er(u)).^2);

    Cp(k) = 1/sqrt(pi)*integral(@(u)exp(-(u-da).^2).*Cpintegrand(u),-Inf,Inf,'ArrayValued',true);
% When x is *NOT* at the numerator of Cp (i.e. Cp is divided by x), 
% it is the heat capacity data *per Tm ion*; otherwise it is the raw Cp
    end
%     if t(k)>=tc-dt
%         Cp(k) = Cp(k) + fnrmtemp_sym(da,d0,t(k));
%     end
end

end



% 
% 
% 
% 
% 
%

