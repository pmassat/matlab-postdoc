 function y = Cp_LFIM_NF(T,Tc,e,c1,c2)
t = T/Tc;% t is the reduced temperature T/Tc
y = zeros(size(t));
tc = 1;% tc used to distinguish two cases: T/Tc<1 and >1
T66 = 1.85/Tc;% temperature of divergence of S66=1/C66 reported by Melcher 1976
T_C = 0;% Curie temperature for a "para-nematic" material

for i=1:length(t)
    y(i) = 0;
    x = OP_TFIM(t(i),0,e);
    gamma = (x+e);
    r = gamma./t(i);% ratio of reduced order_parameter to reduced temperature
    D = 1 - sech(r)^2./t(i);
    if t(i)<tc 
        y(i) = r^2.*sech(r)^2 ./ D;% mean-field heat capacity for non-zero order parameter
    elseif e~=0
        y(i) = r^2.*sech(r)^2 ./ D + Cp_nem_fluct(c1,c2,t(i),T66,T_C);% mean-field heat capacity
    elseif t(i)>tc 
        y(i) = Cp_nem_fluct(c1,c2,t(i),T66,T_C);% contribution of nematic fluctuations to Cp above Tc
    end
end
end