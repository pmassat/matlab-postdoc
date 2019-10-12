 function y = Cp_LFIM_CW2(t,e,t0,c)
% t is the reduced temperature T/Tc
y = zeros(size(t));
tc =1;% tc used to distinguish two cases: T/Tc<1 and >1

for i=1:length(t)
    y(i) = 0;
    x = OP_TFIM(t(i),0,e);
    gamma = (x+e);
    r = gamma./t(i);% ratio of reduced order_parameter to reduced temperature
    D = 1 - sech(r)^2./t(i);
    if t(i)<tc 
        y(i) = r^2.*sech(r)^2 ./ D;% mean-field heat capacity for non-zero order parameter
    elseif e~=0
%         y(i) = r^2.*sech(r)^2 ./ D + c/(t(i)-t0);% mean-field heat capacity + Curie-Weiss divergence above Tc
        y(i) = r^2.*sech(r)^2 ./ D + c*t(i)/(t(i)-t0)^2;% mean-field heat capacity + Curie-Weiss divergence above Tc
    else
%         y(i) = c/(t(i)-t0);% Curie-Weiss divergence above Tc
        y(i) = c*t(i)/(t(i)-t0)^2;% mean-field heat capacity + Curie-Weiss divergence above Tc
    end
end
end