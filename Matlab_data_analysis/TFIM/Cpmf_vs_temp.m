function y = Cpmf_vs_temp(t)
% t is the reduced temperature T/Tc
y = zeros(size(t));
for i=1:length(t)
    if t(i)<=0
        y(i) = 0;
    elseif t(i)<1
    r = order_parameter(t(i))./t(i);% ratio of reduced order_parameter to reduced temperature
    y(i) = r^2.*sech(r)^2 ./ (1 - 1./t(i).*sech(r)^2);% mean-field heat capacity in the ordered phase
    else
        y(i) = 0;% heat capacity is zero above Tc
    end
end
end