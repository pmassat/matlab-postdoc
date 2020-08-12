function y = Cp_TFIM(t,h)
% t is the reduced temperature T/Tc
y = zeros(size(t));

if h==0; tc = 1;
elseif abs(h)<1; tc = h/atanh(h);
else; tc = 0;
end

for i=1:length(t)
    if t(i)<=0
        y(i) = 0;
    elseif t(i)<tc
%         r = order_parameter(t(i))./t(i);% ratio of reduced order_parameter to reduced temperature
        r = order_parameter(t(i))./t(i);% ratio of reduced order_parameter to reduced temperature
        y(i) = r.^2.*sech(r).^2 ./ (1 - 1./t(i).*sech(r).^2);% mean-field heat capacity in the ordered phase
    else
        r = h./t(i);
        y(i) = r.^2.*sech(r).^2;% mean-field heat capacity in the disordered phase
    end
end
end