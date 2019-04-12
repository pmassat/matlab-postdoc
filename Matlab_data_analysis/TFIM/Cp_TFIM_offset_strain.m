function y = Cp_TFIM_offset_strain(t,e,h)
% t is the reduced temperature T/Tc
y = zeros(size(t));

if h==0; tc = 1;
elseif abs(h)<1; tc = h/atanh(h);
else; tc = 0;
end

for i=1:length(t)
    if t(i)<=0
        y(i) = 0;
    elseif t(i)<tc || e>0
        r = order_parameter_offset_strain(t(i),e)./t(i);% ratio of reduced order_parameter to reduced temperature
        y(i) = r^2.*sech(r+e)^2 ./ (1 - 1./t(i).*sech(r+e)^2);% mean-field heat capacity in the ordered phase
    else
        r = h./t(i);
        y(i) = r^2.*sech(r+e)^2;% mean-field heat capacity in the disordered phase
    end
end
end