 function y = Cp_LFIM_MEC_cst_stress(t,sigma,a,t0)
% t is the reduced temperature T/Tc
y = zeros(size(t));
tc =1;

for i=1:length(t)
    y(i) = 0;
    x = OP_TFIM(t(i),0,sigma);
    gamma = (x+sigma);
    r = gamma./t(i);% ratio of reduced order_parameter to reduced temperature
    D = 1 - sech(r)^2./t(i);
    if t(i)<tc
        y(i) = r^2.*sech(r)^2 ./ D;% mean-field heat capacity in the ordered phase
    elseif sigma~=0
        y(i) = r^2.*sech(r)^2 ./ D + a*t(i)/(t(i)-t0)^3;
    else
        y(i) = a*t(i)/(t(i)-t0)^3;
    end
end
end