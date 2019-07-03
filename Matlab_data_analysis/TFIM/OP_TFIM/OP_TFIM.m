function s = OP_TFIM(t,h,e)
%     if h==0; tc = 1;
%     elseif abs(h)<1; tc = h/atanh(h);
%     else; tc = 0;
%     end
    gamma = @(x) sqrt((x+e)^2+h^2);
    options = optimoptions('fsolve','Algorithm','levenberg-marquardt');% empirically better for the calculation without offset strain (e=0)
    fun = @(x) x-(x+e)/gamma(x)*tanh(gamma(x)/t);
%     fun = @(x) sqrt(x^2)-tanh((x+e)./(t));
    s = fsolve(fun,1,options);
end