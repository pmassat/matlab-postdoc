function r = OP_TFIM(t,h,e)
    if h==0; tc = 1;
    elseif abs(h)<1; tc = h/atanh(h);
    else; tc = 0;
    end
    fun = @(x) sqrt(x^2+h^2)-tanh((sqrt((x+e)^2+h^2))./t);
%     fun = @(x) sqrt(x^2)-tanh((x+e)./(t));
    r = fsolve(fun,1);
    s = sqrt(r^2-(h)^2);
end