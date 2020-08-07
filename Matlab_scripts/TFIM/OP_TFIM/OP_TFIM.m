% function s = OP_TFIM(t,h,e)
function [s,opt] = OP_TFIM(t,h,e)
    gamma = @(x) sqrt((x+e).^2+h.^2);
    opt = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
        'StepTolerance',1e-8);
    % LM algorithm is empirically better for the calculation without offset strain (e=0)
    fun = @(x) x-(x+e)./gamma(x).*tanh(gamma(x)./t);
%     fun = @(x) sqrt(x^2)-tanh((x+e)./(t));
    s = fsolve(fun,1,opt);
end   