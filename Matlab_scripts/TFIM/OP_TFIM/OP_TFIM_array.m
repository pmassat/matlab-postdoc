% function s = OP_TFIM(t,h,e)
function [s,opt] = OP_TFIM_array(t,h,e)
    s = zeros(length(t),length(h));
    opt = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
        'StepTolerance',1e-10, 'FunctionTolerance', 1e-6, 'Display', 'off');
    % LM algorithm is empirically better for the calculation without offset strain (e=0)
    for jh=1:length(h)
        gamma = @(x) sqrt((x+e).^2+h(jh).^2);
        for jt=1:length(t)
            fun = @(x) x-(x+e)./gamma(x).*tanh(gamma(x)./t(jt));
            s(jt,jh) = fsolve(fun,1,opt);
        end
    end
end   