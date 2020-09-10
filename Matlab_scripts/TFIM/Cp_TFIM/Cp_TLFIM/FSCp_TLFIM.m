function [F,S,Cp] = FSCp_TLFIM(t,h,e)
    x = OP_TFIM_array(t,h,e);
    gamma = sqrt((x+e).^2 + h.^2);
    F = 0.5*x.^2 - t'.*log(2*cosh(gamma./t'));
    dt = diff(t');
    try
        S = -diff(F)./dt;% entropy
        dtm = 0.5*(dt(1:end-1)+dt(2:end));
        Cp = t(2:end-1)'.*diff(S)./dtm;
%     catch
%         warning(['Unable to compute entropy and heat capacity.'...
%             newline...
%             'Check dimensions of temperature array if that is not expected.'])
    end
end
