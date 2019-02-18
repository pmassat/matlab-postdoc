function [H,eta] = convertVoigtToPseudo(Sigma,Gamma)
%     fwhmGsq = 8.0 * log(2) * SigmaSq;fwhmG = sqrt(fwhmGsq);
    fwhmG = 2*Sigma*sqrt(2*log(2));
    fwhmGsq = fwhmG^2;
    fwhmG4 = fwhmGsq^2;
    fwhmL = 2*Gamma;
    fwhmLsq = Gamma * Gamma;
    fwhmL4 = fwhmLsq * fwhmLsq;

    H = power(fwhmG^5 + 2.69269 * fwhmG^4 * fwhmL +...
          2.42843 * fwhmGsq * fwhmG^3 * fwhmL^2 +...
          4.47163 * fwhmGsq * fwhmLsq * fwhmL + 0.07842 * fwhmG * fwhmL4 +...
          fwhmL4 * fwhmL, 0.2);

    if H == 0
        H = epsilon * 1000.0;
    end

    tmp = fwhmL / H;

    eta = 1.36603 * tmp - 0.47719 * tmp * tmp + 0.11116 * tmp * tmp * tmp;
end
