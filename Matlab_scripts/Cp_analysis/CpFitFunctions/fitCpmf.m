function y = fitCpmf(T,Ttrans,Tmf,Amf)
%%% Fitting mean-field heat capacity jump as a function of temperature below Tmf 
    y = zeros(size(T));
    gamma = pseudospin(T,Ttrans,Tmf);
    y(T<Tmf) = fmf(gamma,Amf);%
    y(T>Ttrans) = 0;
end
