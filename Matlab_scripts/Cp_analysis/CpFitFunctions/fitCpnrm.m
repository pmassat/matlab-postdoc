function y = fitCpnrm(T,Ttrans,Tmf,Tsch,Amf,An,mu,sgm)
%%% Fitting function for heat capacity as a combination of : 
%%% a mean-field jump and *normal* Schottky distribution below Tmf 
%%% a *normal* Schottky distribution above Tsch
    y = zeros(size(T));
    gamma = pseudospin(T,Ttrans,Tmf);
    y(T<Tmf) = fmf(gamma,Amf) + fnrmps(An,mu,sgm,Ttrans,gamma);%
    y(T>Tsch) = fnrmtemp(An,mu,sgm,T(T>Tsch));
end
