function y = fitCpunif(T,Ttrans,Tmf,Tsch,Amf,Au,di,df)
%%% Fitting function for heat capacity as a combination of: 
%%% a mean-field jump and *uniform* Schottky distribution below Tmf 
%%% a *uniform* Schottky distribution above Tsch
    y = zeros(size(T));
    gamma = pseudospin(T,Ttrans,Tmf);
    y(T<Tmf) = fmf(gamma,Amf) + funifps(Au,di,df,Ttrans,gamma);%
    y(T>Tsch) = funiftemp(Au,di,df,T(T>Tsch));
end
