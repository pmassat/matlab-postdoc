function S = averageCpwithH2(relTsep,T,Cp,CpErr,H)
%% Average datapoints taken at any given temperature and field setpoint
% We want to compute the average of data points that are repetitions of the
% same measurement, i.e. with same temperature and field setpoints
% relTsep;% temperature stability reported in the PPMS manual is <0.2\%
% (should be called temperature instability). Hence data points taken 
% within a relative interval of relTsep are considered to be measured at the same temperature 
% relTsep <= 0.002 according to the PPMS manual
Tm = zeros(size(T));
Tsd = zeros(size(T));
Cpm = zeros(size(Cp));
Cpsd = zeros(size(Cp));
CpmErr = zeros(size(CpErr));
for k=1:length(T)
    if T(k)==0
        continue
    elseif length(T(abs(T-T(k))/T(k)<relTsep))>3
% T(abs(T-T(k))/T(k)<relTsep) is the subset of temperatures which have
% a relative separation of relTsep from the temperature of datapoint #k
% If this subset contains more than 3 data points...
        halfrelTsep = relTsep/2;% ... reduce the temperature interval
%             T(abs(T-T(k))<Tsep2)%print out values of temperature which
%             verify the if statement
        Tselect1 = abs(T-T(k))/T(k)<halfrelTsep;
        Tm(k) = mean(T(Tselect1));
        Tsd(k) = std(T(Tselect1));
        Cpm(k) = mean(Cp(Tselect1));
        Cpsd(k) = std(Cp(Tselect1));
        CpmErr(k) = sum(CpErr(Tselect1))/...
            sqrt(length(CpErr(Tselect1)));
        T(Tselect1)=0;
    else
        Tselect2 = abs(T-T(k))/T(k)<relTsep;
        Tm(k) = mean(T(Tselect2));
        Tsd(k) = std(T(Tselect2));
        Cpm(k) = mean(Cp(Tselect2));
        Cpsd(k) = std(Cp(Tselect2));
        CpmErr(k) = sum(CpErr(Tselect2))/...
            sqrt(length(CpErr(Tselect2)));
        T(Tselect2)=0;
    end
end
S.H = H(Tm>0);
S.T = Tm(Tm>0);
S.Tsd = Tsd(Tm>0);
S.Cp = Cpm(Tm>0);
S.CpFullErr = Cpsd(Tm>0) + CpmErr(Tm>0);
end







  