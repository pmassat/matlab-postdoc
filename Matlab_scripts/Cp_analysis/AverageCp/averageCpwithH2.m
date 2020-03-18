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
        Tm(k) = mean(T(abs(T-T(k))/T(k)<halfrelTsep));
        Tsd(k) = std(T(abs(T-T(k))/T(k)<halfrelTsep));
        Cpm(k) = mean(Cp(abs(T-T(k))/T(k)<halfrelTsep));
        Cpsd(k) = std(Cp(abs(T-T(k))/T(k)<halfrelTsep));
        CpmErr(k) = sum(CpErr(abs(T-T(k))/T(k)<halfrelTsep))/...
            sqrt(length(CpErr(abs(T-T(k))/T(k)<halfrelTsep)));
        T(abs(T-T(k))/T(k)<halfrelTsep)=0;
    else
        Tm(k) = mean(T(abs(T-T(k))/T(k)<relTsep));
        Tsd(k) = std(T(abs(T-T(k))/T(k)<relTsep));
        Cpm(k) = mean(Cp(abs(T-T(k))/T(k)<relTsep));
        Cpsd(k) = std(Cp(abs(T-T(k))/T(k)<relTsep));
        CpmErr(k) = sum(CpErr(abs(T-T(k))/T(k)<relTsep))/...
            sqrt(length(CpErr(abs(T-T(k))/T(k)<relTsep)));
        T(abs(T-T(k))/T(k)<relTsep)=0;
    end
end
S.H = H(Tm>0);
S.T = Tm(Tm>0);
S.Tsd = Tsd(Tm>0);
S.Cp = Cpm(Cpm>0);
S.CpFullErr = Cpsd(Cpm>0) + CpmErr(Cpm>0);
end







  