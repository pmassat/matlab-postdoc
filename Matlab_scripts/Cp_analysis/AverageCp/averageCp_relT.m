function S = averageCp_relT(relTsep,T,Cp,CpErr)
%% Average datapoints taken at any given temperature and field setpoint
% We want to compute the average of data points that are repetitions of the
% same measurement, i.e. with same temperature and field setpoints
Tsep = relTsep;% Data points taken within an interval of Tsep are considered to be measured at the same temperature setpoint
% 6mK is an empirical estimate of the dispersion of data points taken consecutively at a given temperature. 
    T2 = repmat(T,1);
    Tm = zeros(length(T),1);
    Tsd = zeros(length(T),1);
    Cpm = zeros(length(Cp),1);
    Cpsd = zeros(length(Cp),1);
    CpmErr = zeros(length(CpErr),1);
    for k = 1:length(T2)
        if T2(k)==0
            continue
        elseif length(T2(abs(T2/T2(k)-1)<Tsep))>3
            Tsep2 = Tsep/2;% reduce the temperature interval
%             T2(abs(T2/T2(k)-1)<Tsep2)%print out values of temperature which
%             verify the if statement
            Tm(k) = mean(T2(abs(T2/T2(k)-1)<Tsep2));
            Tsd(k) = std(T2(abs(T2/T2(k)-1)<Tsep2));
            Cpm(k) = mean(Cp(abs(T2/T2(k)-1)<Tsep2));
            Cpsd(k) = std(Cp(abs(T2/T2(k)-1)<Tsep2));
            CpmErr(k) = sum(CpErr(abs(T2/T2(k)-1)<Tsep2))/...
                sqrt(length(CpErr(abs(T2/T2(k)-1)<Tsep2)));
            T2(abs(T2/T2(k)-1)<Tsep2)=0;
        else
            Tm(k) = mean(T2(abs(T2/T2(k)-1)<Tsep));
            Tsd(k) = std(T2(abs(T2/T2(k)-1)<Tsep));
            Cpm(k) = mean(Cp(abs(T2/T2(k)-1)<Tsep));
            Cpsd(k) = std(Cp(abs(T2/T2(k)-1)<Tsep));
            CpmErr(k) = sum(CpErr(abs(T2/T2(k)-1)<Tsep))/...
                sqrt(length(CpErr(abs(T2/T2(k)-1)<Tsep)));
            T2(abs(T2/T2(k)-1)<Tsep)=0;
        end
    end
    S.T = Tm(Tm>0);
    S.Tsd = Tsd(Tm>0);
    S.Cp = Cpm(Cpm>0);
    S.CpFullErr = Cpsd(Cpm>0) + CpmErr(Cpm>0);
end







  