%% Analyze heat capacity from DR
% This routine is intended at anaLzing Cp data acquired with our DR used in 
% the Dynacool PPMS of the Lee lab
%% To do: 
% * Error bars:
% *     For dT data, use standard deviation
% *     For dCp values, add the SampHCErrmJmoleK field 
%% Import data YTmVO4
Data = ImportCpDR('2018-07-31_TmVO4-LS5228-DR-HC180731.dat');
% Data0 = importCpSharedPPMS_('20180322_TmVO4-LS5228-MP3-Plt-HC1803_Cp.dat');
Data10 = ImportCpDR('DRHC_10YTmVO4-LS5251-DR-HC180511.dat');
Data15 = ImportCpDR('DRHC_15Y-TmVO4-LS5343-DR-HC180721.dat');
Data18 = ImportCpDR('2018-07-24_DR-HC_18Y-TmVO4-LS5348-DR-HC180724.dat');
Data20 = ImportCpDR('DRHC_20YTmVO4-LS5252-DR-HC180519.dat');
Data22 = ImportCpDR('2018-07-23_DR-HC_22Y-TmVO4-LS5344-DR-HC180723.dat');
Data24 = ImportCpDR('2018-07-25_DR-HC_24Y-TmVO4-LS5349-DR-HC180725.dat');
Data30 = ImportCpDR('DRHC_30YTmVO4-LS5335-DR-HC180526.dat');
%% Concatenate them in a cell array
%%
split = {Data; Data10; Data15; Data18; Data20; Data22; Data24; Data30};
L = length(split);
%% Rename variables
%%
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
% function to test if a column exists in a table
for i=1:L% for doped samples measured in DR
    if isTableCol(split{i},'SampleTempKelvin')
        split{i}.Properties.VariableNames{'SampleTempKelvin'} = 'T';% rename the temperature column
    end
    if isTableCol(split{i},'FieldOersted')
        split{i}.Properties.VariableNames{'FieldOersted'} = 'H';% rename the magnetic field column
    end
    if isTableCol(split{i},'SampHCJK')% for samples measured in DR
        split{i}.Properties.VariableNames{'SampHCJK'} = 'Cp';% rename the heat capacity column
    elseif isTableCol(split{i},'SampHCmJmoleK')% for samples measured in shared PPMS
        split{i}.Properties.VariableNames{'SampHCmJmoleK'} = 'Cp';
    elseif isTableCol(split{i},'SampHCJmoleK')% for samples measured in Fisher He4 PPMS
        split{i}.Properties.VariableNames{'SampHCmJmoleK'} = 'Cp';
    end
end
%% Remove NaN rows
%%
for i=1:L
    split{i}(any(isnan(split{i}.T), 2), :) = [];% Remove rows where T is NaN
end
%% Keep only data under zero magnetic field
%%
for i=1:L%for all datasets
    split{i} = split{i}(round(split{i}.H,-1)==0,:);% keep only data at zero field
end
%% Plot parameters for heat capacity
%%
xlblTemp = 'Temperature (K)';
ylblCp = 'C_p (J\cdotmol^{-1}\cdotK^{-1})';
ttlCpY ='Heat capacity of Tm_{1-x}Y_xVO4';
%% Compute molar heat capacity
%%
M = [283.87326% Molar mass of each sample, in g/mol
275.3902538
271.869006
269.4681552
266.3470492
266.2670208
264.6664536
258.6643266];
m = 1e-3*[0.38% mass in g of the sample on which each dataset was measured
0.13
0.24
0.045
0.25
0.59
0.76
0.14];
dpg = 1e-2*[0% Y content for each dataset
10.6
15.3%
18.1%
21.9%
22.0%
24.8%
31.5];
% split{1}.Cpmol = split{1}.Cp *1e-3;% molar heat capacity, in J/mol/K,
% starting from a heat capacity in mJ/mol/K as is measured in the shared PPMS
for i=1:L
    split{i}.Cpmol = split{i}.Cp *1e-6*M(i)/(m(i)*(1-dpg(i)));% molar heat capacity, in J/mol/K
    % starting from a heat capacity measured in microJoules per Kelvin, as is measured in the DR
    % Cpmol is calculated per mole of Tm3+ ions, hence the (1-dpg) factor in the denominator
end
%% Plot the dataset for YTmVO4
%%
plotCpDoping(split,dpg,1,L,ttlCpY)
%% 
% 
% 
% 
% 
% 
% 
% 
%% Average data
%% Sort each dataset by increasing value of temperature
%%
for i=1:L
    srtd{i} = sortrows(split{i},{'T'});
%     [split{i}.T,wo] = sort(split{i}.T);
%     split{i}.Cp = split{i}.Cp(wo);
end
srtd = srtd';
%% Test the temperature scattering between data points supposedL taken at the "same" temperature
% 4mK is the empirical maximal dispersion of data points taken consecutiveL 
% at a given temperature. In order to check this, one can run the following code 
% on the temperature variable T:

nrep = 3;% number of repetitions of the measurement at each temperature
% We want to compute the temperature separation between two data points
% that are NOT repetitions of the same measurement
Tsep = 4e-3;% value of temperature separation to be tested, in Kelvin
% Data points taken within an interval of Tsep are considered to be measured at the same temperature setpoint
for i = 1
    for k = 1:length(srtd{i}.T)-nrep
        if abs(srtd{i}.T(k+nrep)-srtd{i}.T(k))<Tsep
% if the values of temperature for two data points
% measured at two different temperature setpoints are separated by less than Tsep
            srtd{i}.T(k);% display the value of temperature for the first data point
        end
    end
end
%% 
% If this piece of code outputs one or more values AND IF these values obviousL 
% correspond to measurements that were supposed to be measured at different temperature 
% setpoints, THEN reduce Tsep in order to reach a value for which the code does 
% not output anything.
%% Compute average of data points taken
%%
for i = 1:L
    Tm = [];Cpm = [];stdCpm = [];% initialize temporary tables;
    k = 1;% initiliaze loop index k
    while k<length(srtd{i}.T)% loop over k
        ind = [k];% initialize temporary table of indices
        j = k+1;% initialize subloop index 
        while abs(srtd{i}.T(j)-srtd{i}.T(k))<Tsep
            % loop over datapoints within a temperature interval of Tsep
            % i.e. datapoints measured at the same temperature setpoint
            ind = [ind, j];% store index of current datapoint
            j = j+1;% increase subloop index
            if j>= length(srtd{i}.T)% safety condition
                break
            end
        end
        if length(ind)>1% do not include isolated data points
            Tm = [Tm,mean(srtd{i}.T(ind))];% average of temperature for data points
            % taken at the same temperature setpoint
            Cpm = [Cpm,mean(srtd{i}.Cpmol(ind))];% average of heat capacity
            stdCpm = [stdCpm,std(srtd{i}.Cpmol(ind))];% standard deviation of heat capacity
        end
        k = j+1;% once all datapoints measured at the same temperature setpoints are averaged,
        % increase the main loop index in order to jump to the next data point measured
        % at a different temperature setpoint
    end
    avgData(i).T = Tm;% store the full averaged temperature table
    avgData(i).Cp = Cpm;% store the full averaged heat capacity table
    avgData(i).stdCp = stdCpm;% store the full table of standard deviation of heat capacity
end
%% Display name of structure containing average data
%%
disp("Structure containing averaged data is called avgData")
% 
% 
% 
% 
%% Export averaged Cp data to text file
% savepath = 'C:\Users\Pierre\Desktop\Postdoc\YTmVO4\YTmVO4_HeatCapacity\YTmVO4_Cp_anaLsis\2018-10-17_YTmVO4_averaged_Cp\';
% for i=L-3:L
%     dpg(i);
%     expstr(i) = fullfile(savepath, ['AnaLzeCpDRdoping_2018-11-30_YTmVO4_x=',num2str(dpg(i)),'.txt']);
%     fid = fopen(expstr(i), 'wt');
%     fprintf(fid, '%s\n', ['Exported from AnaLzeCpDRdoping_2018-10-10_YTmVO4_VTmAsO4.mlx on ',date]);  % header
%     fclose(fid);
%     dlmwrite(expstr(i),cat(2,avgData(i).T',avgData(i).Cp',avgData(i).stdCp'),'-append','delimiter','\t')
% end
% 