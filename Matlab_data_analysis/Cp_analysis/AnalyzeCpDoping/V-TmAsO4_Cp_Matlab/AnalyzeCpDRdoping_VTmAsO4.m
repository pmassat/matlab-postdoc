%% Analyze heat capacity from DR
% This routine is intended at analyzing Cp data acquired with our DR used in 
% the Dynacool PPMS of the Lee lab
%% To do: 
% * Error bars:
% *     For dT data, use standard deviation
% *     For dCp values, add the SampHCErrmJmoleK field 
%% Import data YTmVO4

AsData0=ImportCpFisherPPMSTable('TmAsO4-LS5341-HC180602_2018-06-02.dat');
AsData10 = importCpSharedPPMS_('2018-09-07_10V-TmAsO4-LS5350-HC180907_CpvsT.dat');
%% Concatenate them in a cell array
%%
splitAs = {AsData0;AsData10};
Las = length(splitAs);
%% Rename variables
%%
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
% function to test if a column exists in a table
for i=1:Las% for doped samples measured in DR
    if isTableCol(splitAs{i},'SampleTempKelvin')
        splitAs{i}.Properties.VariableNames{'SampleTempKelvin'} = 'T';% rename the temperature column
    end
    if isTableCol(splitAs{i},'FieldOersted')
        splitAs{i}.Properties.VariableNames{'FieldOersted'} = 'H';% rename the magnetic field column
    end
    if isTableCol(splitAs{i},'SampHCJK')% for samples measured in DR
        splitAs{i}.Properties.VariableNames{'SampHCJK'} = 'Cp';% rename the heat capacity column
    elseif isTableCol(splitAs{i},'SampHCmJmoleK')% for samples measured in shared PPMS
        splitAs{i}.Properties.VariableNames{'SampHCmJmoleK'} = 'Cp';
    elseif isTableCol(splitAs{i},'SampHCJmoleK')% for samples measured in Fisher He4 PPMS
        splitAs{i}.Properties.VariableNames{'SampHCmJmoleK'} = 'Cp';
    end
    if isTableCol(splitAs{i},'SampHCErrJmoleK')% for samples measured in DR
        splitAs{i}.Properties.VariableNames{'SampHCErrJmoleK'} = 'CpErr';
    elseif isTableCol(splitAs{i},'SampHCErrJK')% for samples measured in DR
        splitAs{i}.Properties.VariableNames{'SampHCErrJK'} = 'CpErr';        
    elseif isTableCol(splitAs{i},'SampHCErrmJmoleK')% for samples measured in DR
        splitAs{i}.Properties.VariableNames{'SampHCErrmJmoleK'} = 'CpErr';        
    end
end
%% Remove NaN rows
%%
for i=1:Las
    splitAs{i}(any(isnan(splitAs{i}.T), 2), :) = [];% Remove rows where T is NaN
end
%% Keep only data under zero magnetic field
%%
for i=1:Las%for all data but the last dataset
    splitAs{i} = splitAs{i}(round(splitAs{i}.H,-1)==0,:);% keep only data at zero field
end
% splitAs{Las} = splitAs{Las}(round(splitAs{Las}.H,-1)==20000,:);% keep only data at 20000 Oe for the second TmAsO4 dataset
%% Plot parameters for heat capacity
%%
xlblTemp = 'Temperature (K)';
ylblCp = 'C$_p$ (J$\cdot$mol$^{-1}\cdot$K$^{-1}$)';
ttlCpAs ='Heat capacity of TmAs$_{1-x}$V$_x$O$_4$';
%% Compute molar heat capacity
%%
Mas = [
307.8534% Molar mass of each sample, in g/mol
305.45539];
mAs = 1e-3*[% mass in g of the sample on which each dataset was measured
1.23
0.35%This value was adjusted (instead of 0.32mg) in order to let the
];
% curve fall on top of that of the parent compound
dpgAs = 1e-2*[0% This is V content in Tm(As,V)O4
10.80];
for i=1:Las
    splitAs{i}.Cpmol = splitAs{i}.Cp *1e-6*Mas(i)/mAs(i);% molar heat capacity, in J/mol/K
    splitAs{i}.CpmolErr = splitAs{i}.CpErr *1e-6*M(i)/(m(i)*(1-dpg(i)));% molar heat capacity, in J/mol/K
    % starting from a heat capacity measured in microJoules per Kelvin, as is measured in the DR
    % Cpmol is calculated per mole of Tm3+ ions, hence the (1-dpg) factor in the denominator
end
%% Plot the dataset for TmAsVO4
%%
plotCpDoping(splitAs,dpgAs,1,Las,ttlCpAs)
xlim([0 inf]);ylim([0 inf]);

%% Average data
%% Sort each dataset by increasing value of temperature
%%
for i=1:Las
    srtdAs{i} = sortrows(splitAs{i},{'T'});
%     [split{i}.T,wo] = sort(split{i}.T);
%     split{i}.Cp = split{i}.Cp(wo);
end
srtdAs = srtdAs';

%% Compute average of data points taken
for i = 1:Las
    avgDataAs(i) = averageCp(6e-3,srtdAs{i}.T,srtdAs{i}.Cpmol,srtdAs{i}.CpmolErr);
end

%% Plot averaged data
figure
for i=1:Las
    errorbar(avgDataAs(i).T,avgDataAs(i).Cp,avgDataAs(i).CpFullErr,'.','MarkerSize',18,'DisplayName',['x = ',num2str(dpgAs(i))])
    hold on
end
xlabel(xlblTemp); ylabel(ylblCp);
title(ttlCpAs);
legend('show');

%% Plot Cp data vs reduced temperature
TmaxAs = [6.04; 5.3];% temperature of the maximum of the Cp jump
figure; %
co = get(gca,'ColorOrder');
for i=1:2
    avgDataAs(i).t = (avgDataAs(i).T-TmaxAs(i))/TmaxAs(i);
    avgDataAs(i).tp = avgDataAs(i).t(avgDataAs(i).t>0);
    avgDataAs(i).tm = avgDataAs(i).t(avgDataAs(i).t<0);
    semilogx(avgDataAs(i).tp,avgDataAs(i).Cp(avgDataAs(i).t>0),'.',...
        'Color',co(i,:),'DisplayName',sprintf('x=%.2f t$>$0',dpg(i)));
    hold on;
    semilogx(-avgDataAs(i).tm,avgDataAs(i).Cp(avgDataAs(i).t<0),'x',...
        'MarkerSize',8,'Color',co(i,:),'DisplayName',sprintf('x=%.2f t$<$0',dpg(i)));
end
legend('show','Location','best')
