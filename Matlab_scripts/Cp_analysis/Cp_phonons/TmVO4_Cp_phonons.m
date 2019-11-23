%% Import data from directory containing it
cd C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2019-03_TmVO4-LS5216-NMR-FIB\2019-03_TmVO4-LS5216-NMR-FIB_analysis;
data(1).tbl = Import_Cp_avg('2019-03-25_TmVO4-LS5216-NMR-FIB-01_avg.txt');
data(2).tbl = Import_Cp_avg('2019-04-04_TmVO4-LS5214-NMR-FIB-02-HC1_avg.txt');
data(3).tbl = Import_Cp_avg('2019-03-29_TmVO4-LS5216-NMR-FIB-03-HC1_avg.txt');

%% Import data
cd C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2019-05-08_TmVO4-LS5210-NMR1905-1-HC1\TmVO4-LS5210-NMR1905-1-HC1_analysis;
data(4).tbl = Import_Cp_avg('2019-05-08_TmVO4-LS5210-NMR1905-1-HC1_avg.txt');

%% Define variables for CFT (Curve Fitting Tool)
R = 8.314;% gas constant, in J/mol/K
for i=1:4
% Important note: the following syntax is #1 in the list of things that one
% should NOT do in Matlab, HOWEVER I do not see any other iterative way of
% creating variables that can be accessed in the CFT.
eval(['fitT' num2str(i) '= data('  num2str(i) ').tbl.T']);
eval(['fitCp' num2str(i) '= data('  num2str(i) ').tbl.Cp /R']);
eval(['fitWghts' num2str(i) '= R./data('  num2str(i) ').tbl.Cp_Err']);
end

%% Phonons fit coefficients before optimizing 
% Fit formula: c0 + (T/Tp)^3
% Fit includes weights
% fitCp1;  Adjusted R-square: 0.9892
Tp1 =       26.86;%  (26.32, 27.4)
c01 =    0.007943;%  (0.005699, 0.01019)
% fitCp2;  Adjusted R-square: 0.8799
Tp2 =       28.77;%  (25.09, 32.45)
c02 =    -0.01208;%  (-0.02855, 0.004376)
% fitCp3;  Adjusted R-square: 0.969
Tp3 =       29.68;%  (28.2, 31.16)
c03 =    0.008906;%  (0.005064, 0.01275)
% fitCp4;  Adjusted R-square: 0.9877
Tp4 =       22.16;%  (21.57, 22.75)
c04 =     0.02448;%  (0.02129, 0.02768)

%% Phonons fit coefficients after optimizing 
% Fit formula: c0 + (T/Tp)^3
% Fit includes weights
% fitCp1;  Adjusted R-square: 0.9906
data(1).Tp =       26.72;%  (26.2, 27.25)
data(1).c0 =     0.00714;%  (0.004807, 0.009472)
% fitCp2: no significant improvement
% fitCp3;  Adjusted R-square: 0.9979
data(3).Tp =       27.79;%  (26.73, 28.85)
data(3).c0 =   -0.000139;%  (-0.005698, 0.00542)
% fitCp4;  Adjusted R-square: 0.999
data(4).Tp =       21.93;%  (21.65, 22.21)
data(4).c0 =      0.0228;%  (0.01953, 0.02606)

%% Subtract phonons contribution from total Cp
rgn = [1,3,4];%length(data);
for i=rgn
    data(i).tbl.Cp_el = data(i).tbl.Cp - R*(data(i).tbl.T ./ data(i).Tp).^3 - R*data(i).c0;
end

%% Define variables for CFT (Curve Fitting Tool)
for i=rgn
% Important note: the following syntax is #1 in the list of things that one
% should NOT do in Matlab, HOWEVER I do not see any other iterative way of
% creating variables that can be accessed in the CFT.
eval(['fitCpel' num2str(i) '= data('  num2str(i) ').tbl.Cp_el /R']);
end

%% Fit parameters of mean-field Ising contribution to Cp
% fitCp1;  Adjusted R-square: 0.9733
data(1).Tc =       2.194;%  (2.186, 2.203)
data(1).sigma =    5e-4;% fixed by hand
% fitCp2;  Adjusted R-square: 0.9805
data(2).Tc = 2.171;%  (2.158, 2.185)
data(2).sigma =    0.002713;%  (0.001533, 0.003894)
% fitCp3;  Adjusted R-square: 0.9911
data(3).Tc =  2.189;%  (2.182, 2.195)
data(3).sigma =    1e-3;% fixed by hand
% fitCp4;  Adjusted R-square: 0.988
data(4).Tc =        2.174;%  (2.164, 2.183)
data(4).sigma =    0.001123;%  (0.0005623, 0.001684)

%% Subtract mean-field Ising + phonons contributions to Cp
for i=rgn
    data(i).tbl.Cpres = data(i).tbl.Cp... 
        - R*Cp_LFIM_phonons(data(i).tbl.T/data(i).Tc,data(i).sigma,data(i).Tp/data(i).Tc,data(i).c0);
end

%% Define variables for CFT (Curve Fitting Tool)
for i=rgn
% Important note: the following syntax is #1 in the list of things that one
% should NOT do in Matlab, HOWEVER I do not see any other iterative way of
% creating variables that can be accessed in the CFT.
eval(['fitCpres' num2str(i) '= data('  num2str(i) ').tbl.Cpres /R']);
end













