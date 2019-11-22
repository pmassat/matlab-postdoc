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








