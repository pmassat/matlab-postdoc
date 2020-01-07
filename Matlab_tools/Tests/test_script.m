%% 
i = 15;
% Use these variables in Curve fitting tool
clear fitT fitCp fitCpErr fitwghts 
fitT = avgData(i).T;
fitCp = avgData(i).Cp/R;
fitCpErr = avgData(i).CpFullErr/R;
fitwghts = 1./fitCpErr;



