% fieldinfo0p94K = import manually data from "field_info.txt"
fieldinfo0p94K.FileName = strcat("HH0_",num2str(fieldinfo0p94K.FileID));
for i=1:length(fieldinfo0p94K.FileName)
    tbl2 = ImportNeutronsLineCut(strcat(fieldinfo0p94K.FileName(i,:),".txt"));
    nDatap94(i).file = fieldinfo0p94K.FileName(i,:);
    nDatap94(i).hh0 = tbl2.hh0;
    nDatap94(i).I = tbl2.I;
    nDatap94(i).dI = tbl2.dI;
    nDatap94(i).field = fieldinfo0p94K.H_T(i);
    nDatap94(i).temp = fieldinfo0p94K.T_K(i);
end
H = extractfield(nDatap94,'field');
%% Data formatting for curve fitting tool analysis
i=67;
hh67 = nDatap94(i).hh0;
I67 = nDatap94(i).I;
dI = nDatap94(i).dI;
%%
istart = length(H);
iend = 37;
%% Fit single Ikeda-Carpenter-pseudo-Voigt peak at high field
% Analysis parameters
hc = -8;% value of h in reciprocal space
% ufb = 0.99; % upper fit boundary = highest value of h-hc for which to include datapoints for the fit
field = H;%cell2mat( arrayfun(@(c) c.field, nDatap94(1:istart).', 'Uniform', 0) );
datExcld = nDatap94(istart).hh0<hc-.7 | nDatap94(istart).hh0>hc+0.55 |...
    (nDatap94(istart).hh0>hc-.35 & nDatap94(istart).hh0<hc-0.15)...
     | (nDatap94(istart).hh0>hc+0.15 & nDatap94(istart).hh0<hc+0.2);% Exclude 
% data points that correspond to other peaks as well as those that are too far away

%% Perform and plot fit
rng = istart:-1:iend;
I1 = 3e5; R1 = 0.1; a1 = 200; b1 = 0.1; g1 = 1e-3; s1 = 6.6e-3;% free parameters initial values
% freePrms1 = {'I',I1;'R',R1;'alpha',a1;'beta',b1;'gamma',g1;'sigma',s1;'x0',hc};%7 free parameters
% freePrms1 = {'I',I1;'R',R1;'beta',b1;'gamma',g1;'x0',hc};% 5 free parameters
freePrms1 = {'I',I1;'gamma',g1;'x0',hc};% 3 free parameters
% Note: the order in which free parameters are defined matters because of
% how the array of initial fitting parameters 'initParams' is defined
for i=rng
    myfit = fitICpV(nDatap94(i).hh0,nDatap94(i).I,hc); 
    myfit.dataExcl = datExcld;
    ap1 = myfit.allParams{1}; ap1('alpha') = 140; ap1('sigma') = 6.6e-3;% 'ap1' is shorter than 'myfit.allParams{1}'
    ap1('R')=0.0; ap1('beta')=0; ap1('I') = I1; ap1('gamma') = g1;
    myfit.freeParams = {freePrms1};
    Nprms = length(vertcat(myfit.freeParams{:}));% total number of free parameters
    fmt = "ENS pattern along [hh0] T=%.2fK H=%.3fT R=%.2f %i params fit";
    label = sprintf(fmt,nDatap94(i).temp,nDatap94(i).field,ap1('R'),Nprms);
    fitStr = ['fit' int2str(Nprms) 'rslt']; gofStr = ['gof' int2str(Nprms)];
    [nDatap94(i).(fitStr), nDatap94(i).(gofStr)] = myfit.compute_fit();
    if mod(i,67)==0% select data to plot
        myfit.plot_fit(nDatap94(i).(fitStr)); title(label);
%         xlim([hc-.8 hc+.6]);
    end
    disp(label); disp(nDatap94(i).(fitStr)); disp(nDatap94(i).(gofStr));
end
%%
np = Nprms;
fitStrnp = ['fit' int2str(np) 'rslt']; gofStrnp = ['gof' int2str(np)];
if ~isfield(tbl,(fitStrnp)); tbl(1).(fitStrnp) = []; end% if structure tbl does not contain any field called (fitStrnp), create that field
flag = 0;% to check whether it is necessary to add a new row to the field or not
for i = 1:numel(tbl)
  if isempty(tbl(i).(fitStrnp)); Ntbl = i; flag = 1; break; end%
  %if a row of tbl.(fitStrnp) is empty, use it to store the following table
end
if ~flag; Ntbl = numel(tbl)+1; end% if all the rows of tbl.(fitStrnp) are non-empty, store table in a new row
tbl(Ntbl).(fitStrnp) = table(field(rng)','VariableNames',{'Field_Oe'});% first column of the table contains magnetic field values
cfn = coeffnames(nDatap94(istart).(fitStrnp));% extract fit coefficient names, i.e. free parameters 
for nc=1:np% for each free parameter
    prm = extract_structure_field(nDatap94(rng),fitStrnp,cfn{nc});% get its value
    cft = cell2mat(arrayfun(@(c) confint(c.(fitStrnp)),nDatap94(rng).','Uniform',0));% get the 95% confidence intervals
    cftm = cft(1:2:end,nc); cftp = cft(2:2:end,nc);% get the lower and upper bounds of that interval, respectively
    prmErrm = prm-cftm;prmErrp = prm-cftp;% calculate corresponding negative and positive absolute errors
    relErrm = abs(prmErrm./prm);relErrp = abs(prmErrp./prm);% and relative errors
    tbl(Ntbl).(fitStrnp).(cfn{nc}) = prm;% store fit parameter value in a new column
    tbl(Ntbl).(fitStrnp).([cfn{nc} '_RelErr']) = relErrm;% and relative error in another new column
end
if ~any(strcmp('R',tbl(1).(fitStrnp).Properties.VariableNames))% if there is no field called 'R' in table tbl.(fitStrnp)
    tbl(Ntbl).(fitStrnp).R = ap1('R')*ones(length(rng),1);% store value of parameter R in another column (this is unnecessary if R is a free parameter)
end
rsquarenp = extract_structure_field(nDatap94(rng),gofStrnp,'rsquare');% extract r^2 value from goodness of fit
tbl(Ntbl).(fitStrnp).Rsquare = rsquarenp;% store r^2 value in a new column
flag = 0;% reset flag for next run

%% Write table to file
fileChar = [fitStr '_R=' sprintf('%2.e',ap1('R')) '.txt'];
fileID = fopen(fileChar,'a');
% fprintf(fileID,'\nValues of free parameters after fit:\n');
% fprintf(fileID,'%s\n',nDatap94(np).tbl);
writetable(tbl(Ntbl).(fitStrnp),fileChar);
fprintf(fileID,'\nInitial values of free parameters:\n');
S = string(myfit.freeParams{1});
K = keys(myfit.allParams{1}); V = values(myfit.allParams{1});
fprintf(fileID,'%s\n',strcat(S(:,1)," = ",S(:,2)));
fprintf(fileID,'\nFixed values of all parameters (disregard free parameters):\n');
fprintf(fileID,'%s\n',strcat(K," = ",string(V)));
fclose(fileID);
%% Extract 7-fit parameters from data structure and create table containing fit parameters 
[gamma,relErrGm,alpha,relErrAm,sigma,relErrSm] = ENS_peak_fit_extract_params(nDatap94(rng),'fit7rslt');
% [gamma7,relErrGm7,alpha7,relErrAm7,sigma7,relErrSm7] = ENS_peak_fit_extract_params(nDatap94(rng),'fit7rslt',1);
% Add third argument (e.g. value 1) to check equality of positive and negative error bars for all three parameters
% Note: the case where they are not equal has NOT been tested, so it may
% raise an error when that happens (although it is unlikely)
rsquare = extract_structure_field(nDatap94(rng),'gof7','rsquare');
fitPrms = {'gamma','relErrG','alpha','relErrA','sigma','relErrS','rsquare'};
vars = cell(1, length(fitPrms));
for iv=1:length(vars); vars{iv} = [fitPrms{iv} int2str(Nprms)]; end
% Check that both the values and errors make sense
nDatap94(Nprms).tbl = table(field(iend:istart),gamma,relErrGm,alpha,relErrAm,sigma,relErrSm,rsquare,...
    'VariableNames',[{'Field_Oe'},vars]);
%% Compute average value of parameter to be fixed in next iterations
% calculate mean of alpha over j data points at highest fields,
% where the data is single-peaked and sample behavior does not change (too much)
% alpha5 only depends on the behavior of the neutron beam and should thus
% not change under applied magnetic field
Navg = 20;
am7 = ones(Navg,1); stdam7 = ones(Navg,1);
sm7 = ones(Navg,1); stdsm7 = ones(Navg,1);% same with sigma, 
% which only depends on the instrument resolution
for j=1:Navg
    am7(j) = mean(alpha(1:j)); sm7(j) = mean(sigma(1:j));
    stdam7(j) = std(alpha(1:j)); stdsm7(j) = std(sigma(1:j));
end

%% Perform fit with 5 free parameters: I, alpha, gamma, sigma, x0
for i=istart%:-1:iend
    label = strcat("T=",num2str(round(nDatap94(i).temp,2)),"K & H=",num2str(nDatap94(i).field),"T");
    myfit = fitICpV(nDatap94(i).hh0,nDatap94(i).I,hc); 
    myfit.dataExcl = datExcld;
    [nDatap94(i).fit5rslt, nDatap94(i).gof5] = myfit.compute_fit({'I',5e5;...
        'R',0.5;'beta',1;'gamma',1e-3;'sigma',1e-2;'x0',hc});
    if mod(i,10)==7% plot data every 10 fields
        myfit.plot_fit(nDatap94(i).fit7rslt);
        title(strcat("ENS pattern cut along [hh0] at ",label));
%         xlim([hc-.8 hc+.6]);
    end
    disp(strcat("Fit ICpV1 5 parameters at ",label)); disp(nDatap94(i).fit5rslt);
end
%% Extract columns from data structure
field = cell2mat( arrayfun(@(c) c.field, nDatap94(rng).', 'Uniform', 0) );
%% Extract 5-fit parameters from data structure
[gamma5,relErrGm5,alpha5,relErrAm5,sigma5,relErrSm5] = ENS_peak_fit_extract_params(nDatap94(rng),'fit5rslt');
% [gamma5,relErrGm5,alpha5,relErrAm5,sigma5,relErrSm5] = ENS_peak_fit_extract_params(nDatap94(rng),'fit5rslt',1);
% Add third argument (e.g. value 1) to check equality of positive and negative error bars for all three parameters
% Note: the case where they are not equal has not been tested, so it may
% raise an error when that happens (although it is unlikely)
rsquare5 = extract_structure_field(nDatap94(rng),'gof5','rsquare');
%% Create table containing fit parameters and corresponding relative standard errors
% Check that both the values and errors make sense
tbl5 = table(field,gamma5,relErrGm5,alpha5,relErrAm5,sigma5,relErrSm5,rsquare5,'VariableNames',...
    {'Field_Oe','gamma5','relErrG','alpha5','relErrA','sigma5','relErrS','rsquare5'})
%% Compute average value of parameter to be fixed in next iterations
% calculate mean of alpha5 over j data points at highest fields,
% where the data is single-peaked and sample behavior does not change (too much)
% alpha5 only depends on the behavior of the neutron beam and should thus
% not change under applied magnetic field
am5 = ones(20,1); stdam5 = ones(20,1);
sm5 = ones(20,1); stdsm5 = ones(20,1);% same with sigma5, 
% which only depends on the instrument resolution
for j=1:20
    am5(j) = mean(alpha5(1:j)); sm5(j) = mean(sigma5(1:j));
    stdam5(j) = std(alpha5(1:j)); stdsm5(j) = std(sigma5(1:j));
end 

%% Recompute fit with fixed sigma parameter, i.e. 4 free parameters
for i=rng
    [nDatap94(i).fit4rslt, nDatap94(i).gof4] = ENS_peak_fit_ICpV_function_4params(...
        nDatap94(i).hh0,nDatap94(i).I,hc,datExcld);
    label = strcat("T=",num2str(round(nDatap94(i).temp,2)),"K & H=",...
        num2str(nDatap94(i).field),"T");
    if mod(i,10)==7% plot data every 5 fields
        ENS_peak_fit_plot(nDatap94(i).hh0,nDatap94(i).I,hc,nDatap94(i).fit4rslt,datExcld)
        title(strcat("ENS pattern cut along [hh0] at ",label))
        xlim([hc-.75 hc+.6]);
    end
    disp(strcat("Fit ICpV1 4 parameters at ",label));
    disp(nDatap94(i).fit4rslt)
end
%% Extract 4-fit parameters from data structure
[gamma4,relErrGn4,alpha4,relErrAn4] = ENS_peak_fit_extract_params(nDatap94(rng),'fit4rslt');
% Add third argument (e.g. value 1) to check equality of positive and negative error bars for all three parameters
rsquare4 = extract_structure_field(nDatap94(rng),'gof4','rsquare');
%% Create table containing fit parameters and corresponding relative standard errors
% Check that both the values and errors make sense
tbl4 = table(field,gamma4,relErrGn4,alpha4,relErrAn4,rsquare4,'VariableNames',...
    {'Field_Oe','gamma4','relErrG','alpha4','relErrA','rsquare4'})
%% Recompute average value of parameter to be fixed in next iterations
am4 = ones(20,1); stdam4 = ones(20,1);
for j=1:20
    am4(j) = mean(alpha4(1:j));
    stdam4(j) = std(alpha4(1:j));
end

%% Recompute fit with fixed alpha parameter, i.e. 3 free parameters
for i=istart%:-1:iend
    [nDatap94(i).fit3rslt, nDatap94(i).gof3] = ENS_peak_fit_ICpV_function_3params(...
        nDatap94(i).hh0,nDatap94(i).I,hc,datExcld);
    label = strcat("T=",num2str(round(nDatap94(i).temp,2)),"K, H=",...
        num2str(nDatap94(i).field),"T, R=0, beta=0 & 3 params fit");
    if mod(i,10)==7% plot data every 10 fields
        ENS_peak_fit_plot(nDatap94(i).hh0,nDatap94(i).I,hc,nDatap94(i).fit3rslt,datExcld)
        title(strcat("ENS [hh0] cut at ",label))
        xlim([hc-.7 hc+.6]);
    end
    disp(strcat("Fit ICpV1 3 parameters at ",label));
    disp(nDatap94(i).fit3rslt)
end
%% Extract 3-fit parameters from data structure
[gamma3,relErrGn3] = ENS_peak_fit_extract_params(nDatap94(rng),'fit3rslt');
% Add third argument (e.g. value 1) to check equality of positive and negative error bars for all three parameters
rsquare3 = extract_structure_field(nDatap94(rng),'gof3','rsquare');
%% Create table containing fit parameters and corresponding relative standard errors
% Check that both the values and errors make sense
tbl3 = table(field,gamma3,relErrGn3,rsquare3,'VariableNames',...
    {'Field_Oe','gamma3','relErrG','rsquare3'})
%% Identify peak maximum and width
xM = ones(istart,1);
fwhm = ones(istart,1);
for i=rng
    I1 = nDatap94(i).fit3rslt.I;
    gamma = nDatap94(i).fit3rslt.gamma;
    x0 = nDatap94(i).fit3rslt.x0;
    fnfit = @(x)-I1*voigtIkedaCarpenter_ord(x,[0,140,0,gamma,6.6e-3,0.05,x0]);%fnfit = -1*[fit function] so that the maximum becomes a minimum
    xM(i) = fminbnd(fnfit,hc-.2,hc+.2);% Identify position of the maximum on interval [hc-0.2 hc+0.2]
    M = -fnfit(xM(i));% compute value of the maximum
    fd = @(x)abs(fnfit(x)+M/2);
    xhm1 = fminbnd(fd,xM(i)-0.2,xM(i));% identify x value for which function equals M/2 on interval [xM-0.2 xM]
    xhm2 = fminbnd(fd,xM(i),xM(i)+0.2);% same on interval [xM xM+0.2]
    fwhm(i) = xhm2 - xhm1;% FWHM of big peak
end
%%

%% 
% 2019-02-16
% 
% Next: 
% # Done. Compute average of alpha5 and sigma5 at high field. alpha5 is associated to 
% the time of flight of neutrons and sigma5 to the instrument resolution; neither 
% of these are expected to change under field
% # Done. Fit at lower fields down to H<~Hc using I, gamma5 and x0 as free
% fit parameters (fix alpha and sigma)
% # Identify the value of field H* under which the fits are not satisfactory anymore 
% # For fields <~ H*, fit using a sum of 2 ICpV
% # Extract physical parameters from fits: position of maximum (numerically) 
% and width (relation between alpha5, gamma5 and sigma5)
%% 2 peak fit
%% Compute 2 peak fit with 3 free parameters for each peak 
% and extract splitting between peaks
xM1 = ones(istart,1); xM2 = ones(istart,1); splitting = zeros(istart,1);
for i=1:istart
    label = strcat("T=",num2str(round(nDatap94(i).temp,2)),"K, H=",...
        num2str(nDatap94(i).field),"T");
    try
        if mod(i,10)==1% plot data every 10 values of field
            [nDatap94(i).fit6rslt, nDatap94(i).gof6] = ENS_peak_fit_ICpV2_function_6params(...
            nDatap94(i).hh0,nDatap94(i).I,hc,datExcld,1);% compute fit and plot
            title(strcat("ENS [hh0] cut at ",label," & 6 params fit"))
            xlim([hc-.75 hc+.6]);
        else
            [nDatap94(i).fit6rslt, nDatap94(i).gof6] = ENS_peak_fit_ICpV2_function_6params(...
            nDatap94(i).hh0,nDatap94(i).I,hc,datExcld);% compute fit but do not plot
        end
        disp(strcat(newline+"Fit ICpV1 6 parameters at ",label));
        disp(nDatap94(i).fit6rslt)
        I1 = nDatap94(i).fit6rslt.I1; I2 = nDatap94(i).fit6rslt.I2;
        gamma1 = nDatap94(i).fit6rslt.gamma1; gamma2 = nDatap94(i).fit6rslt.gamma2;
        x01 = nDatap94(i).fit6rslt.x01; x02 = nDatap94(i).fit6rslt.x02;
        fnfit1 = @(x)-(I1*voigtIkedaCarpenter_ord(x,[0,140,0,gamma1,6.6e-3,0.05,x01]));
        fnfit2 = @(x)-(I2*voigtIkedaCarpenter_ord(x,[0,140,0,gamma2,6.6e-3,0.05,x02]));
    %fnfit = -1*[fit function] so that the maximum becomes a minimum
        xM1(i) = fminbnd(fnfit1,hc-.2,hc+.2);% Identify position of the maximum on interval [hc-0.2 hc+0.2]
        xM2(i) = fminbnd(fnfit2,hc-.2,hc+.2);% Identify position of the maximum on interval [hc-0.2 hc+0.2]
        M1 = -fnfit2(xM1(i));
        M2 = -fnfit2(xM2(i));% compute value of the maximum, to check graphically if the result is correct
        splitting(i) = xM2(i) - xM1(i);
    %     if abs(nDatap94(i).fit6rslt.x02)-abs(nDatap94(i).fit6rslt.x01)>0
    % % if the absolute value of h for the smaller fitted peak is smaller than
    % % that of the bigger peak, the fit is not accurate anymore
    %         break 
    %     end
    catch % if the fitting process raises an error
        break
    end
    counter = i;% keep track of up to which field the fit was successful
end

%%
figure
plot(field,splitting,'.')





