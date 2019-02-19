% fieldinfo0p94K = import manually data from "field_info.txt"
fieldinfo0p94K.FileName = strcat("HH0_",num2str(fieldinfo0p94K.FileID));
for i=1:length(fieldinfo0p94K.FileName)
    tbl2 = ImportNeutronsLineCut(strcat(fieldinfo0p94K.FileName(i,:),".txt"));
%     i = i+length(fieldsOrtho);
    nDatap94(i).file = fieldinfo0p94K.FileName(i,:);
    nDatap94(i).hh0 = tbl2.hh0;
    nDatap94(i).I = tbl2.I;
    nDatap94(i).dI = tbl2.dI;
    nDatap94(i).field = fieldinfo0p94K.H_T(i);
    nDatap94(i).temp = fieldinfo0p94K.T_K(i);
end
%% Extract field column from data structure
field = cell2mat( arrayfun(@(c) c.field, nDatap94(1:istart).', 'Uniform', 0) );
%% Data formatting for curve fitting tool analysis
i=67;
hh67 = nDatap94(i).hh0;
I67 = nDatap94(i).I;
dI = nDatap94(i).dI;
%% Fit single Ikeda-Carpenter-pseudo-Voigt peak at high field
% Analysis parameters
hc = -8;% value of h in reciprocal space
% ufb = 0.99; % upper fit boundary = highest value of h-hc for which to include datapoints for the fit
H = extractfield(nDatap94,'field');
istart = length(H);
iend = 37;
datExcld = nDatap94(istart).hh0<hc-.7 | nDatap94(istart).hh0>hc+0.55 | (nDatap94(istart).hh0>hc-.35 & nDatap94(istart).hh0<hc-0.15)...
     | (nDatap94(istart).hh0>hc+0.15 & nDatap94(istart).hh0<hc+0.2);% Exclude 
% data points that correspond to other peaks as well as those that are too far away

%% Perform fit with all 7 free parameters
for i=istart:-1:iend
    [nDatap94(i).fit7rslt, nDatap94(i).gof7] = ENS_peak_fit_ICpV_function_7params(...
        nDatap94(i).hh0,nDatap94(i).I,hc,datExcld);%ufb);
    label = strcat("T=",num2str(round(nDatap94(i).temp,2)),"K & H=",...
        num2str(nDatap94(i).field),"T & 7 params fit");
    if mod(i,10)==7% plot data every 10 fields
        ENS_peak_fit_plot(nDatap94(i).hh0,nDatap94(i).I,hc,nDatap94(i).fit7rslt,datExcld);%,ufb)
        title(strcat("ENS pattern cut along [hh0] at ",label))
        xlim([hc-.75 hc+.6]);
    end
    disp(strcat("ICpV1: ",label));
    disp(nDatap94(i).fit7rslt);
end
%% Extract 7-fit parameters from data structure
[gamma7,relErrGm7,alpha7,relErrAm7,sigma7,relErrSm7] = ENS_peak_fit_extract_params(nDatap94(istart:-1:iend),'fit7rslt');
% [gamma7,relErrGm7,alpha7,relErrAm7,sigma7,relErrSm7] = ENS_peak_fit_extract_params(nDatap94(istart:-1:iend),'fit7rslt',1);
% Add third argument (e.g. value 1) to check equality of positive and negative error bars for all three parameters
% Note: the case where they are not equal has not been tested, so it may
% raise an error when that happens (although it is unlikely)
rsquare7 = extract_structure_field(nDatap94(istart:-1:iend),'gof7','rsquare');
%% Create table containing fit parameters and corresponding relative standard errors
% Check that both the values and errors make sense
tbl7 = table(field,gamma7,relErrGm7,alpha7,relErrAm7,sigma7,relErrSm7,rsquare7,'VariableNames',...
    {'Field_Oe','gamma7','relErrG','alpha7','relErrA','sigma7','relErrS','rsquare7'})
%% Compute average value of parameter to be fixed in next iterations
% calculate mean of alpha5 over j data points at highest fields,
% where the data is single-peaked and sample behavior does not change (too much)
% alpha5 only depends on the behavior of the neutron beam and should thus
% not change under applied magnetic field
am7 = ones(20,1); stdam7 = ones(20,1);
sm7 = ones(20,1); stdsm7 = ones(20,1);;% same with sigma5, 
% which only depends on the instrument resolution
for j=1:20
    am7(j) = mean(alpha7(1:j)); sm7(j) = mean(sigma7(1:j));
    stdam7(j) = std(alpha7(1:j)); stdsm7(j) = std(sigma7(1:j));
end

%% Perform fit with 5 free parameters: I, alpha, gamma, sigma, x0
for i=istart:-1:iend
    [nDatap94(i).fit5rslt, nDatap94(i).gof5] = ENS_peak_fit_ICpV_function_5params(...
        nDatap94(i).hh0,nDatap94(i).I,hc,datExcld);%ufb);
    label = strcat("T=",num2str(round(nDatap94(i).temp,2)),"K & H=",...
        num2str(nDatap94(i).field),"T");
%     ENS_peak_fit_plot(nDatap94(i).hh0,nDatap94(i).I,hc,nDatap94(i).fit5rslt,datExcld);%,ufb)
%     title(strcat("ENS pattern cut along [hh0] at ",label))
%     xlim([hc-.75 hc+.6]);
    disp(strcat("Fit ICpV1 5 parameters at ",label));
    disp(nDatap94(i).fit5rslt);
end
%% Extract columns from data structure
field = cell2mat( arrayfun(@(c) c.field, nDatap94(istart:-1:iend).', 'Uniform', 0) );
%% Extract 5-fit parameters from data structure
[gamma5,relErrGm5,alpha5,relErrAm5,sigma5,relErrSm5] = ENS_peak_fit_extract_params(nDatap94(istart:-1:iend),'fit5rslt');
% [gamma5,relErrGm5,alpha5,relErrAm5,sigma5,relErrSm5] = ENS_peak_fit_extract_params(nDatap94(istart:-1:iend),'fit5rslt',1);
% Add third argument (e.g. value 1) to check equality of positive and negative error bars for all three parameters
% Note: the case where they are not equal has not been tested, so it may
% raise an error when that happens (although it is unlikely)
rsquare5 = extract_structure_field(nDatap94(istart:-1:iend),'gof5','rsquare');
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
sm5 = ones(20,1); stdsm5 = ones(20,1);;% same with sigma5, 
% which only depends on the instrument resolution
for j=1:20
    am5(j) = mean(alpha5(1:j)); sm5(j) = mean(sigma5(1:j));
    stdam5(j) = std(alpha5(1:j)); stdsm5(j) = std(sigma5(1:j));
end 

%% Recompute fit with fixed sigma parameter, i.e. 4 free parameters
for i=istart:-1:iend
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
[gamma4,relErrGn4,alpha4,relErrAn4] = ENS_peak_fit_extract_params(nDatap94(istart:-1:iend),'fit4rslt');
% Add third argument (e.g. value 1) to check equality of positive and negative error bars for all three parameters
rsquare4 = extract_structure_field(nDatap94(istart:-1:iend),'gof4','rsquare');
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
[gamma3,relErrGn3] = ENS_peak_fit_extract_params(nDatap94(istart:-1:iend),'fit3rslt');
% Add third argument (e.g. value 1) to check equality of positive and negative error bars for all three parameters
rsquare3 = extract_structure_field(nDatap94(istart:-1:iend),'gof3','rsquare');
%% Create table containing fit parameters and corresponding relative standard errors
% Check that both the values and errors make sense
tbl3 = table(field,gamma3,relErrGn3,rsquare3,'VariableNames',...
    {'Field_Oe','gamma3','relErrG','rsquare3'})
%% Identify peak maximum and width
xM = ones(istart,1);
fwhm = ones(istart,1);
for i=istart:-1:iend
    I = nDatap94(i).fit3rslt.I;
    gamma = nDatap94(i).fit3rslt.gamma;
    x0 = nDatap94(i).fit3rslt.x0;
    fnfit = @(x)-I*voigtIkedaCarpenter_ord(x,[0,140,0,gamma,6.6e-3,0.05,x0]);%fnfit = -1*[fit function] so that the maximum becomes a minimum
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





