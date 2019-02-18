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
%%
i=67;
hh067 = nDatap94(i).hh0;
I67 = nDatap94(i).I;
dI = nDatap94(i).dI;
%% Fit single Ikeda-Carpenter-pseudo-Voigt peak at high field
%%
hc = -8;% value of h in reciprocal space
H = extractfield(nDatap94,'field');
istart = length(H);
iend = 37;
%%
for i=istart:-1:iend
% for i=length(fieldinfo0p94K.H_T(fieldinfo0p94K.H_T<0.75))
    [nDatap94(i).fit5rslt, nDatap94(i).gof5] = ENS_peak_fit_ICpV_function_5params(...
        nDatap94(i).hh0,nDatap94(i).I,hc);
    label = strcat("T=",num2str(round(nDatap94(i).temp,2)),"K & H=",...
        num2str(nDatap94(i).field),"T");
%     ENS_peak_fit_plot(nDatap94(i).hh0,nDatap94(i).I,hc,nDatap94(i).fit4rslt)
%     title(strcat("ENS pattern cut along [hh0] at ",label))
%     xlim([hc-.25 hc+.25]);
    disp(strcat("Fit ICpV1 5 parameters at ",label));
    disp(nDatap94(i).fit5rslt);
%     if nDatap94(i).fit5rslt.a1<0 | nDatap94(i).fit5rslt.a2<0
%        nDatap94(i).fit5rslt=[];
%        nDatap94(i).gof=[];
%        break
%     end
end
%% Extract columns from data structure
field = cell2mat( arrayfun(@(c) c.field, nDatap94(istart:-1:iend).', 'Uniform', 0) );
%% Extract 5-fit parameters from data structure
% [alpha5,gamma5,sigma5,rsquare5] = ENS_peak_fit_extract_params(nDatap94(istart),'fit5rslt')
[alpha5,gamma5,sigma5,relErrAm5,relErrGm5,relErrSm5] = ENS_peak_fit_extract_params(nDatap94(istart:-1:iend),'fit5rslt');
% Add third argument (e.g. value 1) to check equality of positive and negative error bars for all three parameters
rsquare5 = extract_structure_field(nDatap94(istart:-1:iend),'gof5','rsquare');
%% Create table containing fit parameters and corresponding relative standard errors
% Check that both the values and errors make sense
tbl5 = table(field,alpha5,relErrAm5,gamma5,relErrGm5,sigma5,relErrSm5,rsquare5,'VariableNames',...
    {'Field_Oe','alpha5','relErrA','gamma5','relErrG','sigma5','relErrS','rsquare5'})
%% Compute average value of parameter to be fixed in next iteration
am = mean(alpha5(1:10));% calculate mean of alpha5 over 10 data points at highest fields,
% where the data is single-peaked and sample behavior does not change (too much)
% alpha5 only depends on the behavior of the neutron beam and should thus
% not change under applied magnetic field
sm = mean(sigma5(1:10));% same with sigma5, which only depends on the instrument resolution
%% Recompute fit with fixed sigma5=sm parameter
for i=istart:-1:iend
    [nDatap94(i).fit4rslt, nDatap94(i).gof4] = ENS_peak_fit_ICpV_function_4params(...
        nDatap94(i).hh0,nDatap94(i).I,hc);
    label = strcat("T=",num2str(round(nDatap94(i).temp,2)),"K & H=",...
        num2str(nDatap94(i).field),"T");
    ENS_peak_fit_plot(nDatap94(i).hh0,nDatap94(i).I,hc,nDatap94(i).fit4rslt)
    title(strcat("ENS pattern cut along [hh0] at ",label))
    xlim([hc-.25 hc+.25]);
    disp(strcat("Fit ICpV1 4 parameters at ",label));
    disp(nDatap94(i).fit4rslt)
end
%% Extract 4-fit parameters from data structure
[gamma4,relErrGn4,alpha4,relErrAn4] = ENS_peak_fit_extract_params(nDatap94(istart:-1:iend),'fit4rslt');
% [alpha4,gamma4,relErrAm4,relErrGm4] = ENS_peak_fit_extract_params(nDatap94(istart:-1:iend),'fit4rslt');
% Add third argument (e.g. value 1) to check equality of positive and negative error bars for all three parameters
rsquare4 = extract_structure_field(nDatap94(istart:-1:iend),'gof4','rsquare');
%% Create table containing fit parameters and corresponding relative standard errors
% Check that both the values and errors make sense
tbl4 = table(field,alpha4,relErrAn4,gamma4,relErrGn4,rsquare4,'VariableNames',...
    {'Field_Oe','alpha4','relErrA','gamma4','relErrG','rsquare4'})
%% Identify peak maximum and width
i = 67;
hh0 = nDatap94(i).hh0;
Ip85 = nDatap94(i).fit5rslt.I
alphap85 = nDatap94(i).fit5rslt.alpha5
gammap85 = nDatap94(i).fit5rslt.gamma5
sigmap85 = nDatap94(i).fit5rslt.sigma5
x0p85 = nDatap94(i).fit5rslt.x0
fnfit = @(x)-Ip85*voigtIkedaCarpenter(x,[gammap85,sigmap85,alphap85,1,0,0.05,x0p85]);%fnfit = -1*[fit function] so that the maximum becomes a minimum
xM = fminbnd(fnfit,hc-.2,hc+.2)% Identify position of the maximum on interval [hc-0.2 hc+0.2]
M = -fnfit(xM)% compute value of the maximum
fd = @(x)abs(fnfit(x)+M/2);
xhm1 = fminbnd(fd,xM-0.2,xM)% identify x value for which function equals M/2 on interval [xM-0.2 xM]
xhm2 = fminbnd(fd,xM,xM+0.2)% same on interval [xM xM+0.2]
%%

%% 
% 2019-02-16
% 
% Next: 
% # compute average of alpha5 and sigma5 at high field. alpha5 is associated to 
% the time of flight of neutrons and sigma5 to the instrument resolution; neither 
% of these are expected to change under field
% # Fit at lower fields down to H<~Hc using I, gamma5 and x0 as free fit parameters
% # Identify the threshold of rsquare5 under which the fits are not satisfactory 
% anymore and the value of field H* under which rsquare5 is below this threshold
% # For fields <~ H*, fit using a sum of ICpV
% # Extract physical parameters from fits: position of maximum (numerically) 
% and width (relation between alpha5, gamma5 and sigma5)
%%
lsplt = 30;
for i=1:length(H)
    if isequal(nDatap94(i).fit5rslt,[])
        lsplt = i-1;
        break
    end
    nDatap94(i).a1 = nDatap94(i).fit5rslt.a1;
    nDatap94(i).b1 = nDatap94(i).fit5rslt.b1;
    nDatap94(i).c1 = nDatap94(i).fit5rslt.c1;
    nDatap94(i).a2 = nDatap94(i).fit5rslt.a2;
    nDatap94(i).b2 = nDatap94(i).fit5rslt.b2;
    nDatap94(i).c2 = nDatap94(i).fit5rslt.c2;
    levelCfd = 2*tcdf(-1,nDatap94(i).gof.dfe);% level of confidence of fit
    cft = confint(nDatap94(i).fit5rslt,1-levelCfd);% confidence interval of fit
    stdC1= (cft(2,2)-cft(1,2))/2;% standard error on parameter b1
    stdC2= (cft(2,5)-cft(1,5))/2;% standard error on parameter b2
    nDatap94(i).splitting = abs(nDatap94(i).fit5rslt.b2-nDatap94(i).fit5rslt.b1)/8;
    nDatap94(i).dsplt = (stdC2+stdC1)/8;
end
%%
ntDatap94 = struct2table(nDatap94);
%%
hh0p85 = nDatap94(67).hh0;
Ip85 = nDatap94(67).I;
%%
Msplit = cell2mat(ntDatap94.splitting);% convert splitting table into array
Mdsplit = cell2mat(ntDatap94.dsplt);% same with error bars
MH = H(1:lsplt);% field array with same length as splitting array
avgtemp = round(mean(ntDatap94.temp),2);% average temperature of the full dataset
%%
A = 5.7e-3; Hc = 0.88;% fit parameters determined from Curve Fitting Tool
fsplit = @(x) A*sqrt(1-(x/Hc)^2);% fitting function
figure
fplot(fsplit,'-r','LineWidth',2);% plot fit
hold on
indsep = 30;% index of last element used for fit
errorbar(MH(1:indsep),Msplit(1:indsep),Mdsplit(1:indsep),'.b','MarkerSize',12);% plot data
plot(MH(indsep+1:end),Msplit(indsep+1:end),'xk','MarkerSize',6);% plot exclude data points
xlim([0 0.9]); ylim([0 6e-3]);
title(strcat("Splitting of [880] ENS peak along [hh0] at ",num2str(avgtemp),"K"))
xlabel("H (T)");ylabel("Relative peak splitting");
%% Fit with single Gaussian above Hc
%%
for i=length(H)%:-1:40
    [xData, yData] = prepareCurveData(nDatap94(i).hh0,nDatap94(i).I);
    dataExcl = nDatap94(i).hh0<hc-.15 | nDatap94(i).hh0>hc+0.2;
    ft = fittype( 'gauss1' );
    excludedPoints = excludedata( xData, yData, 'Indices', dataExcl );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0 -Inf 0];
    opts.StartPoint = [6e6 hc 0.01];
    opts.Exclude = excludedPoints;
    % Fit model to data.
    [nDatap94(i).fit1, nDatap94(i).gof1] = fit( xData, yData, ft, opts );
    % Plot
    figure
    h = plot( fitresult, xData, yData, excludedPoints, '.k');
    legend( h, 'I vs. hh0', 'Excluded I vs. hh0', 'fit Gaussian', 'Location', 'NorthEast' );
    title(strcat("ENS pattern cut along [hh0] at T=",num2str(round(nDatap94(i).temp,2)),...
        "K & H=",num2str(nDatap94(i).field),"T"))
    xlabel("hh0");ylabel("I (a.u.)");
    xlim([hc-.25 hc+.25]);
    grid on
    disp(nDatap94(i).fit1)
end
%%
figure
plot(H(1:lsplt),extractfield(nDatap94,'a1'),'.')
xlabel("H (T)");ylabel("a1");
%%
figure
plot(H(1:lsplt),extractfield(nDatap94,'a2'),'.')
xlabel("H (T)");ylabel("a2");
%%
figure
plot(H(1:lsplt),extractfield(nDatap94,'b1'),'.')
xlabel("H (T)");ylabel("b1");
%%
figure
plot(H(1:lsplt),extractfield(nDatap94,'b2'),'.')
xlabel("H (T)");ylabel("b2");
%%
figure
plot(H(1:lsplt),extractfield(nDatap94,'c1'),'.')
xlabel("H (T)");ylabel("c1");
%%
figure
plot(H(1:lsplt),extractfield(nDatap94,'c2'),'.')
xlabel("H (T)");ylabel("c2");