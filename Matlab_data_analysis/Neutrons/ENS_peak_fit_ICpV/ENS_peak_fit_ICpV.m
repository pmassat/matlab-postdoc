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
hc = -8.0;% center of unsplit peak in reciprocal space
%% Data formatting for curve fitting tool analysis
i=1;
hh1 = nDatap94(i).hh0;
I1dat = nDatap94(i).I;
dI = nDatap94(i).dI;
%% Color plot of intensity in H-(hh0) 2D-map
% Prepare plot
[X,Y] = meshgrid(hh1,H);
Ifull = cell2mat( arrayfun(@(c) c.I', nDatap94(1:length(nDatap94)).', 'Uniform', 0) );% intensity data combined in one big matrix
%%
% Perform plot 
xmargin = 0.1;
Xsel = X>hc-xmargin & X<hc+xmargin;
figure
surf(X(:,Xsel(1,:)),Y(:,Xsel(1,:)),Ifull(:,Xsel(1,:)),'EdgeColor','None')
xlim([hc-xmargin hc+xmargin]); ylim([0 max(H)]); colormap jet;
cb = colorbar; cbl = cb.Label; cbl.Rotation = 0;% horizontal colorbar label
% cbl.Position = [-1 8e6 0]; cbl.Interpreter = 'latex';
xlabel("(h h 0)"); ylabel("H (T)");
sfmt = '$T=%.2f$ K';% string format for surf plot title
txtrow = 60; txtcol = 3;
xs = X(txtrow,Xsel(1,:)); ys = Y(txtrow,Xsel(1,:)); is = Ifull(txtrow,Xsel(1,:));
% text(xs(txtcol),ys(txtcol),is(txtcol),[sprintf(sfmt,nDatap94(i).temp)],'HorizontalAlignment','left','FontSize',16);
ann = annotation('textbox',[0.17 0.82 0.2 0.1],'interpreter','latex','String',{sprintf(sfmt,nDatap94(1).temp)},...
    'FontSize',14,'FontName','Arial','LineStyle','-','EdgeColor','r',...
    'LineWidth',2,'BackgroundColor',[1 1 1],'Color','k');% add annotation
ann.FitBoxToText='on';% fit annotation box to text
caxis('auto');% auto rescale of color scale

%%
istart = 1;
iend = 11;%length(H);
%% Fit single Ikeda-Carpenter-pseudo-Voigt peak at high field
% Analysis parameters
% ufb = 0.99; % upper fit boundary = highest value of h-hc for which to include datapoints for the fit
field = H;%cell2mat( arrayfun(@(c) c.field, nDatap94(1:istart).', 'Uniform', 0) );

%% Perform and plot fit
xc = [-8.03 -7.99];% position of peaks in reciprocal space
lx = length(xc);
rng = istart:1:iend;
I1 = 8e4; R1 = 0.1; a1 = 200; b1 = 0.1; g1 = 1e-3; s1 = 6.6e-3;% free parameters initial values
I2 = 1.8e5; R2 = 0.1; a2 = 200; b2 = 0.1; g2 = 1e-3; s2 = 6.6e-3;% free parameters initial values
% freePrms1 = {'I',I1;'R',R1;'alpha',a1;'beta',b1;'gamma',g1;'sigma',s1;'x0',hc};%7 free parameters
% freePrms1 = {'I1',I1;'R1',R1;'beta1',b1;'gamma1',g1;'x01',hc(1)};% 5 free parameters
% freePrms2 = {'I2',I2;'R2',R2;'beta2',b2;'gamma2',g2;'x02',hc(2)};% 5 free parameters
freePrms1 = {'I1',I1;'x01',xc(1)};% 3 free parameters
freePrms2 = {'I2',I2;'x02',xc(2)};% 3 free parameters
% Note: the order in which free parameters are defined matters because of
% how the array of initial fitting parameters 'initParams' is defined
fmt = 'ENS pattern along [hh0] T=%.2fK H=%.3fT %i params fit';% string format for plot title
for i=11%rng
    nData1 = nDatap94(i).hh0(nDatap94(i).dI>0);
    datExcld = nData1<min(hc)-.7 | nData1>max(hc)+0.55 |...
        (nData1>min(hc)-.35 & nData1<min(hc)-0.15) |...
        (nData1>max(hc)+0.14 & nData1<max(hc)+0.2);%
%     (nData1>max(hc)+0.045 & nData1<max(hc)+0.065);% Exclude 
% data points that correspond to other peaks as well as those that are too far away
% datExcld should be defined on the x interval where obj.dY>0, otherwise
% the number of excluded points will be higher than the number of data points
    myfit = fitICpV(nDatap94(i).hh0,nDatap94(i).I,nDatap94(i).dI,xc); 
    myfit.dataExcl = datExcld;
    ap1 = myfit.allParams{1}; ap1('alpha') = 140; ap1('sigma') = 6.6e-3;% 'ap1' is shorter than 'myfit.allParams{1}'
    ap1('R')=0.0; ap1('beta')=0; ap1('I') = I1; ap1('gamma') = 0;
    ap2 = myfit.allParams{2}; ap2('gamma') = 0;
    myfit.freeParams = {freePrms1,freePrms2};
    Nprms = length(vertcat(myfit.freeParams{:}));% total number of free parameters
    label = sprintf(fmt,nDatap94(i).temp,nDatap94(i).field,Nprms);
    fitStr = ['fit'  int2str(lx) 'ICpV' int2str(Nprms) 'rslt']; 
    gofStr = ['gof'  int2str(lx) 'ICpV' int2str(Nprms)];
    [nDatap94(i).(fitStr), nDatap94(i).(gofStr)] = myfit.compute_fit();
    if mod(i,10)==1% select data to plot
        myfit.plot_fit(nDatap94(i).(fitStr)); title(label);
%         xlim([hc-.8 hc+.6]);
    end
    disp(label); disp(nDatap94(i).(fitStr)); disp(nDatap94(i).(gofStr));
end
%%
np = Nprms;
fitStrnp = fitStr; gofStrnp = gofStr;
if ~exist('tbl','var') || ~isfield(tbl,(fitStrnp)); tbl(1).(fitStrnp) = []; end
% if structure tbl does not contain any field called (fitStrnp), create that field
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
if ~any(strcmp('R',tbl(Ntbl).(fitStrnp).Properties.VariableNames))% if there is no field called 'R' in table tbl.(fitStrnp)
    tbl(Ntbl).(fitStrnp).R = ap1('R')*ones(length(rng),1);% store value of parameter R in another column (this is unnecessary if R is a free parameter)
end
rsquarenp = extract_structure_field(nDatap94(rng),gofStrnp,'rsquare');% extract r^2 value from goodness of fit
tbl(Ntbl).(fitStrnp).Rsquare = rsquarenp;% store r^2 value in a new column
flag = 0;% reset flag for next run

%% Write table to file
fileChar = [fitStr '.txt'];% '_R=' sprintf('%2.e',ap1('R'))
fileID = fopen(fileChar,'a');
% fprintf(fileID,'\nValues of free parameters after fit:\n');
% fprintf(fileID,'%s\n',nDatap94(np).tbl);
writetable(tbl(Ntbl).(fitStrnp),fileChar);
fprintf(fileID,'\nInitial values of free parameters:\n');
S = string(vertcat(myfit.freeParams{:}));
fprintf(fileID,'%s\n',strcat(S(:,1)," = ",S(:,2)));
fprintf(fileID,'\nFixed values of all parameters (disregard free parameters):');
for il=1:lx
    fprintf(fileID,['\nPeak #' int2str(il) '\n']);
    K = keys(myfit.allParams{il}); V = values(myfit.allParams{il});
    fprintf(fileID,'%s\n',strcat(K," = ",string(V)));
end
fclose(fileID);
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

%% Identify peak maximum and width
xM = ones(length(rng),1);
fwhm = ones(length(rng),1);
for i=rng
    I1 = nDatap94(i).fitStr.I1;
%     gamma = nDatap94(i).fitStr.gamma;
    x0 = nDatap94(i).fitStr.x01;
    ftot = @(x)-I1*voigtIkedaCarpenter_ord(x,[0,140,0,gamma,6.6e-3,0.05,x0]);%fnfit = -1*[fit function] so that the maximum becomes a minimum
    xM(i) = fminbnd(ftot,hc-.2,hc+.2);% Identify position of the maximum on interval [hc-0.2 hc+0.2]
    M = -ftot(xM(i));% compute value of the maximum
    fd = @(x)abs(ftot(x)+M/2);
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
xM1 = ones(length(rng),1); xM2 = ones(length(rng),1); splitting = zeros(length(rng),1);
% spltRE = ones(length(rng),1);% xM2re = ones(length(rng),1); xgap = ones(length(rng),1);
% X = linspace(hc-.1,hc+.1,501); d1X = diff(X); d2X = diff(d1X);
for i=rng
    label = sprintf(fmt,nDatap94(i).temp,nDatap94(i).field,Nprms);
    I1 = nDatap94(i).(fitStr).I1; I2 = nDatap94(i).(fitStr).I2;
%         gamma1 = nDatap94(i).fitStr.gamma1; gamma2 = nDatap94(i).fitStr.gamma2;
    x01 = nDatap94(i).(fitStr).x01; x02 = nDatap94(i).(fitStr).x02;
    f1 = @(x)-(I1*voigtIkedaCarpenter_ord(x,[0,140,0,0,0.05,6.6e-3,x01]));
    f2 = @(x)-(I2*voigtIkedaCarpenter_ord(x,[0,140,0,0,0.05,6.6e-3,x02]));
%fnfit = -1*[fit function] so that the maximum becomes a minimum
    xM1(i) = fminbnd(f1,xc(1)-.1,xc(1)+.1);% Identify position of the maximum on interval [hc-0.2 hc+0.2]
    xM2(i) = fminbnd(f2,xc(2)-.1,xc(2)+.1);% Identify position of the maximum on interval [hc-0.2 hc+0.2]
    splitting(i) = -(xM2(i) - xM1(i))/hc;
    
%     ftot = @(x) -f1(x)-f2(x);%full fitting function
%     M1 = ftot(xM1(i));%
%     M2 = ftot(xM2(i));%
%     [xgap(i),fgap] = fminbnd(ftot,xM1(i),xM2(i));%position of gap between peaks, if any
%     if fgap<ftot(xM1(i))% if there is actually a gap between both peaks, 
% %   the value of ftot at its position should be lower than at xM1
%         fd1 = @(x)abs(ftot(x)-0.95*M1);% fd is 0 when ftot = 0.95*fgap, positive elsewhere
%         fd2 = @(x)abs(ftot(x)-0.95*M2);% fd is 0 when ftot = 0.95*fgap, positive elsewhere
%         xM1re(i) = fminbnd(fd1,xM1(i),xgap(i))-xM1(i);% lower value for which fd goes to 0
% %         xMerr2(i) = xM2(i)-fminbnd(fd2,xgap(i),xM2(i));% higher value ---
% % Why does the above expression yield a smaller value than this one: ??
%         xM2re(i) = xM2(i)-(fminbnd(fd2,xM2(i)-.1,xM2(i))+fminbnd(fd2,xM2(i),xM2(i)+.1))/2;% error bar on second peak
%     else
%         f = ftot(X);% calculate values of fit function on interval X
%         d1f = diff(f)./d1X;% first derivative of f
%         d2f = diff(d1f)./(d1X(2:end)-d2X/2);% second derivative of f
%         pks = findpeaks(d2f);% find peaks in second derivative of f
% % if ftot consists of only one peak, its second derivative will have exactly two peaks
%         if length(pks)>2% if its second derivative has more than 2 peaks, 
% % it means that ftot has an inflexion point, i.e. a shoulder, which is related to the
% % existence of a second peak in the data; in this case, use wider error bars
%             fd1 = @(x)abs(ftot(x)-0.9*M1);% fd is 0 when ftot = 0.95*fgap, positive elsewhere
%             fd2 = @(x)abs(ftot(x)-0.9*M2);% fd is 0 when ftot = 0.95*fgap, positive elsewhere
%             xM1re(i) = (fminbnd(fd1,xM1(i)-.1,xM1(i))+fminbnd(fd1,xM1(i),xM1(i)+.1))/2-xM1(i);
% % error bar on first peak equals (x1(fd1=0)+x2(fd1=0))/2-xM1(i), where x1 and x2 are the
% % lower and upper x values where fd1=0, respectively
%             xM2re(i) = xM2(i)-(fminbnd(fd2,xM2(i)-.1,xM2(i))+fminbnd(fd2,xM2(i),xM2(i)+.1))/2;% error bar on second peak
%         else% if there is no clear second peak in the data
%             break
%         end
%     end
end
% estimation of error bars
spltRE = tbl.(fitStr).x01_RelErr + tbl.(fitStr).x02_RelErr;

%%
for i=51
    f = -(f1(X)+f2(X));
    d1f = diff(f)./d1X;
    d2f = diff(d1f)./(d1X(2:end)-d2X/2);
    d3f = diff(d2f)./d1X(2:end-1);
    if mod(i,1)==0% select data to plot
        figure
%         plot(X,f1);
%         plot(X(2:end)-dX/2,d1f);
        plot(X(2:end-1),d2f); hold on
%         plot(X(3:end-1)-d1X(2:end-1)/2,d3f);
        title(label)
%         xlim([hc-.8 hc+.6]);
    end
end
%%
figure
% plot(field,splitting,'.')
errorbar(field,splitting,spltRE(rng),'.','MarkerSize',12)





