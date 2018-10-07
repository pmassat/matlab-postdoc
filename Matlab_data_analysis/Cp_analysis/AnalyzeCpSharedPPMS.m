%% This routine is intended at analyzing Cp data acquired with the shared PPMS in McCullough 015
%% To do: add error bars
% For T data, use the percent temperature increase
% For Cp values, use the SampHCErrmJmoleK field
% For averaged data, compute reduced errors: dxm = 1/n*sqrt(sum(dx^2)) (to be checked)
%% Import data
Data = importCpSharedPPMS_('20180322_TmVO4-LS5228-MP3-Plt-HC1803_Cp.dat');
%% Assign variable names and get rid of NaN rows
% H=[FieldOersted];
T=[Data.SampleTempKelvin];
Cp=[Data.SampHCmJmoleK];
dT = [Data.TempRiseKelvin];

whichPoints = isfinite(T) & isfinite(Cp) & isfinite(dT);
% H=H(whichPoints);
T=T(whichPoints);
Cp=Cp(whichPoints);
dT = dT(whichPoints);

%% Plot the full dataset
figure
% hmax=5500;
scatter(T,Cp,'.')
% xlim([0 hmax])
% ylim([0 3])
% xlabel('Field (Oe)')
xlabel('Temperature (K)')
ylabel('Cp (mJ/K/mol)')
title('Heat capacity of TmVO4-LS5228-MP3-Plt-HC1803')

%% Separate data according to the percent temperature increase of the measurement
clear separatedCpData
% temprise = zeros(length(T));
temprise = dT./T;
septemprise = [0.001,0.0035,0.01];
strtol = [0.001,0.0015,0.003];

for i = 1:length(septemprise)
    wp = abs(temprise-septemprise(i))< strtol(i);
%     separatedCpData(i).H = H(wp);
    separatedCpData(i).T = T(wp);
    separatedCpData(i).Cp = Cp(wp);
    separatedCpData(i).temprise = temprise(wp);
    
    [separatedCpData(i).T,wo] = sort(separatedCpData(i).T);
%     separatedCpData(i).H = separatedCpData(i).H(wo);
    separatedCpData(i).Cp = separatedCpData(i).Cp(wo);
    separatedCpData(i).temprise = separatedCpData(i).temprise(wo);
end

%%
figure
for i=1:length(septemprise)
    plot(separatedCpData(i).T,separatedCpData(i).Cp,'-+')
    hold on
end
xlabel('Temperature (K)')
ylabel('Cp (mJ/K/mol)')
legendCell = cellstr(num2str(100*septemprise', '%.2f %%'));
legend(legendCell)
title('Heat capacity of TmVO4-LS5228-MP3-Plt-HC1803')
hold off

%% Compute average of data points taken within 4.5mK of one another
% 4.5mK is the empirical maximal dispersion of data points taken
% consecutively at a given temperature
[Ts,so] = sort(T);
Cps = Cp(so);
Tsm = [];
Cpsm = [];
for i = 1:3:length(Ts)
ind = [i]
    for j = i+1:i+4
        if abs(Ts(j)-Ts(i))<0.0045
            ind = [ind, j]
        end
        if j>= length(Ts)
            break
        end
    end
Tsm = [Tsm,mean(Ts(ind))]
Cpsm = [Cpsm,mean(Cps(ind))]
end

%% Plot averaged data
figure
plot(Tsm,Cpsm,'.','markers',12)
xlabel('Temperature (K)')
ylabel('Cp (mJ/K/mol)')
title('Heat capacity of TmVO4-LS5228-MP3-Plt-HC1803')

%%
%% Fitting H-T phase diagram -- 1
% tstep=0.05;
% hstep=500;
% [Hg,Tg] = meshgrid(0:hstep:hmax,0.365:tstep:3);% for griddata
% Hgl=0:hstep:hmax;% for gridfit
% Tgl=0.365:tstep:3;% for gridfit
% figure
% Cpg = gridfit(H,T,Cp,Hgl,Tgl);
% %Cpg = griddata(H,T,Cp,Hg,Tg);
% surf(Hg,Tg,Cpg)
% hold on
% scatter3(H,T,Cp,100,'.','.k')
% xlabel('Field (Oe)')
% ylabel('Temperature (K)')
% zlabel('Cp (uJ/K)')

% %%
% 
% steps = -50:50;
% x= tstep*steps;
% s=0.05;
% d1Gaussian = -exp(-x.^2/(2*s^2)).*x./sqrt(s^6*2*pi);
% d2Gaussian = exp(-x.^2/(2*s^2)).*(x.^2 - s^2)/sqrt(s^10*2*pi);
% 
% figure
% d1Cpg = conv2(Cpg,d1Gaussian','same')
% d2Cpg = conv2(Cpg,d2Gaussian','same')
% surf(Hg,Tg,-d1Cpg,'EdgeColor','none')
% xlabel('Field (Oe)')
% ylabel('Temperature (K)')
% zlabel('-dCp/dT (uJ/K^2)')
% 
% 
%%
%% Fitting H-T phase diagram -- 2
%% 
for i=1:length(septemprise)
    separatedCpData(i).f=fitSpline(separatedCpData(i).T,-separatedCpData(i).Cp);
%Is fitSpline a Matlab function? Cannot find it in the help.
%    separatedCpData(i).d2=differentiate(separatedCpData(i).f,separatedCpData(i).T);
    plot(separatedCpData(i).f,'deriv1')
    hold on
end
legend(legendCell)
xlabel('Temperature (K)')
ylabel('-dCp/dT (uJ/K^2)')
hold off

%%
figure
% n = length(fields);
% C=[]
% for i=0:n-1
%     C=[C; 1-i/n i/n 0];
% end
% set(0,'defaultaxescolororder',C) %red to green
for i=1:length(septemprise)
    separatedCpData(i).f=fitSpline(separatedCpData(i).T,separatedCpData(i).Cp);
    separatedCpData(i).d2=differentiate(separatedCpData(i).f,separatedCpData(i).T);
    scatter3(separatedCpData(i).H,separatedCpData(i).T,-separatedCpData(i).d2,'filled')
    hold on
end
xlim([0 hmax])
ylim([0 3])
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('-dCp/dT (uJ/K)')
hold off

%%
T2=[];
H2=[];
d2Cp2=[];
for i=1:length(septemprise)
    T2= [T2; separatedCpData(i).T];
    H2=[ H2; separatedCpData(i).H];
    d2Cp2=[d2Cp2; separatedCpData(i).d2];
end

tstep=0.025;
% hmax=8000;
% hstep=500;
% [Hg,Tg] = meshgrid(0:hstep:hmax,0.365:tstep:3);%for griddata
Hgl=0:hstep:hmax;%for gridfit
Tgl=0.365:tstep:3;%for gridfit
figure
d2Cpgf = gridfit(H2,T2,d2Cp2,Hgl,Tgl);
d2Cpg = griddata(H2,T2,d2Cp2,Hg,Tg);
surf(Hg,Tg,-d1Cpg,'EdgeColor','none')
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('-dCp/dT (uJ/K^3)')
hold on
%scatter3(H2,T2,-d2Cp2,'.k')
xlim([0 hmax])
ylim([0 3])
