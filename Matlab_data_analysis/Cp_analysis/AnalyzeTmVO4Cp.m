DATA(1)=ImportTmVO4Cp('TmVO4_Mosaic_2017-07-28.dat');
%DATA(2)=ImportTmVO4Cp('TmVO4_RF-E_2017-07-14.dat');
%%
% H=[DATA(1).FieldOersted; DATA(2).FieldOersted];
% T=[DATA(1).SampleTempKelvin; DATA(2).SampleTempKelvin];
% Cp=[DATA(1).SampHCJmoleK; DATA(2).SampHCJmoleK];
H=[DATA(1).FieldOersted];
T=[DATA(1).SampleTempKelvin];
Cp=[DATA(1).SampHCJmoleK];

whichPoints = isfinite(H) & isfinite(T) & isfinite(Cp);
H=H(whichPoints);
T=T(whichPoints);
Cp=Cp(whichPoints);

scatter3(H,T,Cp,'.')
xlim([0 6000])
ylim([0 3])
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('Cp (uJ/K)')

%%
tstep=0.05;
[Hg,Tg] = meshgrid(0:500:5000,0.365:tstep:3);% for griddata
Hgl=0:500:5000;% for gridfit
Tgl=0.365:tstep:3;% for gridfit
figure
Cpg = gridfit(H,T,Cp,Hgl,Tgl);
Cpg = griddata(H,T,Cp,Hg,Tg);
surfc(Hg,Tg,Cpg)
hold on
scatter3(H,T,Cp,100,'.','.k')
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('Cp (uJ/K)')

%%

steps = -50:50;
x= tstep*steps;
s=0.05;
d1Gaussian = -exp(-x.^2/(2*s^2)).*x./sqrt(s^6*2*pi);
d2Gaussian = exp(-x.^2/(2*s^2)).*(x.^2 - s^2)/sqrt(s^10*2*pi);

figure
d1Cpg = conv2(Cpg,d1Gaussian','same')
%d2Cpg = conv2(Cpg,d2Gaussian','same')
surf(Hg,Tg,d1Cpg,'EdgeColor','none')

%%
clear separatedCpData

fields = [10 1000 2000 3000 4000 4500 5000];

for i = 1:length(fields)
    wp = abs(H-fields(i))<50;
    separatedCpData(i).H = H(wp);
    separatedCpData(i).T = T(wp);
    separatedCpData(i).Cp = Cp(wp);
    
    [separatedCpData(i).T,wo] = sort(separatedCpData(i).T);
    separatedCpData(i).H = separatedCpData(i).H(wo);
    separatedCpData(i).Cp = separatedCpData(i).Cp(wo);
end

%%
for i=1:length(fields)
    plot(separatedCpData(i).T,separatedCpData(i).Cp,'-+')
    hold on
end
xlabel('Temperature (K)')
ylabel('Cp (uJ/K)')
legendCell = cellstr(num2str(fields', '%-d Oe'));
legend(legendCell)
hold off

%%

x=separatedCpData(4).T;
y=separatedCpData(4).Cp;

%%
for i=1:length(fields)
    separatedCpData(i).f=fitSpline(separatedCpData(i).T,-separatedCpData(i).Cp);
%Is fitSpline a Matlab function? Cannot find it in the help.
%    separatedCpData(i).d2=differentiate(separatedCpData(i).f,separatedCpData(i).T);
    plot(separatedCpData(i).f,'deriv1')
    hold on
end
legend(legendCell)
hold off

%%

for i=1:length(fields)
    separatedCpData(i).f=fitSpline(separatedCpData(i).T,separatedCpData(i).Cp);
    separatedCpData(i).d2=differentiate(separatedCpData(i).f,separatedCpData(i).T);
    scatter3(separatedCpData(i).H,separatedCpData(i).T,separatedCpData(i).d2,'.k')
    hold on
end
xlim([0 6000])
ylim([0 3])
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('dCp/dT (uJ/K)')


%%
T2=[];
H2=[];
d2Cp2=[];
for i=1:length(fields)
    T2= [T2; separatedCpData(i).T];
    H2=[ H2; separatedCpData(i).H];
    d2Cp2=[d2Cp2; separatedCpData(i).d2];
end

tstep=0.025;
[Hg,Tg] = meshgrid(0:500:5000,0.365:tstep:3);%for griddata
Hgl=0:500:5000;%for gridfit
Tgl=0.365:tstep:3;%for gridfit
figure
d2Cpgf = gridfit(H2,T2,d2Cp2,Hgl,Tgl);
d2Cpg = griddata(H2,T2,d2Cp2,Hg,Tg);
surf(Hg,Tg,-d2Cpg,'EdgeColor','none')
xlabel('Field (Oe)')
ylabel('Temperature (K)')
zlabel('-dCp/dT (uJ/K)')
hold on
scatter3(H2,T2,-d2Cp2,'.k')
xlim([0 5000])
ylim([0 3])
