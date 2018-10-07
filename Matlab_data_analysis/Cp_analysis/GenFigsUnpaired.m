%%
figure
traces = [10:10:70 85];
traces = [2 10 30 50 70 85];
%traces=50;
for i = traces
    plot(HCDATA(i).T, HCDATA(i).DSUnpaired,'DisplayName',strcat('x=',num2str(HCDATA(i).xx)))
    hold on
end
legend('show')

xlabel('Temperature')
ylabel('dS per Mol Unpaired Mn (J/Mol K^2)')
xlim([0 2])
%%
figure

traces = [2 10 30 50 70 85];
for i = fliplr(traces)
    plot(HCDATA(i).T, HCDATA(i).HCUnpaired,'+','DisplayName',strcat('x=',num2str(HCDATA(i).xx)))
    hold on
    i
end
legend('show')
%plot(Cp2.Temp,(Cp2.Ba3Mn098V002O8-(Cp2.Ba3Mn2O8*0.98*0.98))./(2*0.98*0.02),'+')
xlabel('Temperature')
ylabel('Heat Capacity per Mol Unpaired Mn (J/Mol K)')
xlim([0 2])

%%
figure
%traces = [50];
for i = traces
    plot(HCDATA(i).T, HCDATA(i).DeltaEntropyUnpaired,'+')
    hold on
end
legend(cellstr(num2str(traces', 'x=0.%-d')))
plot([0 4],8.3145*log(3)*[1 1],'k--')
xlabel('Temperature')
%xlim([0 0.6])
ylabel('Entropy Change per Mol Unpaired Mn (J/Mol K)')
xlim([0 3])   

%%
figure
%traces = [50];
for i = traces
    plot(HCDATA(i).TExtrapolated, HCDATA(i).DeltaEntropyUnpairedExtrapolated,'+')
    hold on
end
legend(cellstr(num2str(traces', 'x=0.%-d')))
plot([0 4],8.3145*log(3)*[1 1],'k--')
xlabel('Temperature')
%xlim([0 0.6])
ylabel('Entropy Change per Mol Unpaired Mn (J/Mol K)')
xlim([0 3])   

%%
traces = [10:10:70 85];
traces = [10:10:30 50:10:70 85];
temps = 0.05:0.05:4;
T3=[];
X3=[];
dS3=[];
HC3=[];
DeltaEntropy3=[];
for i = traces
    T3 = [T3; HCDATA(i).T];
    i
    X3 = [X3; HCDATA(i).x];
    2
    dS3 = [dS3; HCDATA(i).DSUnpaired ];
    3
    HC3= [HC3; HCDATA(i).HCUnpaired ];
    4
    DeltaEntropy3=[DeltaEntropy3; HCDATA(i).DeltaEntropy];
    5
end
figure
scatter3(T3,X3,dS3)
h=gca;
xlabel(h,'Temperature(K)')
ylabel(h,'x(%V)')
zlabel(h,'dS (J/Mol K^2)')

figure
scatter3(T3,X3,HC3)
h=gca;
xlabel(h,'Temperature(K)')
ylabel(h,'x(%V)')
zlabel(h,'Heat Capacity (J/Mol K)')

figure
scatter3(T3,X3,DeltaEntropy3)
h=gca;
xlabel(h,'Temperature(K)')
ylabel(h,'x(%V)')
zlabel(h,'Entropy Change (J/Mol K)')

%%
[Tg,Xg] = meshgrid(0.05:0.01:4,0.1:0.01:0.85);
figure
dSg = griddata(T3,X3,dS3,Tg,Xg);
mesh(Tg,Xg,dSg);
hold on
scatter3(T3,X3,dS3,'r.')
h=gca;
xlabel(h,'Temperature(K)')
ylabel(h,'x(%V)')
zlabel(h,'dS (J/Mol K^2)')
%%
figure
imagesc([min(Tg(:)) max(Tg(:))],[min(Xg(:)) max(Xg(:))],dSg)
h=gca;
xlabel(h,'Temperature(K)')
ylabel(h,'x(%V)')
zlabel(h,'dS (J/Mol K^2)')
%%
figure
scatter3(T3,X3,dS3,'r.')
hold on
imagesc([min(Tg(:)) max(Tg(:))],[min(Xg(:)) max(Xg(:))],dSg)
h=gca;
h.YLim = [0 1];
h.XLim = [0.05 4];
xlabel(h,'Temperature(K)')
ylabel(h,'x(%V)')
zlabel(h,'dS (J/Mol K)')

%%
[Tg,Xg] = meshgrid(0.05:0.01:4,0.1:0.01:0.85);
figure
HCg = griddata(T3,X3,HC3,Tg,Xg);
mesh(Tg,Xg,HCg);
hold on
scatter3(T3,X3,HC3,'r.')
h=gca;
xlabel(h,'Temperature(K)')
ylabel(h,'x(%V)')
zlabel(h,'Heat Capacity (J/Mol K)')

figure
imagesc([min(Tg(:)) max(Tg(:))],[min(Xg(:)) max(Xg(:))],HCg)
h=gca;
xlabel(h,'Temperature(K)')
ylabel(h,'x(%V)')
zlabel(h,'dS (J/Mol K^2)')

figure
scatter3(T3,X3,HC3,'r.')
hold on
imagesc([min(Tg(:)) max(Tg(:))],[min(Xg(:)) max(Xg(:))],HCg)
h=gca;
h.YLim = [0 1];
h.XLim = [0.05 4];
xlabel(h,'Temperature(K)')
ylabel(h,'x(%V)')
zlabel(h,'dS (J/Mol K)')

clear Xg X3 Tg temps T3 i HCg HC3 dS3 dSg