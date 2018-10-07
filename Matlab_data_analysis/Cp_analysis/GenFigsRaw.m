%%
figure
traces = [10:10:70 85];
traces = [10:10:30 50:10:70 85];
for i = traces
    plot(HCDATA(i).T, HCDATA(i).DS)
    hold on
end
legend(cellstr(num2str(traces', 'x=0.%-d')))
xlabel('Temperature')
ylabel('dS per Mol  Mn (J/Mol K^2)')
title('Cp/T per Mol Ba3Mn2-2xV2xO8')
%%
figure
traces = [10:10:30 50:10:70 85];
estimates = [10 20 30 30 60 70 85];
Delta = 11;
subset =1:7;
for j = subset
    i = traces(j);
    plot(HCDATA(i).T, HCDATA(i).HC'-((100-estimates(j))^2/1e4)*getSchottky(Delta,HCDATA(i).T))
    hold on
end
legend([cellstr(num2str(traces(subset)', 'x=0.%-d'))])
xlabel('Temperature')
ylabel('Heat Capacity per Mol Ba3Mn2O8 (J/Mol K)')
title('Cp per Mol Ba3Mn2-2xV2xO8')  
%temps=0.01:0.01:4;
%plot(background.SampleTemp,background.SampleCp,'rd')
%plot(temps, getSchottky(16,temps),'^k')
%%
traces = [10:10:70 85];
traces = [10:10:30 50:10:70 85];
temps = 0.05:0.05:4;
T3=[];
X3=[];
dS3=[];
HC3=[];
for i = traces
    T3 = [T3; HCDATA(i).T'];
    X3 = [X3; HCDATA(i).x'];
    dS3 = [dS3; HCDATA(i).DS ];
    HC3= [HC3; HCDATA(i).HC ];
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

figure
imagesc([min(Tg(:)) max(Tg(:))],[min(Xg(:)) max(Xg(:))],dSg)
h=gca;
xlabel(h,'Temperature(K)')
ylabel(h,'x(%V)')
zlabel(h,'dS (J/Mol K^2)')

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
zlabel(h,'C/ (J/Mol K)')