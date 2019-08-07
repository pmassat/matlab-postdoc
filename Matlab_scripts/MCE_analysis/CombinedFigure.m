%%
figure
imagesc([0 0.8],[.365 3],-d2Cpg)
hold on
xlabel('Field (T)')
ylabel('Temperature (K)')
set(gca,'YDir','normal')
h=colorbar('eastoutside');
h.Label.String = '-d C_p /dT';
%plot(SelectOnlyFastDown.newH,SelectOnlyFastDown.NewD*90+SelectOnlyFastDown.TPuck,'.k','DisplayName','Fast Down')

%plot(SelectOnlySlowDown.newH,SelectOnlySlowDown.NewD*100+SelectOnlySlowDown.TPuck,'.','DisplayName','Slow Down')
%plot(SelectOnlySlowUp.newH,SelectOnlySlowUp.NewD*100+SelectOnlySlowUp.TPuck,'.','DisplayName','Slow Up')
plot([SelectOnlyFastUp6.newH(950:1460)],[SelectOnlyFastUp6.NewD(950:1460)*250+mean(SelectOnlyFastUp6.newT(950:1460)) ],'m.','DisplayName','250 *\partial \Delta T /\partial H')

plot(SelectOnlyFastUp.newH,SelectOnlyFastUp.NewD*120+SelectOnlyFastUp.TPuck,'.r','DisplayName','120*\partial \Delta T /\partial H')
legend('show')