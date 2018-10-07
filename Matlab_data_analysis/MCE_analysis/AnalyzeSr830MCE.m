load('betamodel.mat')
SR830MCE(3)=ImportSR830MCE('Run2_0p5uA.dat',300e-9);
%%
plot(SR830MCE(3).H,1./BetaModel(SR830MCE(3).H,SR830MCE(3).ResPlat))

%%
newH=gConvolve(SR830MCE(3).H,7);
newT=gConvolve(1./BetaModel(SR830MCE(3).H,SR830MCE(3).ResPlat),7);
NewD=diff(newT);
%%
figure
plot(newH(30:(end-30)),NewD(30:(end-29))*50+SR830MCE(3).TPuck(30:(end-30)),'.')
ylim([1 2])

%%
figure;
plot(diff(SR830MCE(3).H));
hold on;
plot(diff(newH));

%%
FilterFastUp=diff(newH)>1.5e-3;
FilterSlowUp=diff(newH)> 0.5e-3 &  diff(newH)< 1.5e-3 ;
FilterFastDown=diff(newH)<-1.5e-3;
FilterSlowDown=diff(newH)< -0.5e-3 &  diff(newH)> -1.5e-3 ;
%%
figure
SelectOnlyFastUp.newH=newH(FilterFastUp);
SelectOnlyFastUp.newH=SelectOnlyFastUp.newH(30:(end-30));
SelectOnlyFastUp.NewD=NewD(FilterFastUp);
SelectOnlyFastUp.NewD=SelectOnlyFastUp.NewD(30:(end-30));
SelectOnlyFastUp.TPuck=SR830MCE(3).TPuck(FilterFastUp);
SelectOnlyFastUp.TPuck=SelectOnlyFastUp.TPuck(30:(end-30));
plot(SelectOnlyFastUp.newH,SelectOnlyFastUp.NewD*50+SelectOnlyFastUp.TPuck,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title('Traces Sweeping up fast in field ')

%%
figure
SelectOnlySlowUp.newH=newH(FilterSlowUp);
SelectOnlySlowUp.newH=SelectOnlySlowUp.newH(30:(end-30));
SelectOnlySlowUp.NewD=NewD(FilterSlowUp);
SelectOnlySlowUp.NewD=SelectOnlySlowUp.NewD(30:(end-30));
SelectOnlySlowUp.TPuck=SR830MCE(3).TPuck(FilterSlowUp);
SelectOnlySlowUp.TPuck=SelectOnlySlowUp.TPuck(30:(end-30));
plot(SelectOnlySlowUp.newH,SelectOnlySlowUp.NewD*100+SelectOnlySlowUp.TPuck,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title('Traces Sweeping up slow in field ')
%%
figure
SelectOnlyFastDown.newH=newH(FilterFastDown);
SelectOnlyFastDown.newH=SelectOnlyFastDown.newH(30:(end-30));
SelectOnlyFastDown.NewD=NewD(FilterFastDown);
SelectOnlyFastDown.NewD=SelectOnlyFastDown.NewD(30:(end-30));
SelectOnlyFastDown.TPuck=SR830MCE(3).TPuck(FilterFastDown);
SelectOnlyFastDown.TPuck=SelectOnlyFastDown.TPuck(30:(end-30));
plot(SelectOnlyFastDown.newH,SelectOnlyFastDown.NewD*50+SelectOnlyFastDown.TPuck,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title('Traces Sweeping Down fast in field ')

%%
figure
SelectOnlySlowDown.newH=newH(FilterSlowDown);
SelectOnlySlowDown.newH=SelectOnlySlowDown.newH(30:(end-30));
SelectOnlySlowDown.NewD=NewD(FilterSlowDown);
SelectOnlySlowDown.NewD=SelectOnlySlowDown.NewD(30:(end-30));
SelectOnlySlowDown.TPuck=SR830MCE(3).TPuck(FilterSlowDown);
SelectOnlySlowDown.TPuck=SelectOnlySlowDown.TPuck(30:(end-30));
plot(SelectOnlySlowDown.newH,SelectOnlySlowDown.NewD*100+SelectOnlySlowDown.TPuck,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title('Traces Sweeping Down slow in field ')

%%
figure
plot(SelectOnlyFastDown.newH,SelectOnlyFastDown.NewD*50+SelectOnlyFastDown.TPuck,'.','DisplayName','Fast Down')
hold on
plot(SelectOnlySlowDown.newH,SelectOnlySlowDown.NewD*100+SelectOnlySlowDown.TPuck,'.','DisplayName','Slow Down')
plot(SelectOnlySlowUp.newH,SelectOnlySlowUp.NewD*100+SelectOnlySlowUp.TPuck,'.','DisplayName','Slow Up')

plot(SelectOnlyFastUp.newH,SelectOnlyFastUp.NewD*50+SelectOnlyFastUp.TPuck,'.','DisplayName','Fast Up')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
legend('show')