load('betamodel.mat')
SR830MCE(8)=ImportSR830MCE('Run8_USINGSRSSeq_.5uA.dat',500e-9);
%%
figure
plot(SR830MCE(8).H,1./BetaModel(SR830MCE(8).H,SR830MCE(8).ResPlat))
hold on
plot(SR830MCE(3).H,1./BetaModel(SR830MCE(3).H,SR830MCE(3).ResPlat))
%%
newH8=gConvolve(SR830MCE(8).H,7);
newT8=gConvolve(1./BetaModel(SR830MCE(8).H,SR830MCE(8).ResPlat),7);
newD8=diff(newT8);
%%
figure
plot(newH8(30:(end-30)),newT8(30:(end-30)),'-')
ylim([0 3])
hold on
plot(newH(30:(end-30)),newT(30:(end-30)),'-')
ylim([0 3])

%%
figure;
plot(diff(SR830MCE(8).H));
hold on;
plot(diff(newH8));

%%
FilterFastUp8=diff(newH8)>0.4e-3;
FilterSlowUp8=diff(newH8)> 0.01e-3 &  diff(newH8)< 0.4e-3 ;
FilterFastDown8=diff(newH8)<-0.4e-3;
FilterSlowDown8=diff(newH8)< -0.01e-3 &  diff(newH8)> -0.4e-3 ;
%%
figure
SelectOnlyFastUp8.newH=newH8(FilterFastUp8);
SelectOnlyFastUp8.newH=SelectOnlyFastUp8.newH(30:(end-30));
SelectOnlyFastUp8.NewD=NewD8(FilterFastUp8);
SelectOnlyFastUp8.NewD=SelectOnlyFastUp8.NewD(30:(end-30));
SelectOnlyFastUp8.TPuck=SR830MCE(3).TPuck(FilterFastUp8);
SelectOnlyFastUp8.TPuck=SelectOnlyFastUp8.TPuck(30:(end-30));
plot(SelectOnlyFastUp8.newH,SelectOnlyFastUp8.NewD,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title('Traces Sweeping up fast in field ')

%%
figure
SelectOnlySlowUp8.newH=newH(FilterSlowUp8);
SelectOnlySlowUp8.newH=SelectOnlySlowUp8.newH(30:(end-30));
SelectOnlySlowUp8.NewD=NewD(FilterSlowUp8);
SelectOnlySlowUp8.NewD=SelectOnlySlowUp8.NewD(30:(end-30));
SelectOnlySlowUp8.TPuck=SR830MCE(3).TPuck(FilterSlowUp8);
SelectOnlySlowUp8.TPuck=SelectOnlySlowUp8.TPuck(30:(end-30));
plot(SelectOnlySlowUp8.newH,SelectOnlySlowUp8.NewD*100+SelectOnlySlowUp8.TPuck,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title('Traces Sweeping up slow in field ')
%%
figure
SelectOnlyFastDown8.newH=newH(FilterFastDown8);
SelectOnlyFastDown8.newH=SelectOnlyFastDown8.newH(30:(end-30));
SelectOnlyFastDown8.NewD=NewD(FilterFastDown8);
SelectOnlyFastDown8.NewD=SelectOnlyFastDown8.NewD(30:(end-30));
SelectOnlyFastDown8.TPuck=SR830MCE(3).TPuck(FilterFastDown8);
SelectOnlyFastDown8.TPuck=SelectOnlyFastDown8.TPuck(30:(end-30));
plot(SelectOnlyFastDown8.newH,SelectOnlyFastDown8.NewD*50+SelectOnlyFastDown8.TPuck,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title('Traces Sweeping Down fast in field ')

%%
figure
SelectOnlySlowDown8.newH=newH(FilterSlowDown8);
SelectOnlySlowDown8.newH=SelectOnlySlowDown8.newH(30:(end-30));
SelectOnlySlowDown8.NewD=NewD(FilterSlowDown8);
SelectOnlySlowDown8.NewD=SelectOnlySlowDown8.NewD(30:(end-30));
SelectOnlySlowDown8.TPuck=SR830MCE(3).TPuck(FilterSlowDown8);
SelectOnlySlowDown8.TPuck=SelectOnlySlowDown8.TPuck(30:(end-30));
plot(SelectOnlySlowDown8.newH,SelectOnlySlowDown8.NewD*100+SelectOnlySlowDown8.TPuck,'.')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title('Traces Sweeping Down slow in field ')

%%
figure
plot(SelectOnlyFastDown8.newH,SelectOnlyFastDown8.NewD*50+SelectOnlyFastDown8.TPuck,'.','DisplayName','Fast Down')
hold on
plot(SelectOnlySlowDown8.newH,SelectOnlySlowDown8.NewD*100+SelectOnlySlowDown8.TPuck,'.','DisplayName','Slow Down')
plot(SelectOnlySlowUp8.newH,SelectOnlySlowUp8.NewD*100+SelectOnlySlowUp8.TPuck,'.','DisplayName','Slow Up')

plot(SelectOnlyFastUp8.newH,SelectOnlyFastUp8.NewD*50+SelectOnlyFastUp8.TPuck,'.','DisplayName','Fast Up')
xlabel('Field(Oe)')
ylabel('He3 Temperature + dT_{Platform}/dH')
legend('show')