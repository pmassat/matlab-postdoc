load('betamodel.mat')
dataname = "Run3\_0p5uA.dat";
date = " 31 Jul 2017";
SR830MCE(3)=ImportSR830MCE('Run3_0p5uA.dat',5e-7);
%%
plot(SR830MCE(3).H,1./BetaModel(SR830MCE(3).H,SR830MCE(3).ResPlat))% R depends on T and H
% R and H are known in the experiment; we want to know T
% BetaModel results from the fit of 1/T as a function of R and H using the
% calibration file 'thermometer_calibration_raw.txt'
xlabel('Field (T)')
ylabel('Temperature (K)')
title('T-H MCE traces ')
%%
sgm=10;%number of points for smoothing
newH=gConvolve(SR830MCE(3).H,sgm);%smooth the field data
newT=gConvolve(1./BetaModel(SR830MCE(3).H,SR830MCE(3).ResPlat),sgm);%smooth the temperature data
NewD=diff(newT);%compute the derivative of the smoothed temperature
%%
figure
plot(newH(30:(end-30)),NewD(30:(end-29))*50+SR830MCE(3).TPuck(30:(end-30)),'-')%plot dT/dH versus H
%ylim([0.3 1.5])
xlabel('Field (T)')
ylabel('dT/dH + T_{puck}')
title('MCE effect')
%%
figure;
yyaxis left;
plot(diff(SR830MCE(3).H),'DisplayName','Raw');
hold on;
plot(diff(newH),'-g','DisplayName','Smoothed');
xlabel('Measurement time t (s)')
ylabel('dH/dt (T/s)')
title('Sweep rate and puck temperature')
legend('show')

yyaxis right;
plot(SR830MCE(3).TPuck,'DisplayName','Puck temperature')
ylabel('He3 temperature (K)')
%%
ref_sweep = 0.6e-3;%Reference sweep rate of magnetic field (in Tesla/s) used 
%to discriminate between high and low sweep rates
FilterFastUp = diff(newH)>ref_sweep;
FilterSlowUp = diff(newH)>0 &  diff(newH)<ref_sweep;
FilterFastDown = diff(newH)<-ref_sweep;
FilterSlowDown = diff(newH)<0 &  diff(newH)>-ref_sweep;
%%
figure
SelectOnlyFastUp.newH=newH(FilterFastUp);
SelectOnlyFastUp.NewD=NewD(FilterFastUp);
SelectOnlyFastUp.TPuck=SR830MCE(3).TPuck(FilterFastUp);
% SelectOnlyFastUp.newH=SelectOnlyFastUp.newH(30:(end-30));
% SelectOnlyFastUp.NewD=SelectOnlyFastUp.NewD(30:(end-30));
% SelectOnlyFastUp.TPuck=SelectOnlyFastUp.TPuck(30:(end-30));
plot(SelectOnlyFastUp.newH,SelectOnlyFastUp.NewD*50+SelectOnlyFastUp.TPuck,'.')
xlabel('Field (T)')
ylabel('He3 Temperature + dT_{Platform}/dH')
maintitle = dataname + date + newline;
title(maintitle + "Sweeping up in field at 8 Oe/s")
% ylim([0.3 1.4])
% breakyaxis([0.5 0.9])

%%
figure
SelectOnlyFastDown.newH=newH(FilterFastDown);
SelectOnlyFastDown.NewD=NewD(FilterFastDown);
SelectOnlyFastDown.TPuck=SR830MCE(3).TPuck(FilterFastDown);
plot(SelectOnlyFastDown.newH,SelectOnlyFastDown.NewD*50+SelectOnlyFastDown.TPuck,'.')
xlabel('Field (T)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title(maintitle + "Sweeping down in field at 8 Oe/s")
% ylim([0.3 1.4])
% breakyaxis([0.5 0.9])

%%
figure
plot(SelectOnlyFastDown.newH,SelectOnlyFastDown.NewD*50+SelectOnlyFastDown.TPuck,'.','DisplayName','8 Oe/s down')
hold on
plot(SelectOnlyFastUp.newH,SelectOnlyFastUp.NewD*50+SelectOnlyFastUp.TPuck,'.','DisplayName','8 Oe/s up')
xlabel('Field (T)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title(maintitle + "Comparing up and down sweeps")
legend('show')
