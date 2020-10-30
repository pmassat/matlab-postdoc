load('betamodel.mat')
dataname = "Run4\_0p5uA.dat";
date = " 31 Jul 2017";
SR830MCE(3)=ImportSR830MCE('Run4_0p5uA.dat',5e-7);
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
fast = 14e-4;%Fast sweep rate of magnetic field (in Tesla/s)
slow = 6e-4;%Slow sweep rate of magnetic field (in Tesla/s)
ref_sweep = (fast+slow)/2%Reference sweep rate of magnetic field (in Tesla/s) used 
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
%pts_ctff = 50;%Number of points to cut off when plotting the data in order 
%to avoid the Gaussian aberrations at the edges of the plots
% SelectOnlyFastUp.newH=SelectOnlyFastUp.newH(pts_ctff:(end-pts_ctff));
% SelectOnlyFastUp.NewD=SelectOnlyFastUp.NewD(pts_ctff:(end-pts_ctff));
% SelectOnlyFastUp.TPuck=SelectOnlyFastUp.TPuck(pts_ctff:(end-pts_ctff));
plot(SelectOnlyFastUp.newH,SelectOnlyFastUp.NewD*50+SelectOnlyFastUp.TPuck,'.')
xlabel('Field (T)')
ylabel('He3 Temperature + dT_{Platform}/dH')
maintitle = dataname + date + newline;
title(maintitle + sprintf("Sweeping up in field at %i Oe/s",fast*10^4))

%%
figure
SelectOnlyFastDown.newH=newH(FilterFastDown);
SelectOnlyFastDown.NewD=NewD(FilterFastDown);
SelectOnlyFastDown.TPuck=SR830MCE(3).TPuck(FilterFastDown);
plot(SelectOnlyFastDown.newH,SelectOnlyFastDown.NewD*50+SelectOnlyFastDown.TPuck,'.')
xlabel('Field (T)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title(maintitle + sprintf("Sweeping down in field at %i Oe/s",fast*10^4))

%%
figure
plot(SelectOnlyFastDown.newH,SelectOnlyFastDown.NewD*50+SelectOnlyFastDown.TPuck,'.',...
'DisplayName',sprintf('%i Oe/s down',fast*10^4))
hold on
plot(SelectOnlyFastUp.newH,SelectOnlyFastUp.NewD*50+SelectOnlyFastUp.TPuck,'.',...
'DisplayName',sprintf('%i Oe/s up',fast*10^4))
xlabel('Field (T)')
ylabel('He3 Temperature + dT_{Platform}/dH')
title(maintitle + "Comparing up and down sweeps")
legend('show')
