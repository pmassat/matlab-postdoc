load('betamodel.mat')
dataname = "Run7\_0p5uA\_2017-08-03.dat";
date = " 03 Aug 2017";
SR830MCE(3)=ImportSR830MCE('Run7_0p5uA_2017-08-03.dat',5e-7);

%%
figure
plot(SR830MCE(3).H,1./BetaModel(SR830MCE(3).H, SR830MCE(3).ResPlat))% R depends on T and H
% R and H are known in the experiment; we want to know T
% BetaModel results from the fit of 1/T as a function of R and H using the
% calibration file 'thermometer_calibration_raw.txt'
xlabel('Field (T)')
ylabel('Temperature (K)')
title('T-H MCE traces ')

%%
sgm=10;%number of points for smoothing
newH=gConvolve(SR830MCE(3).H,sgm);%smooth the field data
newT=gConvolve(1./BetaModel(SR830MCE(3).H, SR830MCE(3).ResPlat),sgm);%smooth the temperature data
dT = diff(newT);%compute the derivative of the smoothed temperature
dH = diff(newH);

%%
figure
plot(newH(30:(end-30)),dT(30:(end-29))*50 + SR830MCE(3).TPuck(30:(end-30)),'-')%plot dT/dH versus H
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
fast = 28e-4;%Fast sweep rate of magnetic field (in Tesla/s)
slow = 14e-4;%Slow sweep rate of magnetic field (in Tesla/s)
ref_sweep = (fast+slow)/2;%Reference sweep rate of magnetic field (in Tesla/s) used 
%to discriminate between high and low sweep rates
SweepSel(2).sweeprate = fast;
SweepSel(4).sweeprate = slow;
SweepSel(1).sweeprate = -fast;
SweepSel(3).sweeprate = -slow;
SweepSel(2).Filter = diff(newH)>ref_sweep;
SweepSel(4).Filter = diff(newH)>1e-3 &  diff(newH)<ref_sweep;
SweepSel(1).Filter = diff(newH)<-ref_sweep;
SweepSel(3).Filter = diff(newH)<-1e-3 &  diff(newH)>-ref_sweep;

%% Create structure containing smoothed data separated by sweeprate
for i = 1:4
SweepSel(i).H = newH(SweepSel(i).Filter);
SweepSel(i).T = newT(SweepSel(i).Filter);
SweepSel(i).dT = dT(SweepSel(i).Filter);
SweepSel(i).dH = dH(SweepSel(i).Filter);
SweepSel(i).TPuck = SR830MCE(3).TPuck(SweepSel(i).Filter);
SweepSel(i).d2T = diff(SweepSel(i).dT);
SweepSel(i).dHmid = (SweepSel(i).dH(1:end-1) + SweepSel(i).dH(2:end))/2;
SweepSel(i).TPuckMid = (SweepSel(i).TPuck(1:end-1) + SweepSel(i).TPuck(2:end))/2;
SweepSel(i).Hmid = (SweepSel(i).H(1:end-1) + SweepSel(i).H(2:end))/2;
end
% 
% FastDown.newH = newH(FilterFastDown);
% FastDown.newT = newT(FilterFastDown);
% FastDown.dT = dT(FilterFastDown);
% FastDown.dH = dH(FilterFastDown);
% FastDown.TPuck = SR830MCE(3).TPuck(FilterFastDown);
% 
% SlowUp.newH = newH(FilterSlowUp);
% SlowUp.newT = newT(FilterSlowUp);
% SlowUp.dT = dT(FilterSlowUp);
% SlowUp.dH = dH(FilterSlowUp);
% SlowUp.TPuck = SR830MCE(3).TPuck(FilterSlowUp);
% 
% SlowDown.newH = newH(FilterSlowDown);
% SlowDown.newT = newT(FilterSlowDown);
% SlowDown.dT = dT(FilterSlowDown);
% SlowDown.dH = dH(FilterSlowDown);
% SlowDown.TPuck = SR830MCE(3).TPuck(FilterSlowDown);
% 
% %% Second order differentials
% FastUp.d2T = diff(FastUp.dT);
% FastDown.d2T = diff(FastDown.dT);
% SlowUp.d2T = diff(SlowUp.dT);
% SlowDown.d2T = diff(SlowDown.dT);
% FastDown.dHmid = (FastDown.dH(1:end-1)+FastDown.dH(2:end))/2;

%% MCE labels
mceystr = '$T_{\mathrm{bath}} + \Delta T_{\mathrm{MCE}}$';
dmceystr = '$T_{\mathrm{bath}} + dT_{\mathrm{Platform}}/dH$';
d2mceystr = '$T_{\mathrm{bath}} + d^{2}T_{\mathrm{Platform}}/dH^{2}$';

%%
figure
%pts_ctff = 50;%Number of points to cut off when plotting the data in order 
%to avoid the Gaussian aberrations at the edges of the plots
% FastUp.newH=FastUp.newH(pts_ctff:(end-pts_ctff));
% FastUp.dT=FastUp.dT(pts_ctff:(end-pts_ctff));
% FastUp.TPuck=FastUp.TPuck(pts_ctff:(end-pts_ctff));
plot(FastUp.newH, FastUp.dT*50 + FastUp.TPuck,'.')
xlabel('Field (T)')
ylabel(dmceystr)
maintitle = dataname + date + newline;
title(maintitle + sprintf("Sweeping up in field at %i Oe/s",fast*10^4))

%%
figure
plot(FastDown.newH, FastDown.dT*50 + FastDown.TPuck,'.')
xlabel('Field (T)')
ylabel(dmceystr)
title(maintitle + sprintf("Sweeping down in field at %i Oe/s",fast*10^4))

%%
figure
%pts_ctff = 50;%Number of points to cut off when plotting the data in order 
%to avoid the Gaussian aberrations at the edges of the plots
% SlowUp.newH=SlowUp.newH(pts_ctff:(end-pts_ctff));
% SlowUp.dT=SlowUp.dT(pts_ctff:(end-pts_ctff));
% SlowUp.TPuck=SlowUp.TPuck(pts_ctff:(end-pts_ctff));
plot(SlowUp.newH, SlowUp.dT*50 + SlowUp.TPuck,'.')
xlabel('Field (T)')
ylabel(dmceystr)
maintitle = dataname + date + newline;
title(maintitle + sprintf("Sweeping up in field at %i Oe/s",slow*10^4))

%%
figure
plot(SlowDown.newH, SlowDown.dT*50 + SlowDown.TPuck,'.')
xlabel('Field (T)')
ylabel(dmceystr)
title(maintitle + sprintf("Sweeping down in field at %i Oe/s",slow*10^4))

%% Plot MCE traces
figure; hold on
for i = 1:2
plot(SweepSel(i).H, SweepSel(i).T,'.',...
    'DisplayName',sprintf('%i Oe/s', SweepSel(i).sweeprate*10^4))
end
xlabel('Field (T)')
ylabel(mceystr)
title(maintitle + "Comparing up and down sweeps")
legend('show')

%% Plot 1st derivative of MCE traces
figure
plot(FastDown.newH, 0.1*FastDown.dT ./ FastDown.dH + FastDown.TPuck,'.',...
'DisplayName',sprintf('%i Oe/s down',fast*10^4))
hold on
plot(FastUp.newH, 0.1*FastUp.dT ./ FastUp.dH + FastUp.TPuck,'.',...
'DisplayName',sprintf('%i Oe/s up',fast*10^4))
plot(SlowDown.newH, 0.1*SlowDown.dT ./ SlowDown.dH + SlowDown.TPuck,'.',...
'DisplayName',sprintf('%i Oe/s down',slow*10^4))
plot(SlowUp.newH, 0.1*SlowUp.dT ./ SlowUp.dH + SlowUp.TPuck,'.',...
'DisplayName',sprintf('%i Oe/s up',slow*10^4))
xlabel('Field (T)')
ylabel(dmceystr)
title(maintitle + "Comparing up and down sweeps")
legend('show')

%% Plot 2nd derivative of MCE traces
figure; hold on
for i = 1:2
% plot(SweepSel(i).Hmid, SweepSel(i).TPuckMid,'.',...
% plot(SweepSel(i).Hmid, SweepSel(i).d2T ./ (10*SweepSel(i).dHmid).^2 + SweepSel(i).TPuckMid,'.',...
plot(SweepSel(i).Hmid, SweepSel(i).d2T*1000 + SweepSel(i).TPuckMid,'.',...
    'DisplayName',sprintf('%i Oe/s',SweepSel(i).sweeprate*10^4))
end
xlabel('Field (T)')
ylabel(d2mceystr)
ylim([0.35 .7])
title(maintitle + "Comparing up and down sweeps")
legend('show')
