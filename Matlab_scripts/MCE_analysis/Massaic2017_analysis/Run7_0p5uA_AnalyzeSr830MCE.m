cd 'C:\Users\Pierre\Desktop\Postdoc\TmVO4\TmVO4_heat-capacity\2017-07_TmVO4_Cp_MCE\2017-07-28--31\Massaic_MCE\'
load('betamodel.mat')
dataname = "Run7\_0p5uA\_2017-08-03.dat";
date = " 03 Aug 2017";
SR830MCE(7)=ImportSR830MCE('Run7_0p5uA_2017-08-03.dat',5e-7);

%%
figure
plot(SR830MCE(7).H,1./BetaModel(SR830MCE(7).H, SR830MCE(7).ResPlat))% R depends on T and H
% R and H are known in the experiment; we want to know T
% BetaModel results from the fit of 1/T as a function of R and H using the
% calibration file 'thermometer_calibration_raw.txt'
xlabel('Magnetic field (T)')
ylabel('Temperature (K)')
title('MCE in assembly of TmVO4 needles')

%%
sgm=10;%number of points for smoothing
newH=gConvolve(SR830MCE(7).H,sgm)*1e4;%smooth the field data
newT=gConvolve(1./BetaModel(SR830MCE(7).H, SR830MCE(7).ResPlat),sgm);%smooth the temperature data
dT = diff(newT);%compute the derivative of the smoothed temperature
dH = diff(newH);

%%
xstr = '$H$ (Oe)';

%%
figure
plot(newH(30:(end-30)),dT(30:(end-29))*50 + SR830MCE(7).TPuck(30:(end-30)),'-')%plot dT/dH versus H
%ylim([0.3 1.5])
xlabel(xstr)
ylabel('$dT/dH + T_{\mathrm{puck}}$')
title('Derivative of MCE in assembly of TmVO4 needles')

%%
figure;
yyaxis left;
plot(diff(SR830MCE(7).H*1e4),'DisplayName','Raw');
hold on;
plot(diff(newH),'-g','DisplayName','Smoothed');
xlabel('Measurement time $t$ (s)')
ylabel('$dH/dt$ (Oe/s)')
title('Sweep rate and puck temperature')
legend('show')

yyaxis right;
plot(SR830MCE(7).TPuck,'DisplayName','Puck temperature')
ylabel('He3 temperature (K)')

%%
fast = 28;%Fast sweep rate of magnetic field (in Oe/s)
slow = 14;%Slow sweep rate of magnetic field (in Oe/s)
ref_sweep = (fast+slow)/2;%Reference sweep rate of magnetic field (in Tesla/s) used 
%to discriminate between high and low sweep rates
SweepSel(1).sweeprate = fast;
SweepSel(3).sweeprate = slow;
SweepSel(2).sweeprate = -fast;
SweepSel(4).sweeprate = -slow;
SweepSel(1).Filter = diff(newH)>ref_sweep;
SweepSel(3).Filter = diff(newH)>slow/2 &  diff(newH)<ref_sweep;
SweepSel(2).Filter = diff(newH)<-ref_sweep;
SweepSel(4).Filter = diff(newH)<-slow/2 &  diff(newH)>-ref_sweep;

%% Create structure containing smoothed data separated by sweeprate
for i = 1:4
SweepSel(i).H = newH(SweepSel(i).Filter);
SweepSel(i).T = newT(SweepSel(i).Filter);
SweepSel(i).dT = dT(SweepSel(i).Filter);
SweepSel(i).dH = dH(SweepSel(i).Filter);
SweepSel(i).TPuck = SR830MCE(7).TPuck(SweepSel(i).Filter);
SweepSel(i).d2T = diff(SweepSel(i).dT);

% The following lines are wrong since SweepSel(i).Filter already
% selects data in the range [2:end]
% SweepSel(i).Hmid = (SweepSel(i).H(1:end-1) + SweepSel(i).H(2:end))/2;
% SweepSel(i).dHmid = (SweepSel(i).dH(1:end-1) + SweepSel(i).dH(2:end))/2;
% SweepSel(i).TPuckMid = (SweepSel(i).TPuck(1:end-1) + SweepSel(i).TPuck(2:end))/2;
end

%% global variables
ut = unique(round(SweepSel(i).TPuck,1));
mceystr = '$T_{\mathrm{bath}} + \Delta T_{\mathrm{MCE}}$ (K)';
dmceystr = '$T_{\mathrm{bath}} + dT_{\mathrm{MCE}}/dH$';

%% Plot MCE derivative at each sweeprate in individual figures
for i=1:2
    figure
    plot(SweepSel(i).H, 1e3*SweepSel(i).dT./SweepSel(i).dH + SweepSel(i).TPuck, '.')
    xlabel(xstr)
    ylabel(dmceystr)
    maintitle = dataname + date + newline;
    title(maintitle + sprintf("Sweeping magnetic field at %i Oe/s", SweepSel(i).sweeprate))
end

%% Plot 1st derivative of MCE traces
figure; hold on
for i = 1:2
    plot(SweepSel(i).H, 1e3.*SweepSel(i).dT./SweepSel(i).dH + SweepSel(i).TPuck,'.',...
        'DisplayName',sprintf('%i Oe/s down', SweepSel(i).sweeprate))
end
xlabel(xstr)
ylabel(dmceystr)
title(maintitle + "Comparing up and down sweeps")
legend('show')

%%
d2mceystr = ['$T$ (K) $+\, d^2\Delta T/dH^{2}$ (K/kOe$^2$)'];

%% Plot 2nd derivative of MCE traces
figure; hold on
clr = lines(4);

for i = 1:2
SweepSel(i).d2max = zeros(size(ut));
SweepSel(i).hmax = zeros(size(ut));
SweepSel(i).hctd2 = cell(size(ut));
Tpuck = SweepSel(i).TPuck(1:end-1);
dH2 = SweepSel(i).dH(1:end-1);
d2mce = SweepSel(i).d2T ./ (1e-3*dH2).^2 + Tpuck;

% plot(SweepSel(i).Hmid, SweepSel(i).TPuckMid,'.',...
% plot(SweepSel(i).Hmid, SweepSel(i).d2T*1000 + SweepSel(i).TPuckMid,'.',...
p(i) = plot(SweepSel(i).H(1:end-1), d2mce, '.', 'Color', clr(i,:), 'MarkerSize', 12,...
    'DisplayName', sprintf('%i Oe/s',SweepSel(i).sweeprate));

    % Identify the maximum of the second derivative and highlight it on the plot
    for jut = 1:length(ut)
        THselup = Tpuck>ut(jut)-.05 & Tpuck<ut(jut)+.05 &...
                    SweepSel(i).H(1:end-1)>4000 & SweepSel(i).H(1:end-1)<5500;
        SweepSel(i).d2max(jut) = max(d2mce(THselup));
        SweepSel(i).hmax(jut) = SweepSel(i).H(d2mce==SweepSel(i).d2max(jut));
        SweepSel(i).hctd2{jut} = text(SweepSel(i).hmax(jut), SweepSel(i).d2max(jut),...
            sprintf('%.0f',SweepSel(i).hmax(jut)));
        SweepSel(i).hctd2{jut}.FontSize = 16;
        SweepSel(i).hctd2{jut}.VerticalAlignment = 'bottom';
%         SweepSel(i).hctd2{jut}.BackgroundColor = 'w';
        if mod(i,2)==0
            SweepSel(i).hctd2{jut}.HorizontalAlignment = 'left';
        else
            SweepSel(i).hctd2{jut}.HorizontalAlignment = 'right';
        end
    end
    
    plot(SweepSel(i).hmax, SweepSel(i).d2max, '.k', 'MarkerSize', 12)

end

xlabel(xstr)
ylabel(d2mceystr)
ylim([0.35 .68])
ax=gca; ax.YTick = .4:.1:.7;
% title(maintitle + "Comparing up and down sweeps")
legend([p(:)])

%% Plot MCE traces together
figure; hold on

for i = 1:2
    SweepSel(i).hct = cell(size(ut));
    SweepSel(i).td2max = zeros(size(ut));
    p(i) = plot(SweepSel(i).H, SweepSel(i).T, '.', 'Color', clr(i,:), 'MarkerSize', 12, ...
        'DisplayName', sprintf('%i Oe/s', SweepSel(i).sweeprate));

    % Highlight the critical field, as determined from the second derivative
    for jut = 1:length(ut)
        SweepSel(i).td2max(jut) = SweepSel(i).T(SweepSel(i).H==SweepSel(i).hmax(jut));
        SweepSel(i).hct{jut} = text(SweepSel(i).hmax(jut), SweepSel(i).td2max(jut),...
            sprintf('%.0f',SweepSel(i).hmax(jut)));
        SweepSel(i).hct{jut}.FontSize = 16;
        if mod(i,2)==0
            SweepSel(i).hct{jut}.VerticalAlignment = 'top';
        else
            SweepSel(i).hct{jut}.HorizontalAlignment = 'right';
            SweepSel(i).hct{jut}.VerticalAlignment = 'bottom';
        end
    end

    plot(SweepSel(i).hmax, SweepSel(i).td2max, '.k', 'MarkerSize', 12)
end

xlabel(xstr)
ylabel(mceystr)
% title(maintitle + "Comparing up and down sweeps")
legend([p], 'Location', 'southeast')
ax=gca; ax.YTick = .4:.1:.7;

%%
% formatFigure
% printPDF([todaystr '_TmVO4-needles_MCE'])
% printPDF([todaystr '_TmVO4-needles_MCE_d2T-dH'])

