%% Change to relevant directory
cd C:\Users\Pierre\Desktop\Postdoc\Software\Matlab\Matlab_scripts\TFIM\Cp_TFIM

%% Initialize variables for longitudinal field dependence 
e = logspace(-3,-1,3);% longitudinal field
t = linspace(1e-2,1.5,150);
F = zeros(length(e),length(t));
S = diff(F,1,2);
Cp = diff(S,1,2);

%% Compute free energy for various values of longitudinal field
he = 0;
for ie=1:length(e)
    [F(ie,:),S(ie,:),Cp(ie,:)] = FSCp_TLFIM(t,he,e(ie));
end

%% Plot Cp
figure; hold on
for ie=1:length(e)
    plot(t(2:end-1),Cp(ie,:),'DisplayName',sprintf('%.2g',e(ie)))
end
xlabel('$T$ (K)')
ylabel('$C_p/R$')
title(['Heat capacity in the TFIM, $h=$ ' sprintf('%.2g',0)])
lgd = legend('Location','best');
lgd.Title.String = '$e$';


%% Initialize variables for longitudinal field dependence 
factor = 0.97/0.92;
h = [0.,0.55,0.69,0.84,0.92]*factor;% transverse field
% Fh = zeros(length(e),length(t));
% Sh = diff(F,1,2);
% Cphe = diff(S,1,2);

%% Compute free energy for various values of transverse field
eh = e(1);
% for ih=1:length(h)
%     [Fh(ih,:),Sh(ih,:),Cphe(ih,:)] = FSCp_TLFIM(t,h(ih),eh);
% end
[Fh,Sh,Cphe] = FSCp_TLFIM(t,h,eh);

%% Plot Cp
figure; hold on
for ih=1:length(h)
    plot(t(2:end-1)',Cphe(:,ih),'DisplayName',sprintf('%.2g',h(ih)))
end
xlabel('$T$ (K)')
ylabel('$C_p/R$')
title(['Heat capacity in the TFIM, $e=$ ' sprintf('%.2g',e(1))])
lgd = legend('Location','best');
lgd.Title.String = '$h$';

%% Export figure
formatFigure;
% printPNG([todaystr '_FSCp_TLFIM_e=1e-3']);










%% Outdated as of 2020-08-06

%% Compute free energy for various values of longitudinal strains
% F = freeEnergy_TLFIM_compute(t,h,e);% free energy
%% Compute free energy for a normal distribution of fields
% Fn = freeEnergy_TLFIM_normpdf(Ttlf(10),h,0,rhsgm*h);

%% Compute entropy
S = Entropy_TLFIM(t,F);
% Sn = Entropy_TLFIM(Ttlf,Fn);

%%
Cptlf = Cp_TLFIM_array(t,S);

%% Plot free energy
TLFIM_plot(t,F,h,e);
title('Free energy');

%% Plot entropy
TLFIM_plot(t(2:end),S,h,e);
title('Entropy')

%% Plot Cp
TLFIM_plot(t(2:end-1),Cptlf,h,e);
title({'$C_p$ in transverse ($h$) and' 'longitudinal ($e$) fields Ising model'})
xlabel('$T (K)$')
ylabel('$C_p/R$')

%% Export figure
% formatFigure;
printPDF('2019-07-26_Cp_TFIM_vs_t_h_e')
