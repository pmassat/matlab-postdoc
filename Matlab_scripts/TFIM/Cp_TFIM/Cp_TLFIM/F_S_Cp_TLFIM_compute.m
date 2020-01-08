%% Change to relevant directory
cd C:\Users\Pierre\Desktop\Postdoc\Software\Matlab\Matlab_data_analysis\TFIM\Cp_TFIM

%%
% factor = 0.97/0.92;
% h = [0.,0.55,0.69,0.84,0.92]*factor;% transverse field
% e = 1.5e-3*(1+7*h.^3);% longitudinal field
e = 1.5e-3*(1+0*h);% longitudinal field
Ttlf = linspace(1e-2,1.5,150);
% F = zeros(length(T),length(h),length(e));

%% Compute free energy at a given field
F = freeEnergy_TLFIM_compute(Ttlf,h,e);% free energy

%% Compute free energy for a normal distribution of fields
Fn = freeEnergy_TLFIM_normpdf(Ttlf(10),h,0,rhsgm*h);

%% Compute entropy
S = Entropy_TLFIM(Ttlf,F);
Sn = Entropy_TLFIM(Ttlf,Fn);

%%
Cptlf = Cp_TLFIM_array(Ttlf,S);

%% Plot free energy
TLFIM_plot(Ttlf,F,h,e);
title('Free energy');

%% Plot entropy
TLFIM_plot(Ttlf(2:end),S,h,e);
title('Entropy')

%% Plot Cp
TLFIM_plot(Ttlf(2:end-1),Cptlf,h,e);
title({'$C_p$ in transverse ($h$) and' 'longitudinal ($e$) fields Ising model'})
xlabel('$T (K)$')
ylabel('$C_p/R$')

%% Export figure
% formatFigure;
printPDF('2019-07-26_Cp_TFIM_vs_t_h_e')
