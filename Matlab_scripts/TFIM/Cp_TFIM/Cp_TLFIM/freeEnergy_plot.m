%% Change to relevant directory
cd C:\Users\Pierre\Desktop\Postdoc\Software\Matlab\Matlab_data_analysis\TFIM\Cp_TFIM

%%
h = [0.,0.6,0.9];% transverse field
e = [1e-5,1e-3,1e-2];% longitudinal field
T = linspace(1e-3,2,2000);
% F = zeros(length(T),length(h),length(e));
F = freeEnergy_compute(T,h,e);% free energy

%%
TLFIM_plot(T,F,h,e);
title('Free energy');

%% Update the following...

S = Entropy_TLFIM(T,F);
%%
TLFIM_plot(T(2:end),S,h,e);
title('Entropy')

%%
Cp = Cp_TLFIM(T,S);

%%
TLFIM_plot(T(2:end-1),Cp,h,e);
title({'$C_p$ in transverse ($h$) and' 'longitudinal ($e$) fields Ising model'})
xlabel('$T (K)$')
ylabel('$C_p/R$')

%% Export figure
% formatFigure;
printPDF('2019-07-26_Cp_TFIM_vs_t_h_e')
