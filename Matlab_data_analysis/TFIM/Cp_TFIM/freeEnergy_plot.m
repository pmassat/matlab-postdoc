%% Change to relevant directory
cd C:\Users\Pierre\Desktop\Postdoc\Software\Matlab\Matlab_data_analysis\TFIM\OP_TFIM

%%
h = [0.,0.6];%,0.9];% transverse field
e = 1e-5;%[1e-5,1e-3]% longitudinal field
T = linspace(1e-3,2,1999);
% F = zeros(length(T),length(h),length(e));
F = freeEnergy_compute(T,h,e);% free energy

%%
TLFIM_plot(T,F,h,e);

%% Update the following...
% %%
% figure; hold on
% plot(T,F)
% title('Free energy')
% %%
% diff1f = diff(F);
% dT = diff(T);
% S = -diff1f./dT;% entropy
% %%
% figure; hold on
% plot(T(2:end),S)
% title('Entropy')
% %%
% diff2f = diff(S);
% Cp = T(3:end).*diff2f./dT(2:end);
% %% 
% figure; hold on
% plot(T(3:end),Cp)
% title('Cp')
% 
% %% Export figure
% formatFigure;
% printPDF('2019-07-02_OP_TFIM_vs_t_h_e')
