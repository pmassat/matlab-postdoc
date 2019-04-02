function y = random_strains_phase_boundary(delta)
tc = @(t) random_strains_phase_boundary_equation(delta,t);
y = fzero(tc,[-1 2]);% when tc goes to zero, fzero needs to go to negative values to output a value

% %% Plot phase boundary
% td = @(d) random_strains_phase_boundary(d);
% figure;
% fplot(@(d)td(d),[1e-3 1.25],'LineWidth',2)
% title("T$_D$ vs spread in strain distribution");
% xlabel('$\Delta_0$/($x$(Tm)$\cdot T_D$($\Delta_0$=0))'); 
% ylabel('$T_D$/($x$(Tm)$\cdot T_D$($\Delta_0$=0))');

end