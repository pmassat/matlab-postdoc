function plot_Cp_avg_w_fits(rng, avgData, Cpnum, Tc, label, varargin)

%% Parse function arguments
    defaultTnumStr = 't_single_h_w_e';
    defaultCpnumStr = 'single_h_w_e';
%     validLiteral = @(x) isstring(x) | ischar(x);
    
    p = inputParser;
    addParameter(p, 'TnumStr', defaultTnumStr);
    addParameter(p, 'CpnumStr', defaultCpnumStr);

    parse(p, varargin{:});
    
    TnumStrstr = p.Results.TnumStr;
    CpnumStrstr = p.Results.CpnumStr;

%% Actual function
figure; 
hold on
range = [1,rng];
clr = lines(length(range));
eb = cell(size(range));
maxTplot = 3.2;%
R = 8.314;

plot(Cpnum(1).t_single_h_w_e*Tc, Cpnum(1).single_h_w_e);

for i=range(2:end)
	plot(Cpnum(i).(TnumStrstr)*Tc, Cpnum(i).(CpnumStrstr), 'Color', clr(range==i,:));
%     fp = fplot(@(t) Cp_TFIM(t/Tc0rf,uhrf(i)*rescaling/(Hc0*1e4)),[0 maxTplot],'LineWidth',2);
end

for i=range
    eb{range==i} = errorbar(avgData(i).T, avgData(i).Cpelr, avgData(i).CpFullErr/R,...
        '.','MarkerSize',18,'DisplayName', num2str(label(i),'%.2f'),...
        'Color',clr(range==i,:),'LineWidth',2);
end

xlabel('$T$ (K)'); ylabel('$C_p/R$');%ylabel('C$_p$ (JK$^{-1}$mol$^{-1}$)');
xlim([0 maxTplot])
ylim([0 1.5])
% title('Heat capacity of TmVO4-RF-E')
lgd = legend([eb{:}]); 
lgd.Title.String = '$H/H_{c,0}^{\mathrm{SC}}$';
lgd.FontSize = 22;
lgd.Location = 'northwest';
ax = gca; ax.YMinorTick = 'on';% Add minor ticks on Y axis
ax.FontSize = 22;
grid on;%
hold off
end

