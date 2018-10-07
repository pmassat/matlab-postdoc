function plotCpField(mycellarray,lgdarray,varargin)
%% Optional arguments
nvarargs = length(varargin);
optargs = {'Heat capacity of TmVO_4', 'C_p (J\cdotmol^{-1}\cdotK^{-1})'};
optargs(1:nvarargs) = varargin;
[ttl, ylbl] = optargs{:};% The optional argument defines the label of the y axis of the plot

%% Actual plot
close 
figure
for i=1:length(mycellarray)
    plot(mycellarray{i}.T,mycellarray{i}.Cpmol,'.','MarkerSize',18,'DisplayName',[num2str(lgdarray(i)),' Oe'])
    hold on
end
xlabel('Temperature (K)')
ylabel(ylbl)
title(ttl)
legend('show')