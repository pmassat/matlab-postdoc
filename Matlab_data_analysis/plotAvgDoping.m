function plotAvgDoping(x,lgdarray,varargin)
%% Optional arguments
nvarargs = length(varargin);
optargs = {'Heat capacity of TmVO_4',true,'C_p (J\cdotmol^{-1}\cdotK^{-1})'};
optargs(1:nvarargs) = varargin;
[ttl,yerrbar,ylbl] = optargs{:};% The optional argument defines the label of the y axis of the plot

%% Actual plot
close 
figure
if yerrbar==true
for i=1:length(mystruct)
    errorbar(mystruct(i).T,mystruct(i).Cp,mystruct(i).stdCp,'.','MarkerSize',18,'DisplayName',['x = ',num2str(lgdarray(i))])
    hold on
end
else
    for i=1:length(mystruct)
        plot(mystruct(i).T,mystruct(i).Cp,'.','MarkerSize',18,'DisplayName',['x = ',num2str(lgdarray(i))])
        hold on
    end
end
xlabel('Temperature (K)')
ylabel(ylbl)
title(ttl)
legend('show')