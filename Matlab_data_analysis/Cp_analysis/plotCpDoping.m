function plotCpDoping(mycellarray,lgdarray,varargin)
%% Optional arguments
nvarargs = length(varargin);
optargs = {1,length(mycellarray),'Heat capacity of TmVO$_4$', 'C$_p$ (J$\cdot$mol$^{-1}\cdot$K$^{-1}$)'};
optargs(1:nvarargs) = varargin;
[istart,iend,ttl, ylbl] = optargs{:};% The optional argument defines the label of the y axis of the plot

%% Actual plot
close 
figure
for i=istart:iend
    plot(mycellarray{i}.T,mycellarray{i}.Cpmol,'.','MarkerSize',18,'DisplayName',['x = ',num2str(lgdarray(i))])
    hold on
end
xlabel('Temperature (K)')
ylabel(ylbl)
title(ttl)
legend('show')