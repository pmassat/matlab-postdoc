function printSVG(filestr)
    h = gcf;% get current figure
%     set(h,'Units','Inches');
% %     pos = get(h,'Position');
% %     set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%     dim = [6 5];% dimensions, in the units set above
    dim = get(h,'Position');% match figure dimensions
    h.PaperSize = dim(3:4);% with paper dimensions
%     InSet = get(gca, 'TightInset');% and expand axes to fill figure
% %     set(gca, 'Position', [InSet(1:2), 0.85-InSet(1)-InSet(3), 0.94-InSet(2)-InSet(4)]);% for color maps
% % Note: this line resizes the figure to make it fill the page; however, in
% % the case of a color map, it does not take into account the colorbar, thus
% % one needs to adjust the values of elements 3 and 4 manually
%     set(gca, 'Position', [InSet(1:2), 0.99-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);% for all other plots
% %     print(h,'-fillpage',filestr,'-dpdf','-r0');
% %     set(h,'Color','none');% set figure background color to 'none'
    h.Color = 'none';
    print(h,filestr,'-dsvg','-r0');
%       saveas(h,filestr,'pdf');
    h.Color = 'w';% set figure background color back to 'white' for use in Matlab
end