function printPDF(filestr)
    h = gcf;% get current figure
    set(h,'Units','Inches');
%     pos = get(h,'Position');
%     set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    set(h,'Color','none');% set figure background color to 'none'
    dim = [6 5];% dimensions, in the units set above
    h.Position(3:4) = dim;% match figure dimensions
    h.PaperSize = dim;% with paper dimensions
    InSet = get(gca, 'TightInset');% and expand axes to fill figure
    set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
%     print(h,'-fillpage',filestr,'-dpdf','-r0');
      print(h,filestr,'-dpdf','-r0');
    h.Color = 'w';% set figure background color back to 'white' for use in Matlab
end