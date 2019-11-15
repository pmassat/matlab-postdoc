% 
% figure;
% ax = axes('LineStyleOrder',{'-*',':','--o'});
% hold(ax,'on');
% for i=1:3
%     ax.ColorOrderIndex = i;
%     ax.LineStyleOrderIndex = i;
%     plot(ax,1+10*i:10+10*i,rand(1,10))
% end

fileID = fopen('fexp.m','w');
fprintf(fileID,str);
fclose(fileID);
