% for i=1:4
%     for j=3:6
%     if (j<5 && i>2) || i<=2
%         disp([i j])
%     end
%     end
% end
x = [0,.1,.2,.3];
y = [1,.8,.5,0];
figure
area(x,y,'FaceColor',[0 .5 0]);% dark green [0 .5 0]; light green [0 1 0]
xlim([])