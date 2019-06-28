function ht = CbTixAbsolute(ax,cb,xoffset,yoffset)
%% Position tick labels absolutely
% ax is the axis containing the colormap
%   use ax = gca prior to calling function for current axis
% cb is the colorbar of a 2D colormap
% xoffset is a user-defined x axis offset; choose 0 for no offset 
% yoffset is a user-defined y axis offset; choose 0 for no offset 
% 
%% Define variables
cb.TickLabels = {};% Delete default tick labels
%%Grab some position properties for later use
ax_pos = get(ax,'position');% relative axes coordinates and dimensions in figure
ax_limx = get(ax,'XLim');% x axis limits
ax_limy = get(ax,'YLim');% y axis limits
cb_pos = get(cb,'position');
cb_lim = get(cb,'limits');% extremal values of colorbar
%%Grab colorbar ticks
tix=cb.Ticks;% colorbar tick values
%%Create cell array with ticklabels
tixlabels=sprintfc('%.4g',tix)';% convert tix to array of strings
% 
%% Normalized coordinates for ticks, i.e. relative position in figure
if cb_pos(4)>cb_pos(3)% Vertical colorbar: height is bigger than width
    yi = cb_pos(2) + cb_pos(4)*(tix(1)-cb_lim(1))/(cb_lim(2)-cb_lim(1));% y coordinate of first tick
    yf = cb_pos(2)+cb_pos(4) - cb_pos(4)*(cb_lim(2)-tix(end))/(cb_lim(2)-cb_lim(1));% coordinate of last tick
    ny = linspace(yi,yf,numel(tix))';% y coordinates of cb tick labels in relative units
    nx = repmat(cb_pos(1)+cb_pos(3),numel(tix),1);% remove cb_pos(3) in order to position tick labels to the left of the colorbar
else% Horizontal colorbar: height is smaller than width
    xi = cb_pos(1) + cb_pos(3)*(tix(1)-cb_lim(1))/(cb_lim(2)-cb_lim(1));% x coordinate of first tick
    xf = cb_pos(1)+cb_pos(3) - cb_pos(3)*(cb_lim(2)-tix(end))/(cb_lim(2)-cb_lim(1));% coordinate of last tick
    nx = linspace(xi,xf,numel(tix))';% x coordinates of cb tick labels in relative units
    ny = repmat(cb_pos(2)+cb_pos(4),numel(tix),1);% remove cb_pos(4) in order to position tick labels under the colorbar
end
% 
%% Convert from normalized (nx,ny) to axes (x,y) coordinates, i.e. position in axes coordinates
% This is necessary because the 'text' function takes axes coordinates as input
ky = diff(ax_limy)/ax_pos(4);% ratio between x axis length in x axis units and x axis length in relative units
kx = diff(ax_limx)/ax_pos(3);% ratio between y axis length in y axis units and y axis length in relative units
y = ax_limy(1) + ky.*(ny-ax_pos(2));% y axis coordinates of cb tick labels starting from beginning of y axis
x = ax_limx(1) + kx.*(nx-ax_pos(1));% x axis coordinates of cb tick labels starting from beginning of x axis
% 
%% Write out ticklabels
if cb_pos(4)>cb_pos(3)% Vertical colorbar
    ht=text(x+xoffset,y+yoffset,tixlabels,...
        'rotation',0,'verticalalignment','middle','horizontalalignment','left');
else% Horizontal colorbar
    ht=text(x+xoffset,y+yoffset,tixlabels,...
        'rotation',0,'verticalalignment','bottom','horizontalalignment','center');
end
end