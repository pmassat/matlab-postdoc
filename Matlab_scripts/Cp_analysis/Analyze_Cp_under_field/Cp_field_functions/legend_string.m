function lgd_str = legend_string(param_index, p1, p2)
% create legend string based on parameter index param_index, which determines whether
% to plot MFD's at various fields at constant temperature, or at various
% temperatures and constant field; 
% p1 and p2 can be the variable names for
% either magnetic field or temperature.
if param_index==1
    lgd_str = sprintf('%.2g', p1);
elseif param_index==2
    lgd_str = sprintf('%.2g', p2);
end
