%% Batch rename files in current folder
% Get all files with relevant name pattern in current folder
filesrn = dir('HH*');
% Loop through each
for id = 1:length(filesrn)
%     and rename them
      movefile(filesrn(id).name, strcat(filesrn(id).name,".txt"));
end
%% Exctract field from multi-dimensional structure to a matrix
% <https://www.mathworks.com/matlabcentral/answers/268762-multi-struct-to-matrix 
% https://www.mathworks.com/matlabcentral/answers/268762-multi-struct-to-matrix>
% 
% Also implemented in 'ENS_peak_fit_ICpV.mlx'
% This does pretty much the same as the extractfield() function, but is
% also applicable to fitresults stored in structures, which the
% extractfield() function cannot do
s(1).x = [1 2 3; 4 5 6; 7 8 9]; s(1).x 
s(2).x = [11 12 13; 14 15 16; 17 18 19]; s(2).x 
cell2mat( arrayfun(@(c) c.x(2,:), s(1:length(s)).', 'Uniform', 0) )
% Note: on the Matlab Answers webpage, the previous line of code is missing the 's' before in "s(1:length(s)).'"
% Output: ans = [4 5 6; 14 15 16]

%% Count the number of elements in a field 'f1' of a structure
% https://www.mathworks.com/matlabcentral/answers/841-how-to-find-number-of-arrays-in-a-structure-filed
% Also implemented in 'ENS_peak_fit_ICpV.mlx'
function N1 = fieldCount(inputStruct)
N1 = 0;
for i = 1: numel(inputStruct)
  if(~isempty(inputStruct(i).f1)))
      N1 = N1 + 1;
  end
end

%% Check if a variable exists in a table column
% https://www.mathworks.com/matlabcentral/answers/313776-how-to-check-whether-a-column-exist-in-a-table
if any(strcmp('Variable',(tableName).Properties.VariableNames));% if there is a field called 'Variable' in table(tableName)
    % do whatever
end

%% Surface plot and colorbar handling
% see 'ENS_peak_fit_ICpV.m'

%% Loop over the fieldnames of a structure
% https://www.mathworks.com/matlabcentral/answers/341454-how-to-loop-over-the-structure-fields-and-get-the-type-of-data
fn = fieldnames(mystruct);
for k=1:numel(fn)
    if( isnumeric(mystruct.(fn{k})) )
        % do stuff
    end
end

%% Flip the values of a colormap for negative intensities
colormap(flipud(jet));% in this case, color map jet
% see https://www.mathworks.com/matlabcentral/answers/103691-how-can-i-invert-the-distribution-of-colors-in-a-colormap-in-matlab-8-1-r2013a?s_tid=gn_loc_drop

%% Delete files or objects
% https://www.mathworks.com/help/matlab/ref/delete.html
delete filename
delete filename1 ... filenameN
delete(obj)

%% Custom legend for two plots called p1 and p2
legend([p1,p2],'p1','p2','Location','best');








