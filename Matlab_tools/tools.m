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

