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
% Also implemented in ENS_peak_fit_ICpV.mlx
%%
s(1).x = [1 2 3; 4 5 6; 7 8 9]; s(1).x 
s(2).x = [11 12 13; 14 15 16; 17 18 19]; s(2).x 
cell2mat( arrayfun(@(c) c.x(2,:), s(1:length(s)).', 'Uniform', 0) )
% Note: on the Matlab Answers webpage, the previous line of code is missing the 's' before in "s(1:length(s)).'"
% Output: ans = [4 5 6; 14 15 16]