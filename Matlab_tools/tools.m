%% Edit Matlab path manually
% The list of directories that Matlab uses is defined in the file
% 'pathdef.m'

%% Log Command Window text to file
% diary
% diary filename
% diary off
% diary on

%% Running for loop containing subsections
% In order to run a for loop that contains subsections like this
% one, the cursor has to be outside the loop when executing the run command
for idx=1:4
    stuff = idx
    %% Subsection
    more_stuff = -idx
end

%% Parse function arguments
% See https://www.mathworks.com/help/matlab/ref/inputparser.html

%% Get today's date and convert it to string
todaystr = datestr(datetime('today'),29);% datetime outputs a 'datetime' variable type, not string

%% Read file and replace regular expression
test1 = fileread('filename.m');% Read file and store its content as string into variable test1
str = regexprep(test1,'[\\\n\r]+','\n');% Replace all strings with a '\' sign followed by newlines with a single newline

%% Save string to file
fileID = fopen('fexp.m','w');
fprintf(fileID,str);
fclose(fileID);

%% Generate array of...
e = 1.2e-3*(1:2:9);%... arithmetic sequence
e = 1.5e-3*2.^(0:4);%... geometric sequence

%% Batch rename files in current folder
filesrn = dir('Dixon1980 (*).JPG');% Get all files with relevant name pattern in current folder
for id = 1:length(filesrn)% Loop through each item of the file structure filesrn
%     movefile(filesrn(id).name, strcat(filesrn(id).name,".txt"));% add string at the end of filename
    movefile(filesrn(id).name,regexprep(filesrn(id).name,' \((\d+)\)','_$1'));
    % replace Windows batch naming system (i.e. 'filename (filenumber).extension') with
    % underscore file numbering, i.e. 'filename_filenumber.extension'.
% Details of regular expression: 
% * the initial space identifies the space in the Windows filename
% * '\(' is the string that identifies a bracket character (same for '\)')
% * '\d+' identifies  a sequence of one or more digits inside the brakets
% * the extra pair of brackets identifies the sequence of digits as a token
%   that is then used in the replacing expression using the string '$1',
%   which guarantees that the final filename will have the same file number
%   as the initial one.
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

%% Check if a variable exists in a table column
% https://www.mathworks.com/matlabcentral/answers/313776-how-to-check-whether-a-column-exist-in-a-table
% if any(strcmp('Variable',(tableName).Properties.VariableNames));% if there is a field called 'Variable' in table(tableName)
%     % do whatever
% end

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

%% Annotation
annttl = annotation('textbox',[0.55 0.775 0.2 0.1],'interpreter','latex',...
    'String',{'Tm$_{1-x}$Y$_x$VO$_4$'},'LineStyle','-','EdgeColor','k',...
    'FitBoxToText','on','LineWidth',1,'BackgroundColor','w','Color','k');% add annotation

%% Find all text objects in figure
a = findall(gcf,'Type','text');

%% Call a function as argument of another function
f1 = @(x,y) x(y) ;
f2 = @(y) sin(y) ;
f1(f2,pi/2)

%% Differentiate a matrix with respect to n-th dimension
% Matlab doc: https://www.mathworks.com/help/matlab/ref/diff.html
% Y = diff(X,n,dim) is the nth difference calculated along the dimension specified by dim. The dim input is a positive integer scalar.
X = [1 3 5;7 11 13;17 19 23];
Y = diff(X,1,1);

%% Default line colors array is called lines
% Calling lines(n) returns a n by 3 array of RGB colors
cmap = lines(n);% Returns a nx3 array 
% then call color in j-th row as cmap(j,:)
% See https://www.mathworks.com/help/matlab/ref/colormap.html for more
% colormaps

%% Change the line style order for plots
% Note: the default line style order is {'-'} i.e. only one line style
% See Matlab doc: https://www.mathworks.com/help/matlab/ref/matlab.graphics.axis.axes-properties.html#budumk7_sep_shared-LineStyleOrder
ax = axes('LineStyleOrder',{'-',':','--','-.'});%

%% Find light object
l = findobj(gcf,'Type','Light');

%% Multi-line string
% Use command "newline"
['two line' newline 'test']

%% Sequence "derivative" 
% Such that the "derivative" is calculated at the same point as in the
% original dataset
diffT = diff(T);
dT = zeros(size(T));
dT(2:end-1) = (diffT(2:end)+diffT(1:end-1))/2;
dT(1) = diffT(1)/2; dT(end) = diffT(end)/2;

%% Lighting properties 
% https://www.mathworks.com/help/matlab/creating_plots/lighting-overview.html
% Light intensity set by the following properties of surface and patch objects:
% AmbientStrength, DiffuseStrength, SpecularStrength

%% Plot fit result on arbitrary scale
% Rather than using plot(cfitobj,...), evaluate the cfit objects on the data range you care about:
xr = x(x>xmin & x<xmax);
y1 = cfit1(xr);
y2 = cfit2(xr);
plot(xr,y1,xr,y2);
% See https://www.mathworks.com/matlabcentral/answers/161244-how-to-plot-a-cfit-object-in-the-borders-of-the-original-data

%% Count the number of elements in a structure field
% https://www.mathworks.com/matlabcentral/answers/841-how-to-find-number-of-arrays-in-a-structure-filed
% Also implemented in 'ENS_peak_fit_ICpV.mlx'
function N1 = fieldCount(inputStruct)
N1 = 0;
for i = 1: numel(inputStruct)
  if(~isempty(inputStruct(i).f1))% field is called f1 here
      N1 = N1 + 1;
  end
end
end
