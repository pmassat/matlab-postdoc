function S = importfielddistrib_csv_old(filename, varargin)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   MAGFIELDDISTRIB = IMPORTFILE(FILENAME)
%   Reads data from text file FILENAME for the default selection.
%
%   MAGFIELDDISTRIB = IMPORTFILE(FILENAME,
%   STARTROW, ENDROW) Reads data from rows STARTROW through ENDROW of text
%   file FILENAME.
%
% Example:
%   magfielddistrib = importfile('2020-07-21_TmVO4-RF-E_zero-temp-magnetization_mag-field_distrib.txt', 1, 13671);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2020/07/21 18:49:31

%% Parse function arguments
    defaultNumFields = 1;
    defaultDomainNum = 2;
    defaultStartRow = 10;
    defaultEndRow = inf;
    validLiteral = @(x) isstring(x) || ischar(x);
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    
    p = inputParser;
    addRequired(p, 'filename', validLiteral);
    addOptional(p, 'numfields', defaultNumFields, validScalarPosNum);
    addParameter(p, 'startRow', defaultStartRow, validScalarPosNum);
    addParameter(p, 'endRow', defaultEndRow, validScalarPosNum);
    addParameter(p, 'domainNum', defaultDomainNum, validScalarPosNum);

    parse(p, filename, varargin{:});
    
    numfields = p.Results.numfields;
    startRow = p.Results.startRow;
    endRow = p.Results.endRow;
    domainNum = p.Results.domainNum;

%% Initialize variables.
% if nargin<=3
%     startRow = 10;
%     endRow = inf;
% end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
% formatSpec = '%6s%16s%15s%22s%11s%12s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%s%[^\n\r]';
data_format = repmat('%s',[1,numfields+3]);
formatSpec = strcat(data_format,'%[^\n\r]');

%% Open the text file.
fileID = fopen(filename,'r');

%% (Edit) Read header
headerFormatSpec = '%s';
header = textscan(fileID, headerFormatSpec, numfields+1, 'TextType', 'string',...
    'Delimiter', {sprintf(',if(dom==%d,mfnc.Hz,0) (Oe) @ ', domainNum)},...
    'HeaderLines', startRow(1)-2, 'ReturnOnError', false, 'EndOfLine', '\r\n');
frewind(fileID);
hdr = cell(1,numfields);% header contains field values
for col=1:numfields
    hdr{col} = header{1}(col+1);
end

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1,...
    'Delimiter', ',', 'WhiteSpace', '', 'TextType', 'string',...
    'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1,...
        'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines',...
        startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=1:numfields+3
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Exclude rows with non-numeric cells
% I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),2); % Find rows with non-numeric cells
% raw(I,:) = [];
I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),1); % Find columns with non-numeric cells
raw(:,I) = {NaN};

%% Create output variable
magfielddistrib = cell2mat(raw);

%% Combine field distribution header and array into structure
for col=length(hdr):-1:1
   S(col).T_Bext = hdr{col};
   S(col).mfd = magfielddistrib(:,col+3);% First 3 columns contain space coordinates, which we are not interested in here
end
end
