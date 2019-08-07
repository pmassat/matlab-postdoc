function [TimeStampSeconds,Comment,SystemStatusCode,PuckTempKelvin,SystemTempKelvin,FieldOersted,PressureTorr,SampleTempKelvin,TempRiseKelvin,SampHCmJmoleK,SampHCErrmJmoleK,AddendaHCJK,AddendaHCErrJK,TotalHCJK,TotalHCErrJK,FitDeviationChiSquare,TimeConsttau1seconds,TimeConsttau2seconds,SampleCouplingPercent,DebyeTempKelvin,DebyeTempErrKelvin,CalCorrectionFactor,ThermResistOhms,HtrResistOhms,PuckResistOhms,WireCondWK,MeasTimeseconds,TempSquaredK2,SampHCTempmJmoleKK,AddendaOffsetHCJK] = importCpSharedPPMS(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [TIMESTAMPSECONDS,COMMENT,SYSTEMSTATUSCODE,PUCKTEMPKELVIN,SYSTEMTEMPKELVIN,FIELDOERSTED,PRESSURETORR,SAMPLETEMPKELVIN,TEMPRISEKELVIN,SAMPHCMJMOLEK,SAMPHCERRMJMOLEK,ADDENDAHCJK,ADDENDAHCERRJK,TOTALHCJK,TOTALHCERRJK,FITDEVIATIONCHISQUARE,TIMECONSTTAU1SECONDS,TIMECONSTTAU2SECONDS,SAMPLECOUPLINGPERCENT,DEBYETEMPKELVIN,DEBYETEMPERRKELVIN,CALCORRECTIONFACTOR,THERMRESISTOHMS,HTRRESISTOHMS,PUCKRESISTOHMS,WIRECONDWK,MEASTIMESECONDS,TEMPSQUAREDK2,SAMPHCTEMPMJMOLEKK,ADDENDAOFFSETHCJK]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [TIMESTAMPSECONDS,COMMENT,SYSTEMSTATUSCODE,PUCKTEMPKELVIN,SYSTEMTEMPKELVIN,FIELDOERSTED,PRESSURETORR,SAMPLETEMPKELVIN,TEMPRISEKELVIN,SAMPHCMJMOLEK,SAMPHCERRMJMOLEK,ADDENDAHCJK,ADDENDAHCERRJK,TOTALHCJK,TOTALHCERRJK,FITDEVIATIONCHISQUARE,TIMECONSTTAU1SECONDS,TIMECONSTTAU2SECONDS,SAMPLECOUPLINGPERCENT,DEBYETEMPKELVIN,DEBYETEMPERRKELVIN,CALCORRECTIONFACTOR,THERMRESISTOHMS,HTRRESISTOHMS,PUCKRESISTOHMS,WIRECONDWK,MEASTIMESECONDS,TEMPSQUAREDK2,SAMPHCTEMPMJMOLEKK,ADDENDAOFFSETHCJK]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [TimeStampSeconds,Comment,SystemStatusCode,PuckTempKelvin,SystemTempKelvin,FieldOersted,PressureTorr,SampleTempKelvin,TempRiseKelvin,SampHCmJmoleK,SampHCErrmJmoleK,AddendaHCJK,AddendaHCErrJK,TotalHCJK,TotalHCErrJK,FitDeviationChiSquare,TimeConsttau1seconds,TimeConsttau2seconds,SampleCouplingPercent,DebyeTempKelvin,DebyeTempErrKelvin,CalCorrectionFactor,ThermResistOhms,HtrResistOhms,PuckResistOhms,WireCondWK,MeasTimeseconds,TempSquaredK2,SampHCTempmJmoleKK,AddendaOffsetHCJK] = importfile('20180322_TmVO4-LS5228-MP3-Plt-HC1803_Cp.dat',23, 270);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/03/27 15:09:20

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 23;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: text (%s)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: text (%s)
%	column24: text (%s)
%   column25: text (%s)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
%	column30: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
TimeStampSeconds = dataArray{:, 1};
Comment = cellstr(dataArray{:, 2});
SystemStatusCode = dataArray{:, 3};
PuckTempKelvin = dataArray{:, 4};
SystemTempKelvin = dataArray{:, 5};
FieldOersted = dataArray{:, 6};
PressureTorr = dataArray{:, 7};
SampleTempKelvin = dataArray{:, 8};
TempRiseKelvin = dataArray{:, 9};
SampHCmJmoleK = dataArray{:, 10};
SampHCErrmJmoleK = dataArray{:, 11};
AddendaHCJK = dataArray{:, 12};
AddendaHCErrJK = dataArray{:, 13};
TotalHCJK = dataArray{:, 14};
TotalHCErrJK = dataArray{:, 15};
FitDeviationChiSquare = dataArray{:, 16};
TimeConsttau1seconds = dataArray{:, 17};
TimeConsttau2seconds = dataArray{:, 18};
SampleCouplingPercent = dataArray{:, 19};
DebyeTempKelvin = dataArray{:, 20};
DebyeTempErrKelvin = dataArray{:, 21};
CalCorrectionFactor = dataArray{:, 22};
ThermResistOhms = cellstr(dataArray{:, 23});
HtrResistOhms = cellstr(dataArray{:, 24});
PuckResistOhms = cellstr(dataArray{:, 25});
WireCondWK = dataArray{:, 26};
MeasTimeseconds = dataArray{:, 27};
TempSquaredK2 = dataArray{:, 28};
SampHCTempmJmoleKK = dataArray{:, 29};
AddendaOffsetHCJK = dataArray{:, 30};

