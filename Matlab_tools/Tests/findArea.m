function a = findArea(width,varargin)
% This is the example given at 
% See https://www.mathworks.com/help/matlab/ref/inputparser.html
% to show how the inputParser class works, in order to create
% custom arguments in functions
   defaultHeight = 1;
   defaultUnits = 'inches';
   defaultShape = 'rectangle';
   expectedShapes = {'square','rectangle','parallelogram'};

   p = inputParser;
   validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
   addRequired(p,'width',validScalarPosNum);
   addOptional(p,'height',defaultHeight,validScalarPosNum);
   addParameter(p,'units',defaultUnits,@isstring);
   addParameter(p,'shape',defaultShape,...
                 @(x) any(validatestring(x,expectedShapes)));
   parse(p,width,varargin{:});
   
   a = p.Results.width*p.Results.height; 
end