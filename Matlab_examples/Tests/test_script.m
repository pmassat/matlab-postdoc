% classdef test_script %< handle
%     properties
%         Value = 1;
%         val2;
%     end
%     methods
%         function obj = test_script
%             obj.val2 = obj.Value;
%         end
%         function y = computeValue(obj, val)
%             obj.val2 = val;
%             y = obj.Value + obj.val2;
%         end
%         function r = roundOff(obj)
%          r = round([obj.Value],2);
%         end
%         function r = multiplyBy(obj,n)
%          r = [obj.Value] * n;
%         end
%         function r = multRound(obj,A)
%           r = ftest(roundOff(obj),A);
%         end
%     end
% end
% 
 strcat("this is a"+newline+"sample text ",...
     "used for a test")