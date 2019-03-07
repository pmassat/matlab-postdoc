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

% generate data
x = 0:.1:10;
y = x.*x + randn(size(x));
w = linspace(0.1, 10.1,length(x));
x = x(:);
y = y(:);
w = w(:);
%plot data
figure
errorbar(x,y,1./w,'.');
%fit
ft = fittype('poly2');
cf = fit(x,y,ft,'Weight',w);
% Plot fit
hold on
plot(cf,'fit',0.95);
