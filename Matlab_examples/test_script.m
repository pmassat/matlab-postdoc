classdef test_script
   properties
      Value
   end
   methods
      function r = roundOff(obj)
         r = round([obj.Value],2);
      end
      function r = multiplyBy(obj,n)
         r = [obj.Value] * n;
      end
      function r = multRound(obj,A)
          r = ftest(roundOff(obj),A);
      end
   end
end

function y = ftest(x,A)
    y = A*x;
end