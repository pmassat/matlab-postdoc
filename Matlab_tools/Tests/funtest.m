function y = funtest(arg1, varargin)

p = inputParser;

addOptional(p, 'arg2', 0);

parse(p, varargin{:});

y = arg1 + p.Results.arg2;
end

% function y = funtest(x,param)
%     a = param(1);b = param(2);
% %     a=3;b=2;
%     y = a.*(x-b).^2;
% end
