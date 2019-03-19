function y = order_parameter(t)
%%% calculate the order parameter of the MF FQ transition (which equals the
%%% normalized pseudospin) as a function of reduced temperature t = T/Tc
    y = fzero(@(x) x-tanh(x/t), [1e-3 1]);
    % if 0 is included in the range of fzero, y will always be zero, since
    % tanh(0) = 0
end