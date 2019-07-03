function y = order_parameter_offset_strain(t,e)
%%% calculate the order parameter of the MF FQ transition (which equals the
%%% normalized pseudospin) as a function of reduced temperature t = T/Tc
y = fzero(@(x) x-tanh((x+e)/t), 1);
% Use 1 as starting point for fzero
% If starting from 0, y will always be zero, since tanh(0) = 0
end