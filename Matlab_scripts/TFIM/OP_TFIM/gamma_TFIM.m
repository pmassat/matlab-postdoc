function y = gamma_TFIM(t,h,e)
% See notations of Stinchcombe 1973
x = OP_TFIM(t,h,e);
y = sqrt((x+e)^2+h^2);
end