function y = dOPdT_randFields(t,sz)
% Numerical derivation of the order parameter in random fields
% Yields the same result as the analytical computation performed by 'order_parameter_derivative_random_strains_offset.m'
% t is the reduced temperature T/Tc
T = t;% transpose temperature array into column vector

dT = diff(T);
diff1sz = diff(sz,1,2);% first order difference along 1st dimension (i.e. between rows) of F
[m,~] = size(sz);% in case F is a multi-column array, e.g. calculated for various values of a parameter
y = zeros(size(diff1sz));
for j = 1:m% loop over all columns of F
    y(j,:) = diff1sz(j,:)./dT;% dsz/dT
end

end
