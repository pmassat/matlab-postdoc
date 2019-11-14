function y = op_fit(t,h,e)
y = zeros(size(t));
for i=1:length(t)
    if t(i)>=0
    y(i) = OP_TFIM(t(i),h,e);
    end
end
end