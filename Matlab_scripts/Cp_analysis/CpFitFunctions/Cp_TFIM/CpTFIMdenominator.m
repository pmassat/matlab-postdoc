function y = CpTFIMdenominator(t,e)
    r = order_parameter_offset_strain(t,e)./t;
    y = (1 - sech(r)^2./t);
end