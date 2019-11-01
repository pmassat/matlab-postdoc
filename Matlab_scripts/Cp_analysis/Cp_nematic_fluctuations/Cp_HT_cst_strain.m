function y = Cp_HT_cst_strain(c1,c2,t,T1,T2)
y = zeros(size(t));
for i=1:length(t)
    if t(i)>1
    y(i) = t(i) * ( 7.1*c1/(t(i)-T1)^3 + c2/(t(i)-T2)^2 );
    end
end
end