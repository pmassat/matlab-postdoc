function y = Cp_NF_expr(c1,c2,T,T1,T2)
y = c1/(T-T1) + c2*7.1/((T-T1)^2*(T-T2)) + c2*(T+1.7)/((T-T1)*(T-T2)^2);
end