function y = Cp_nem_fluct(c1,c2,T,T1,T2)
y = zeros(size(T));
for i=1:length(T)
    y(i) = c1/(T(i)-T1) + c2*7.1/((T(i)-T1)^2*(T(i)-T2)) + c2*(T(i)+1.7)/((T(i)-T1)*(T(i)-T2)^2);
end
end