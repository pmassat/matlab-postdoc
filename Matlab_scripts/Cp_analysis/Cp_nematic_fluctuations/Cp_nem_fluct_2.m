function y = Cp_nem_fluct_2(c1,c2,t,a,T1,T2)
y = zeros(size(t));
for i=1:length(t)
    if t(i)>1
    y(i) = c1/(t(i)-T1) + c2*7.1/((t(i)-T1)^2*(t(i)-T2)) + c2*(t(i)+1.7)/((t(i)-T1)*(t(i)-T2)^2);
    elseif t(i)<1
    y(i) = t(i)^a;
    end
end
end