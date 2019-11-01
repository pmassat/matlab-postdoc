function y = Cp_nem_fluct(c1,c2,t,T1,T2)
y = zeros(size(t));
for i=1:length(t)
    if t(i)>1
%     y(i) = c1/(t(i)-T1) + c2*7.1/((t(i)-T1)^2*(t(i)-T2)) + c2*(t(i)+1.7)/((t(i)-T1)*(t(i)-T2)^2);
    y(i) = t(i)/(t(i)-T1) * ( -7.1*c1/(t(i)-T1)^2 + c2*7.1/((t(i)-T1)*(t(i)-T2)) + c2*(t(i)+1.7)/(t(i)-T2)^2 );
%     elseif t(i)<1
%     y(i) = exp((t(i)-1)/T0);
    end
end
end