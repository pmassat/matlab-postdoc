 function y = Cp_TLFIM(t,h,sigma)% This function was copied from 'Cp_TFIM_offset_strain.m' on 2019-12-04
% t is the reduced temperature T/Tc
y = zeros(size(t));

if h==0; tc = 1;
elseif abs(h)<1; tc = h/atanh(h);
else; tc = 0;
end
% tc =1;

for i=1:length(t)
    y(i) = 0;
    x = OP_TFIM(t(i),h,sigma);
    gamma = sqrt((x+sigma)^2+h^2);
    r = gamma./t(i);
    D = 1 - sech(r)^2./t(i);
    if t(i)<tc% || e~=0
        y(i) = r^2.*sech(r)^2 ./ D;% mean-field heat capacity in the ordered phase
%         r1 = order_parameter_offset_strain(tr,e)./tr;
%         % ratio of reduced order_parameter to reduced temperature
% % Note: the order parameter depends on the ratio t(i)/tc, not on t(i)
%         y(i) = r1^2.*sech(r1)^2 ./ (1 - 1./tr.*sech(r1)^2);% mean-field heat capacity in the ordered phase
    else
%         r11 =0;
%         if abs(e)>0
%             r11 = order_parameter_offset_strain(t(i),e)./t(i);
%         else; r11=0;
%         end
%         if e==0
%         y(i) = r2^2.*sech(r2)^2;% mean-field heat capacity in the disordered phase
%         else
%         r3 = (order_parameter_offset_strain(t(i),e))./t(i);% ratio of reduced order_parameter to reduced temperature
%         Cord = r3^2.*sech(r3+e)^2 ./ (1 - 1./t(i).*sech(r3+e)^2);% mean-field heat capacity in the ordered phase
%         y(i) = r^2.*gc.*sech(r)^2 ./ (D*gc + e.*h^2./gamma);% mean-field heat capacity in the ordered phase
        y(i) = r^2.*sech(r)^2;% mean-field heat capacity in the disordered phase            
%         end
    end
end
end