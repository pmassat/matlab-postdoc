 function y = Cp_TFIM_offset_strain(t,e,h)
% t is the reduced temperature T/Tc
y = zeros(size(t));

if h==0; tc = 1;
elseif abs(h)<1; tc = h/atanh(h);
else; tc = 0;
end
% tc =1;

for i=1:length(t)
    y(i) = 0;
    tr = t(i)/tc;
    if t(i)<tc || (h==0 && abs(e)>0)
        r1 = order_parameter_offset_strain(t(i),e)./t(i);
        % ratio of reduced order_parameter to reduced temperature
% Note: the order parameter depends on the ratio t(i)/tc, not on t(i)
        y(i) = r1^2.*sech(r1)^2 ./ (1 - 1./t(i).*sech(r1)^2);% mean-field heat capacity in the ordered phase
%         r1 = order_parameter_offset_strain(tr,e)./tr;
%         % ratio of reduced order_parameter to reduced temperature
% % Note: the order parameter depends on the ratio t(i)/tc, not on t(i)
%         y(i) = r1^2.*sech(r1)^2 ./ (1 - 1./tr.*sech(r1)^2);% mean-field heat capacity in the ordered phase
%     else
%         r11 =0;
% %         if abs(e)>0
% %             r11 = order_parameter_offset_strain(t(i),e)./t(i);
% %         else; r11=0;
% %         end
%         r2 = r11 + h./t(i);
% %         if e==0
% %         y(i) = r2^2.*sech(r2)^2;% mean-field heat capacity in the disordered phase
% %         else
% %         r3 = (order_parameter_offset_strain(t(i),e))./t(i);% ratio of reduced order_parameter to reduced temperature
% %         Cord = r3^2.*sech(r3+e)^2 ./ (1 - 1./t(i).*sech(r3+e)^2);% mean-field heat capacity in the ordered phase
%         y(i) = y(i) + r2^2.*sech(r2)^2;% mean-field heat capacity in the disordered phase            
% %         end
    end
end
end