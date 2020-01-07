function y = CpTFIM_normpdf(t,hm,sgm)% Normal PDF of Schottky anomaly as a function of temperature T
    y = integral(@(h)Cp_TFIM(t,h).*normpdf(h,hm,sgm),0,inf,'ArrayValued',true);
end
