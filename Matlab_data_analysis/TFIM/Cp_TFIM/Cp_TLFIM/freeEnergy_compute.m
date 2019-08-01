function F = freeEnergy_compute(x,h,e)
%% Plot
F = zeros(length(x),length(h),length(e));
for jh = 1:length(h)
    for je = 1:length(e)
        for j=1:length(x)
            F(j,jh,je) = freeEnergy_TLFIM(x(j),h(jh),e(je));% free energy
        end
    end
end
