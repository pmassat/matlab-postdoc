function F = freeEnergy_TLFIM_compute(x,h,e)
%% Plot
% F = zeros(length(x),length(h),length(e));
% for jh = 1:length(h)
%     for je = 1:length(e)
%         for j=1:length(x)
%             F(j,jh,je) = freeEnergy_TLFIM(x(j),h(jh),e(je));% free energy
%         end
%     end
% end

F = zeros(length(x),length(h));
for jh = 1:length(h)
    for j=1:length(x)
        F(j,jh) = freeEnergy_TLFIM(x(j),h(jh),e(jh));% free energy
    end
end
end
