function M = magnetization_TFIM(t,h)
    M = zeros(size(h));
    hc = critical_field(t);
    for i=1:length(h)
        if h(i)<hc
            M(i) = h(i);
        else
            M(i) = tanh(h(i)/t);
        end
    end
end