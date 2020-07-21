function hc = critical_field(t)
    if t<1
        hc = fzero(@(x) x-tanh(x/t),[1e-3 1]);
    else
        hc = 0;
    end
end