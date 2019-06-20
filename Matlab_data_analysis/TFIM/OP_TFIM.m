function y = OP_TFIM(t,h)
    if h==0; tc = 1;
    elseif abs(h)<1; tc = h/atanh(h);
    else; tc = 0;
    end
    r = fzero(@(x) x-tanh(x./(t/tc)), 1);
    y = r*sqrt(1-h^2);
end