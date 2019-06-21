function s = OP_TFIM(t,h,e)
    if h==0; tc = 1;
    elseif abs(h)<1; tc = h/atanh(h);
    else; tc = 0;
    end
    r = fzero(@(x) x-e-tanh((x)./(t)), 1);
    s = sqrt(r^2-(h)^2);
end