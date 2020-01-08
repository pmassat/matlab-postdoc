%% 
H = 0:0.1:10;
f = zeros(size(H));
t = .1;
e = 0;
for k=1:length(H)
    f(k) = freeEnergy_TLFIM(t,H(k),e).*normpdf(H(k),h,rhsgm*h);
end

%%
figure; plot(H,f);

%%
freeEnergy_TLFIM_normpdf(t,h,e,rhsgm*h)

%%
func = @(x)freeEnergy_TLFIM(t,x,e).*normpdf(x,h,rhsgm*h);
figure;
fplot(func, [0 10]);

%%
I = integral(func,0,10,'ArrayValued',true);