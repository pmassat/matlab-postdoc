
T = linspace(0,2,201);
h = 0;
e = 0.01; 
op = repmat(T,1);
noise = (rand(201,1)-0.5)*0.1;
for i=1:length(T)
    op(i) = OP_TFIM(T(i),h,e)+noise(i);
end

%%
plot(T,op);