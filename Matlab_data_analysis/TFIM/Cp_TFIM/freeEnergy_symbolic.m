syms Gam beta gam h x
%% This code is unnecessary since the numerical one works
% J = 1;
% Equations are taken from Stinchcombe 1973(III) and adapted so as to have
% kB*T_c = 1 and Gam/Gam_c = 1;
% Equation labels from the paper are indicated as comments
Eq = gam == tanh(beta*gam);% Equation (2.2) in the paper
% Eq1 = sin(theta) == Gam/gam;% Equation (2.7)
% Eq2 = gam^2 == Gam^2 + (h + J*R*cos(theta))^2;% Equation (2.6)
% solve([Eq1,Eq2],[gam,theta])% Warning: Unable to find explicit solution.
% gam(Gam,h,x) = sqrt(Gam^2 + (h + x)^2);% Equation (2.6)
% Eq = x == (x+h)/gam(Gam,h,x)*tanh(beta*gam(Gam,h,x));
% solve([Eq],[x]);% Warning: Unable to find explicit solution.
% x(Gam,gam,h) = sqrt(tanh(beta*gam)^2-Gam^2)-h;
% bf(beta,Gam,h,gam) = 0.5*beta*x(Gam,gam,h)^2 - log(2*cosh(beta*gam));% Equation (2.5)
bf(beta,Gam,h,gam) = - log(2*cosh(beta*gam));% Equation (2.5)
d1bf = diff(bf,beta);
d2bf = diff(d1bf,beta);
Cp(beta,Gam,h,gam) = -beta^2*d2bf;
