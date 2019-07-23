syms Gam beta gam h theta
J = 1;
% Equations are taken from Stinchcombe 1973(III) and adapted so as to have
% kB*T_c = 1 and Gam/Gam_c = 1;
% Equation labels from the paper are indicated as comments
R(beta,gam) = tanh(beta*gam);% Equation (2.2) in the paper
Eq1 = sin(theta) == Gam/gam;% Equation (2.7)
Eq2 = gam^2 == Gam^2 + (h + J*R*cos(theta))^2;% Equation (2.6)
% solve([Eq1,Eq2],[gam,theta])% Warning: Unable to find explicit solution.
bf(beta,Gam,h) = 2*beta*J*R^2*cos(theta)^2 - log(2)*cosh(beta*gam);% Equation (2.5)
d1bf = diff(bf,beta);
d2bf = diff(d1bf,beta);
Cp(beta, Gam, h) = -beta^2*d2bf;
% Does not work for numerical values. Try to combine with computation of
% OP...