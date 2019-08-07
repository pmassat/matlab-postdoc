function ytot = IkedaCarpenterPV(t,params)
    alpha = params(3);
    beta = params(4);
    gamma = params(1);
    I = params(6);
    R = params(5);
    sigmaSquared = params(2);
    t0 = params(7);
%     alpha = 1/3;
%     beta = 1/31.9;
%     gamma = 1;
%     I = 1;
%     R = 0.18;
%     sigmaSquared = 1/(8*log(2));
%     t0 = 1;
%     syms t
% 
%    manipulate({@(x,param)(param(1)*x+param(2)),[-10:10],[2 0]},{[-10 10],[-10 10]}, {1,'Slope',1 10}, {2,'Shift',-5 5})
%    manipulate({IkedaCarpenter),[-10:10],[2 0]},{[-10 10],[-10 10]}, {1,'Slope',1 10}, {2,'Shift',-5 5})

    k = 0.05;
    diff = t - t0;
%      diff = xValues[i] - X0;
%      R = exp(-81.799 / (m_waveLength[i] * m_waveLength[i] * kappa));
%      alpha = 1.0 / (alpha0 + m_waveLength[i] * alpha1);

    [W,eta] = convertVoigtToPseudo(sigmaSquared,gamma);
    %Width and weight factor of pseudo-Voigt

    a_minus = alpha * (1 - k);
    a_plus = alpha * (1 + k);
    x = a_minus - beta;
    y = alpha - beta;
    z = a_plus - beta;
    
    Nu = 1 - R * a_minus / x;
    Nv = 1 - R * a_plus / z;
    Ns = -2 * (1 - R * alpha / y);
    Nr = 2 * R * alpha * alpha * beta * k * k / (x * y * z);
    
    u = a_minus * (a_minus * sigmaSquared - 2 * diff) / 2.0;
    v = a_plus * (a_plus * sigmaSquared - 2 * diff) / 2.0;
    s = alpha * (alpha * sigmaSquared - 2 * diff) / 2.0;
    r = beta * (beta * sigmaSquared - 2 * diff) / 2.0;
 
    someConst = 1 / sqrt(2.0 * sigmaSquared);% sigmaSquared should be strictly positive
    
    yu = (a_minus * sigmaSquared - diff) * someConst;
    yv = (a_plus * sigmaSquared - diff) * someConst;
    ys = (alpha * sigmaSquared - diff) * someConst;
    yr = (beta * sigmaSquared - diff) * someConst;

    zs = (-alpha * diff + 1i* 0.5 * alpha * gamma);
    zu = (1 - k) * zs; 
    zv = (1 + k) * zs;
    zr = (-beta * diff + 1i* 0.5 * beta * gamma);

    N = 0.25 * alpha * (1 - k * k) / (k * k);
        
    yg = (Nu * exp(u).*erfc(yu) +...
           Nv * exp(v).*erfc(yv) +...
           Ns * exp(s).*erfc(ys) +...
           Nr * exp(r).*erfc(yr));
    
    yl = - 2.0 / pi *...
              (Nu * imag(exp(zu).*expint(zu)) +...
               Nv * imag(exp(zv).*expint(zv)) +...
               Ns * imag(exp(zs).*expint(zs)) +...
               Nr * imag(exp(zr).*expint(zr)));

    ytot = I * N *...
         ((1 - eta) * yg +...
          eta * yl);
%     %%
%     if select<2
%         y = ytot;
%     elseif select>2
%         y = yl;
%     else
%         y = yg;
%     end
end 


