function val = voigtIkedaCarpenter_ord(t,params)
a = params(2);
b = params(3);
R = params(1);
gamma = params(4);
sigma = params(5);
k = params(6);
t0 = params(7);
xi = t - t0;
%VoigtIkedaCarpenter Calculates the convolution of a pseudo-Voigt
%profile with the Ikeda-Carpenter function
%   xi - position
%   a - I-C parameter
%   b - I-C parameter
%   R - I-C parameter
%   gamma - Voigt parameter (Lorentzian part)
%   sigma - Voigt parameter (Gaussian part)
%   k - I-C approximation parameter

    function Gamma = pseudoVoigtFWHM(gamma, sigma)
        fL = 2*gamma; 
        fG = 2*sigma*sqrt(2*log(2));
        Gamma = (fG^5 + 2.69269*fG^4*fL + 2.42843*fG^3*fL^2 + ...
            4.47163*fG^2*fL^3 + 0.07842*fG*fL^4 + fL^5)^(1/5);
    end
        
    function eta = pseudoVoigtEta(gamma, sigma)
        fL = 2*gamma;
        f = pseudoVoigtFWHM(gamma, sigma);
        fLbyf = fL/f;
        eta = 1.36603*fLbyf - 0.47719*fLbyf^2 + 0.11116*fLbyf^3;
    end

Gamma = pseudoVoigtFWHM(gamma, sigma);
eta = pseudoVoigtEta(gamma, sigma);

sigmaSq = Gamma^2/(8*log(2));
gWidth = sqrt(2*sigmaSq);

am = a*(1 - k);
ap = a*(1 + k);

x = am - b;
y = a - b;
z = ap - b;

zs = -a*xi + 1i*a*Gamma/2;
zu = (1 - k)*zs;
zv = (1 + k)*zs;
zr = -b*xi + 1i*b*Gamma/2;

u = am*(am*sigmaSq - 2*xi)/2;
v = ap*(ap*sigmaSq - 2*xi)/2;
s = a*(a*sigmaSq - 2*xi)/2;
r = b*(b*sigmaSq - 2*xi)/2;

n = (1/4)*a*(1 - k^2)/k^2;
nu = 1 - R*am/x;
nv = 1 - R*ap/z;
ns = -2*(1 - R*a/y);
nr = 2*R*a^2*b*k^2/(x*y*z);

yu = (am*sigmaSq - xi)/gWidth;
yv = (ap*sigmaSq - xi)/gWidth;
ys = (a*sigmaSq - xi)/gWidth;
yr = (b*sigmaSq - xi)/gWidth;

val = n*((1 - eta)*(nu*exp(u).*erfc(yu) + nv*exp(v).*erfc(yv)...
    + ns*exp(s).*erfc(ys) + nr*exp(r).*erfc(yr))...
    - eta*2/pi*(imag(nu*expint(zu).*exp(zu)) + imag(nv*expint(zv).*exp(zv))...
    + imag(ns*expint(zs).*exp(zs)) + imag(nr*expint(zr).*exp(zr))));

end

