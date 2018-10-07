scatter3(Magneticfield,PuckTemperature,Thermometerresistance,'k.')
%%
beta = 1./PuckTemperature;
% General model:
%      f(x,y) = a + c/y^(2) +c5*x/y+ p10*x + p01*y + p20*x^2 + p11*x*y+ d *y^5 
%                     + p02*y^8 + p30*x^-3 + p21*x^3*y+ p3333*x*y^2
% Coefficients (with 95% confidence bounds):
%        a =         911  (893.8, 928.2)
%        c =       777.7  (775, 780.3)
%        c5 =     -0.1172  (-0.12, -0.1145)
%        d =      0.9592  (0.6972, 1.221)
%        p01 =      -142.9  (-154.7, -131.1)
%        p02 =    -0.03734  (-0.04425, -0.03044)
%        p10 =      0.2268  (0.2175, 0.2361)
%        p11 =     -0.1491  (-0.1562, -0.142)
%        p20 =   2.335e-06  (1.88e-06, 2.79e-06)
%        p21 =   -1.17e-10  (-1.405e-10, -9.353e-11)
%        p30 =   0.0007512  (0.0004269, 0.001075)
%        p3333 =     0.03014  (0.02847, 0.0318)
% 
% Goodness of fit:
%   SSE: 5.558e+05
%   R-square: 0.9998
%   Adjusted R-square: 0.9998
%   RMSE: 20.47
%%
 scatter3(Magneticfield,PuckTemperature,(TModel(Magneticfield,PuckTemperature)-Thermometerresistance)./Thermometerresistance)
 
 %%
 
 scatter3(Magneticfield,PuckTemperature,(1./BetaModel(Magneticfield,Thermometerresistance)-PuckTemperature)./PuckTemperature)
 
 %%
 
PPMSMCE(1)=ImportPPMSMCE('LogPPMSMCE_3.dat');
PPMSMCE(2)=ImportPPMSMCE('LogPPMSMCE_4_15-OePerSec.dat');
PPMSMCE(3)=ImportPPMSMCE('LogPPMSMCE_5_80-OePerSec.dat');
PPMSMCE(4)=ImportPPMSMCE('LogPPMSMCE_5_5-OePerSec.dat');
%%
for j = 1:4
    plot(PPMSMCE(j).H,PPMSMCE(j).RPlatform)
    hold on
end

%%
for j = 1:4
    PPMSMCE(j).TPlat=1./BetaModel(PPMSMCE(j).H,PPMSMCE(j).RPlatform);
    plot(PPMSMCE(j).H,PPMSMCE(j).TPlat)
    hold on
end
xlabel('Field (Oe)')
ylabel('Temperature (K)')
%%
PPMSPartitions(2).Up =[621 2707; 7651 10071;35000 38000];
PPMSPartitions(2).Down =  [3712 6130; 10900 12934; 38400 41500];
for j = 2
    for k = 1:length(PPMSPartitions(2).Up)
        plot(PPMSMCE(j).H(PPMSPartitions(2).Up(k,1):PPMSPartitions(2).Up(k,2)),PPMSMCE(j).TPlat(PPMSPartitions(2).Up(k,1):PPMSPartitions(2).Up(k,2)))
        hold on
    end
        
    
end
ylabel('Field (Oe)')