function rho=densityMagneticAll(J,vec,tth,pph,label,l)

if strcmp('Ce3+',label)
    MAXI=5;
else
    MAXI=6;
end

rho=zeros(size(tth));

gn=GNions(label);

if nargin==6
    range=l;
else
   range=1:2:MAXI;
end

for p=range
    factor=(2*p+1)*gn(p+1)*RedMatEle(J,p);
    for q=-p:p                
        rho=rho+factor*(Wpq(J,vec,p,q)*Zpq(p,q,tth,pph));
    end
end
    