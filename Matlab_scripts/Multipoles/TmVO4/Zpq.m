function Zlm=Zpq(p,q,tth,pph)

Zlm=sqrt((factorial(p-q))/(factorial(p+q))) * ...
        Pleg(p,q,cos(pph)).*exp(1i*q*tth);
