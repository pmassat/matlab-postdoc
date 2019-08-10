function w=Wpq(J,vec,p,q)

w=0;

MM=-J:1:J;

for i=1:length(MM)
    for j=1:length(MM)
        wig = Wigner3j( [J,J,p], [-MM(j) MM(i) q] );
        w=w+(-1)^(J+MM(i)-p)*wig*vec(j)'*vec(i);
    end
end