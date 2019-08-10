function pol=Pleg(l,m,X)

val=legendre(l,X);

if l==0
    pol=val;
elseif m<0
    pol=(-1)^m*factorial(l+m)/factorial(l-m)*permute(val(abs(m)+1,:,:),[2,3,1]);
elseif m>=0
    pol=permute(val(m+1,:,:),[2,3,1]);
end
