%% This code is for modeling DyVO4, a Jahn-Teller system


%% Below is the CEF Hamiltonian generation and definitions of Stevens operators

J=15/2;

%Defining Jplus, Jminus and Jz in the mJ representation
Jplus=zeros(2*J+1);
Jminus=zeros(2*J+1);
for k=1:(2*J+1)
    for i=1:(2*J+1)
        if i==k
            Jz(i,k)=i-(J+1);
        end
        if i==k+1
            Jplus(i,k)=sqrt((J-(k-(J+1)))*(J+k-(J+1)+1));
        end
        if i==k-1
            Jminus(i,k)=sqrt((J+(k-(J+1)))*(J-(k-(J+1))+1));
        end 
    end
end

JT=Jz^2+Jplus*Jminus-Jz;
Jx=(Jplus+Jminus)/2;
Jy=(Jplus-Jminus)/(2*sqrt(-1));

%Making sure numbers are input as diagonal matrices
D=diag(ones((2*J+1),1));
X=diag(J*(J+1)*ones((2*J+1),1));

%Saving Jz^4 and Jz^6 as matrices to save space
Jz4=Jz*Jz*Jz*Jz;
Jz6=Jz*Jz*Jz*Jz*Jz*Jz;
Jpm4=Jplus*Jplus*Jplus*Jplus+Jminus*Jminus*Jminus*Jminus;
Jpm2=Jplus*Jplus+Jminus*Jminus;

%Defining Stevens operators
Pxy=-0.25*sqrt(-1)*(Jplus*Jplus-Jminus*Jminus);
O22=0.5*(Jplus*Jplus+Jminus*Jminus);
O02=3*Jz*Jz-X;
O04=35*Jz4-(30*X-25*D)*Jz*Jz+3*X^2-6*X;
O44=0.5*(Jplus*Jplus*Jplus*Jplus+Jminus*Jminus*Jminus*Jminus);
O06=231*Jz6-(315*X-735*D)*Jz4+(105*X*X-525*X+294*D)*Jz*Jz-5*X^3+40*X^2-60*X;
O46=0.25*(Jpm4*(11*Jz*Jz-X-38*D)+(11*Jz*Jz-X-38*D)*Jpm4);
O24=0.25*(Jpm2*(7*Jz*Jz-X-5*D)+(7*Jz*Jz-X-5*D)*Jpm2);

%For TmAg 

% see https://obelix.physik.uni-bielefeld.de/~schnack/molmag/material/Lueken-kurslan_report.pdf
% for alphaJ, betaJ, gammaJ (table 23)

% see https://www.researchgate.net/publication/243535033_Anomalies_in_thermal_expansion_of_DyVO_4_induced_by_quadrupolar_ordering
% for V02 etc.
V02=-114;
V04=50;
V44=973;
V06=-59;
V46=182;
alphaj=-2/315;
betaj=-8/135135;
gammaj=4/3864861;
B02=alphaj*V02;
B04=betaj*V04;
B44=betaj*V44;
B06=gammaj*V06;
B46=gammaj*V46;
 
%Crystal field Hamiltonian for tetragonal symmetry with material dependent
%coefficients 
Hcef=B02*O02+B04*O04+B44*O44+B06*O06+B46*O46;


[V,D]=eig(Hcef);
[sorteig,P]=sort(eig(Hcef),'descend');
realeig=sorteig-min(sorteig);
sortstates=V(:,P);
gsdoublet0=sortstates(:,2*J:2*J+1);

T = 1;
expect(O22,Hcef,T);

%not sure what this does
% for i=1:2
%     for j=1:2
%     OM(i,j)=gsdoublet0(:,j)'*O22*gsdoublet0(:,i);
%     PM(i,j)=gsdoublet0(:,j)'*Pxy*gsdoublet0(:,i);
%     end
% end
% OM;
% PM;

Kg = 10^-3;
Bg=32;
Eg=.000;

%Morin
Bg=32;
Gg=18e-3;
Kg=10e-3;
Cg=14.2e4; %Morin
Gxy=10e-3; %high limit from Morin
Bxy=100;
EA1g=0;

%% Below is the rest of the Hamiltonian

Hz = 0;
%    Exy1=Exy(j);
%    j
%    for i=1:L
    T1= 1:.2:20 ;%=T1(i);
    %delta is the convergence allowance
    delta=0.000001;
    l=2;

 %   M(1)=0;
 %   M(2)=0.001;
 %   RelaxstrainQ(i,j)=0;
 expectation_values = zeros(size(T1));
  expectation_final = zeros(size(T1));
 jj = 1;
for Ti = T1
    clear Qf
    jj = jj+1;
        Qf(1)=0;
    Qf(2)=0.001;
    for l=2:10%while abs(Qf(l)-Qf(l-1))>delta && l<100
            %Hzeeman=-gz*ub*Hz*Jz
        Qf;
            
            strain=(Bg/(2*Cg))*Qf(l) + Eg;
            Hme=-Bg*(strain)*O22;
            Hel = Cg * strain^2*eye(2*J+1);
            %Quadrupolar coupling only for O22
            %Hqq=-Kg*Qf(l)*O22; 
            
            %Hgtotal=-Gg*Qf(l)*O22;
            
            
            %Magneto-elastic coupling only for O2_2
            %Hme2=-Bxy*sqrt(2)*Exy1*Pxy;
            
            %HmeA=-B02*EA1g*O02;
            %Mean field correction to energy
            Ecorr=-(0.5*Kg*Qf(l)^2)*eye(2*J+1);%+0.5*n1*(expect(JT,Hcef,T))^2;%+0.5*n2*(2*Mz^2-Mx^2-My^2)+0.5*Ka*(expect(O02,T))^2;
            
            %Total Hamiltonian
            Htotal=Hcef+Hme+Ecorr+Hel;%+Hme+Hqq+Ecorr
        
            
            %abs(Qf(l)-Qf(l-1))
            Qf(l+1)=expect(O22,Htotal,Ti);
            
           % Qfa(l)=expect(O02,Htotal,T);
          %  Qt(l)=expect(Pxy,Htotal,T);
          Ti;
          
    end
    %Htotal=Hcef+Hme+Ecorr;
    %Hfinal{i}=Hcef+0.5*Cg*Sponstrain(i)^2*eye(2*J+1)-%-Bxy*sqrt(2)*Exy1*Pxy-B02*EA1g*O02; %Kg*RelaxstrainQ(i)^2/2
    %+0.5*(Bg^2/14.2e4)*(Qf(l)^2)*eye(2*J+1)
%finalexpectationvalues= expect(O22,Htotal,Ti);
expectation_values(jj) = Qf(end);
Hfinal = Hcef-0.25*(Bg^2/Cg)*Qf(end)^2*eye(2*J+1)-Bg*(Eg)*O22;
expectation_final(jj) = expect(O22,Hfinal,Ti);
end
expectation_values(:,1) = [];
expectation_final(:,1) = [];
%T1 = [0 T1];
expectation_values = real(expectation_values);
plot(T1,expectation_values,'DisplayName','original expectation values')
% hold on
% plot(T1,expectation_final,'DisplayName','final expectation values')
% hold off
legend
