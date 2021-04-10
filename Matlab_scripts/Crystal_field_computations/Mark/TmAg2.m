function[]=TmAg2(T,strains,coeff)
%[RelaxstrainQ,T2,Exy]
%Tq=Quadrupolar T
%Qmatrix= Relaxed spontaneous Q at each Temp
%T1 is temperature vector
%MaxstrainQ= the max value of Q
%gsenergy are the groundstate energies of doublet and singlet at each temp
%gs... are elements of those quadrupoles at each temp among the groundstate
%doublets and singlet
%%
J=6;

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

%Morin
Bg=32;
Gg=18e-3;
Kg=10e-3;
Cg=14.2e4; %Morin
Gxy=10e-3; %high limit from Morin
Bxy=100;
EA1g=0;

%For TmAg
V02=9;
V04=2;
V44=-395;
V06=-16;
V46=230;
alphaj=1/99;
betaj=8/(81*5*11*11);
gammaj=-5/(81*7*11*11*13);
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

for i=1:2
    for j=1:2
    OM(i,j)=gsdoublet0(:,j)'*O22*gsdoublet0(:,i);
    PM(i,j)=gsdoublet0(:,j)'*Pxy*gsdoublet0(:,i);
    end
end
OM
PM
%%
%Temperatures


%Strain values
% Energyrelaxed=zeros(2*J+1,L);
% Sponstrain=zeros(1,L);
% 
% 
% tic 
%  %For loop runs through each temperature
% 
%  
% for j=1:length(strains) 
%     Exy1=Exy(j);
%     j
%     for i=1:L
%     T=T1(i);
%     %delta is the convergence allowance
%     delta=0.000001;
%     l=2;
%     Qf(1)=0;
%     Qf(2)=0.001;
%     M(1)=0;
%     M(2)=0.001;
%     RelaxstrainQ(i,j)=0;
%     
%     %While loop guesses Hamiltonian until the guess converges
%     while abs(Qf(l)-Qf(l-1))>delta && l<100
%             Hzeeman=-gz*ub*Hz*Jz
%         
%             Hme=-Bg*(Eg)*O22;
%             
%             %Quadrupolar coupling only for O22
%             %Hqq=-Kg*Qf(l)*O22; 
%             
%             Hgtotal=-Gg*Qf(l)*O22;
%             
%             
%             %Magneto-elastic coupling only for O2_2
%             %Hme2=-Bxy*sqrt(2)*Exy1*Pxy;
%             
%             %HmeA=-B02*EA1g*O02;
%             %Mean field correction to energy
%             %Ecorr=(0.5*Kg*Qf(l)^2)*eye(2*J+1);%+0.5*n1*(expect(JT,Hcef,T))^2;%+0.5*n2*(2*Mz^2-Mx^2-My^2)+0.5*Ka*(expect(O02,T))^2;
%             
%             %Total Hamiltonian
%             Htotal=Hcef+Hgtotal+Hme2;%+Hme+Hqq+Ecorr
%         
%             l=l+1;
%             
%             Qf(l)=expect(O22,Htotal,T);
%             
%            % Qfa(l)=expect(O02,Htotal,T);
%           %  Qt(l)=expect(Pxy,Htotal,T);
%     end
% %For each temp i calculates the spontaneous Q and strain and energy spectrum
%     %RelaxstrainQ(i)=Qf(l);
%    % FixedstrainQa=Qfa(l);
%     Sponstrain(i)=-Bg/Cg*Qf(l);
%     Hfinal{i}=Hcef+0.5*Cg*Sponstrain(i)^2*eye(2*J+1)-%-Bxy*sqrt(2)*Exy1*Pxy-B02*EA1g*O02; %Kg*RelaxstrainQ(i)^2/2
%     
%     %Hfinalrelaxed=Hcef;
% %     Energyrelaxed(:,i)=sort(eig(Hcef));
% %     
% %     [V,D]=eig(Hfinalrelaxed);
% %     [sorteig,P]=sort(real(eig(Hfinalrelaxed)),'descend');
% %     realeig=sorteig-min(sorteig);
% %     gsenergy(i,:)=sorteig(1:2*J+1);
% %     sortstates=V(:,P);
% %    % gstriplet=sortstates(:,2*J-1:2*J+1);
% %     gscell{i}=sortstates;
% %     for p=1:3
% %         for j=1:3
% %             gsO22{i}(p,j)=(conj(gstriplet(:,j)))'*O22*gstriplet(:,p);
% %             gsPxy{i}(p,j)=(conj(gstriplet(:,j)))'*Pxy*gstriplet(:,p);
% %             gsO02{i}(p,j)=(conj(gstriplet(:,j)))'*O02*gstriplet(:,p);
% %         end
% %     end
%  
%     %[realentropyr(i),realFr(i),realgroundr(i)]= calcentropy(Hfinalrelaxed,T1(i));
%     end  
% end
%toc


% figure(50)
% plot(T1,RelaxstrainQ)

% figure(100)
% plot(T1,Sponstrain)
% %   
%  maxstrainQ=max(RelaxstrainQ);
%  minstrainQ=min(RelaxstrainQ);
%  %Finds Tq by finding when the spontaneous Q is half of its max 
%  Tq=max(T1(RelaxstrainQ>0.005*(maxstrainQ-minstrainQ)+minstrainQ));
%  
%  %Concludes there is no Tq is only really small Q is present
%  if maxstrainQ<threshold
%      Tq=0;
%  end

% for z=1:L 
%  Qmat=gsO22{z};
%  Qenergy(z)=Qmat(1,1);
%  Q2energy(z)=Qmat(2,2);
%  Qmix(z)=Qmat(1,2);
%  
% end
% figure(50)
% set=gsenergy(1,3);
% plot(T1,gsenergy(:,1)-set,T1,gsenergy(:,2)-set,T1,gsenergy(:,3)-set)
% figure(200)
% plot(T1,Qenergy,T1,Q2energy)
% 
% figure(300)
% plot(T1,Qmix)
 % hold off
%logT=log(T1);
% logQ=log(ChiQ);
% maxQ=max(ChiQ);
% maxindex=find(ChiQ==maxQ(1));
% L=length(maxindex);
% Curiefit_quadrupolar_suscept=polyfit(logT(1:maxindex(L)),logQ(1:maxindex(L)),1);
% InCuriefit=polyval(Curiefit_quadrupolar_suscept,logT);
% 
% %Displays exp(B)*T^A fit
% Curiefit_quadrupolar_suscept
% 
% options = optimset('Algorithm','levenberg-marquardt');
% X1=[2 7];
% [Curie,res]=lsqnonlin(@fit1,X1,[],[],options,T1(1:M-2),ChiQ(1:M-2));
% TQ=Curie(2)
% A=Curie(1)

% %Figure 1 plots the quad suscept vs temp
% figure (1)
% plot(T1,ChiQ,T1,fit1(Curie,T1,ChiQ-ChiQ));
% title('Quadrupolar Susceptibility');
% ylim([0,9999])
% %Figure 2 plots the loglog quad suscept vs temp and its curie fit
% 
% figure (2)
% plot(logT(1:maxindex),logQ(1:maxindex(L)),logT(1:maxindex(L)),InCuriefit(1:maxindex(L)));
% title(horzcat('Log of Quadrupolar Susceptibility, T(-',num2str(Curiefit_quadrupolar_suscept(1)),')'));
% xlabel('Log(T)')
% ylabel ('Log(chi_Q)')
% 
% 
% 
% figure (3)
% plot(T1,ChiQ,T1,exp(Curiefit_quadrupolar_suscept(2))*T1.^(Curiefit_quadrupolar_suscept(1)),T1,fit(Curie,T1,ChiQ-ChiQ))
% title(num2str(TQ))
% ylim([0 1000])
%expectw(O22,O2_2,Hcef,1)


% figure(4)
% plot(T1,NematicdRdE)
% title('Nematic Resistivity')

% figure (5)
% plot(T1,UnstrainQ)
% title('Unstrained Q moment')

% figure(6)
%     for p=1:3
%         plot(T1,Energy(p,:))
%         title(horzcat('CEF  2J+1=',num2str(2*J+1)))
%         hold on
%     end
%     hold off
end



function diff= fit1(x,X,Y)
    A=x(1);
    B=x(2);
    diff=A*1./(X-B)-Y;
  end


function [exvalue]=expectw(X,Y,Hcef,T)
kb=1;
%exvalue=trace(expm(-1/(kb*T)*Hcef)*X);
[CEF,E]=eig(Hcef);
E=diag(E);
[OrderE,P]=sort(E);
OrderE=OrderE-OrderE(1);

for i=1:length(P)
   Eigen{i}=CEF(:,P(i));
end

Z=sum(exp(-1*OrderE./(kb*T)));

Pe=exp(-1/(kb*T)*OrderE)./Z;


w=zeros(length(Pe));
Sum=0;
for i=1:length(Pe)
    for j=1:length(Pe)
        if OrderE(i)==OrderE(j)
            w(i,j)=1;
        else
        w(i,j)=((OrderE(i)-OrderE(j))/T)/(1-exp((OrderE(j)-OrderE(i))/T));
        end
    Sum=Sum+w(i,j)*Pe(i)*(Eigen{i}'*X*Eigen{j})*(Eigen{j}'*Y*Eigen{i});
    end
end
exvalue=Sum;


end

function [exvalue]=expect(X,Hcef,T)
%kb=1;
%exvalue=trace(expm(-1/(kb*T)*Hcef)*X);
[CEF,E]=eig(Hcef);
E=diag(E);
[OrderE,P]=sort(E);
L=length(P);
%OrderE=OrderE ; %-OrderE(1);
%groundstate=OrderE(1);
Eigen=zeros(L);
for i=1:L
   Eigen(:,i)=CEF(:,P(i));
end

Z=sum(exp(-1/(T)*OrderE));

Pe=exp(-1/(T)*OrderE)/Z;

exvalue=sum(Pe.*diag(Eigen'*X*Eigen));


%Sum=0;
% for i=1:length(Pe)
%     Sum=Sum+Pe(i)*Eigen(:,i)'*X*Eigen(:,i);
% end
% exvalue=Sum;


end

function[entropy, freeenergy,ground]=calcentropy(Hcef,T)
   [energy,Z,ground]=expect(Hcef,Hcef,T);
   freeenergy=-T*log(Z);
   entropy=(energy-freeenergy)/T;
    
end

function [Qrot]=rotateQ(Jz,Q,theta)

R=expm(-1i*theta*Jz);
Rdagger=expm(1i*theta*Jz);
Qrot=Rdagger*Q*R;


end



function [entropyfit,HC]=findHC(entropy,temp,smoothing)
Y=entropy;
X=temp;

[xData, yData] = prepareCurveData( X, Y );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = smoothing;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
entropyfit=feval(fitresult,temp);

HC_T=differentiate(fitresult,temp);

HC=HC_T.*temp';
end
























