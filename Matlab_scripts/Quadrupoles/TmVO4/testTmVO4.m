clear all

N=150;

theta=0:2*pi/N:2*pi;
phi=0:pi/N:pi;

tth=repmat(theta,length(phi),1);
pph=repmat(phi',1,length(theta));

    J=6;

    ff=1;
    
    % %      6   5    4 3   2   1   0    1   2   3   4   5   6 
    vec1=[    0   0.89 0 0   0 -0.42 0  0    0   0.19 0   0   0 ];
    vec2=vec1(end:-1:1);
    
    vec=vec1+exp(pi*1i*-0.5)*vec2;
    vec=vec/norm(vec);
    label='Tm3+';
    
    rhoeTot=densityElectronicAll(J,vec,tth,pph,label);
    rhoe0=densityElectronicAll(J,vec,tth,pph,label,[0]);
    rhoeb=densityElectronicAll(J,vec,tth,pph,label,[0,2]);
    rhoe2=densityElectronicAll(J,vec,tth,pph,label,[2]);
    rhoe4=densityElectronicAll(J,vec,tth,pph,label,[4]);
    rhoe6=densityElectronicAll(J,vec,tth,pph,label,[6]);
    
    Xet=abs(rhoeTot).*sin(pph).*cos(tth);
    Yet=abs(rhoeTot).*sin(pph).*sin(tth);
    Zet=abs(rhoeTot).*cos(pph);
    
    Xe0=abs(rhoe0).*sin(pph).*cos(tth);
    Ye0=abs(rhoe0).*sin(pph).*sin(tth);
    Ze0=abs(rhoe0).*cos(pph);
    
        Xeb=abs(rhoeb).*sin(pph).*cos(tth);
    Yeb=abs(rhoeb).*sin(pph).*sin(tth);
    Zeb=abs(rhoeb).*cos(pph);
    
    Xe2=ff*abs(rhoe2).*sin(pph).*cos(tth);
    Ye2=ff*abs(rhoe2).*sin(pph).*sin(tth);
    Ze2=ff*abs(rhoe2).*cos(pph);
    
    Xe4=ff*abs(rhoe4).*sin(pph).*cos(tth);
    Ye4=ff*abs(rhoe4).*sin(pph).*sin(tth);
    Ze4=ff*abs(rhoe4).*cos(pph);
    
    Xe6=ff*abs(rhoe6).*sin(pph).*cos(tth);
    Ye6=ff*abs(rhoe6).*sin(pph).*sin(tth);
    Ze6=ff*abs(rhoe6).*cos(pph);
    
    
 subplot(2,3,1); 
    surf(Xet,Yet,Zet,real(rhoeTot));    
    axis equal
    shading flat
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    caxis([11.3 12.7])
    
    camlight
    camlight
    lighting phong
 %   axis off
    colorbar
    set(gca,'fontsize',14);
    
    
     subplot(2,3,2); 
    surf(Xe0,Ye0,Ze0,real(rhoe0));    
    axis equal
    shading flat
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    caxis([11.3 12.7])
    
    camlight
    camlight
    lighting phong
 %   axis off
    colorbar
    set(gca,'fontsize',14);
    
         subplot(2,3,3); 
    surf(Xe2,Ye2,Ze2,real(rhoe2));    
    axis equal
    shading flat
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    caxis([11.3 12.7]-12)
    
    camlight
    camlight
    lighting phong
  %  axis off
    colorbar
    set(gca,'fontsize',14);
    
         subplot(2,3,4); 
    surf(Xe4,Ye4,Ze4,real(rhoe4));    
    axis equal
    shading flat
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    caxis([11.3 12.7]-12)
    
    camlight
    camlight
    lighting phong
   % axis off
    colorbar
    set(gca,'fontsize',14);
    
         subplot(2,3,5); 
    surf(Xe6,Ye6,Ze6,real(rhoe6));    
    axis equal
    shading flat
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    caxis([11.3 12.7]-12)
    
    camlight
    camlight
    lighting phong
   % axis off
    colorbar
    set(gca,'fontsize',14);
    
    
             subplot(2,3,6); 
    surf(Xeb,Yeb,Zeb,real(rhoeb));    
    axis equal
    shading flat
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    caxis([11.3 12.7])
    
    camlight
    camlight
    lighting phong
   % axis off
    colorbar
    set(gca,'fontsize',14);
    
    
    
 %%

NN=60;
fin=[0,1,0];
mid=[0.6,0.6,0.6];
ini=[1,0,0];

R=[ini(1):(mid(1)-ini(1))/NN:mid(1) mid(1):(fin(1)-mid(1))/NN:fin(1)];
G=[ini(2):(mid(2)-ini(2))/NN:mid(2) mid(2):(fin(2)-mid(2))/NN:fin(2)];
B=[ini(3):(mid(3)-ini(3))/NN:mid(3) mid(3):(fin(3)-mid(3))/NN:fin(3)];
col=[R; G; B];


colormap(col')


    

