%% File info
% Main script for computation of multipoles in TmVO4.
% Derived from 'testTmVO4'.
% All scripts originally written by Nicolas Gauthier from Shen group and
% shared with me in August 2019.

%% Change working directory
cd 'C:\Users\Pierre\Desktop\Postdoc\Software\Matlab\Matlab_scripts\Multipoles'

%% Compute and plot multipoles
clear all

N=150;

theta=0:2*pi/N:2*pi;
phi=0:pi/N:pi;

tth=repmat(theta,length(phi),1);
pph=repmat(phi',1,length(theta));

    J=6;

    ff=1;
    
    % %      6   5    4 3   2   1   0    1   2   3   4   5   6 
    vec1=[   0   0.89 0 0   0 -0.42 0    0   0  0.19 0   0   0 ];% From Hodges et al. 1983, p.5
    vec2=vec1(end:-1:1);% vec1 and vec2 are eigenvectors of Jz and hence Jz^2
% One can show that eigenvectors of O_2^2 = Jx^2 - Jy^2 = J+^2 + J-^2 are
% vec1 +- vec2 (with appropriate normalization)
% and eigenvectors of Pxy = JxJy + JyJx = -i/4.(J+^2 - J-^2) are 
% vec1 +- i.vec2 (with appropriate normalization)
    
%     vec=vec1+exp(pi*1i*-0.5)*vec2;% For wavefunction in the orthorhombic phase
% %     Change exp(pi*1i*0.5) to exp(pi*1i*-0.5) to toggle between two possible wavefunctions
    vec = vec1;% For wavefunction in the tetragonal phase
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
    
    
    plt(1).ax = subplot(2,3,1); 
    plt(1).surf = surf(Xet,Yet,Zet,real(rhoeTot));    
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
    plt(1).cb = colorbar;
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
    
    plt(3).ax = subplot(2,3,3); 
    plt(3).surf = surf(Xe2,Ye2,Ze2,real(rhoe2));    
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
    plt(3).cb = colorbar;
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
    

%% Duplicate a given subplot into an individual figure
fig = figure;
idx = 1;
ax = plt(idx).ax;
% axf = subplot(1,1,1);
% axf = copyobj(ax,fig); %copy axes to figure
copies = copyobj([plt(idx).cb,plt(idx).ax],fig); 
cb = copies(1); axf = copies(2);
sf = axf.Children(end);% allow control of plot, e.g. sf.Visible = 'off';
sf.DiffuseStrength = 1;% Intensity of the diffuse component of the light reflected from the object. 
sf.AmbientStrength = .5;% Intensity of the ambient component of the light reflected from the object. 
sf.SpecularStrength = .75;% Intensity of the specular component of the light reflected from the object. 
% set(hNew, 'pos', [0.23162 0.2233 0.72058 0.63107])
cb.Visible = 'off';% sf.Visible = 'off';
% hNew.SortMethod='ChildOrder';
% l = findobj(gcf,'Type','Light'); l.delete;
axf.XTick=[]; axf.YTick=[]; axf.ZTick=[];
% axf.XLabel.String=''; axf.YLabel.String=''; axf.ZLabel.String='';
axf.FontSize = 48;% Increase fontsize, especially for axis labels
view(2);

% % % colorbar options % % % 
% cb.Location = 'south'; % Orient colorbar horizontally by locating it at the bottom of the plot
% cb.Position(1) = 0; cb.Position(3:4) = [1 .1];% Set position and shape of colorbar
% cb.Ticks = [11.5 12 12.5]; cb.TickLength = .02; cb.TickLabelInterpreter = 'latex';% modify ticks display
% cb.Label.String = '$\rho$ (a.u.)'; cb.Label.Interpreter = 'latex'; cb.FontSize = 24;% format colorbar label

formatFigure;

%% Figure Post-processing 
lgt = findobj(gcf,'Type','Light');% find light objects
lgt(:).delete;% delete all light objects
lgt_origin = 'headlight';% Origin of camlight
camlight(lgt_origin);% light figure twice with headlight (no obvious difference when using 'infinite' option)
formatFigure;

%% Change surf colors

NN=60;
fin=[0,1,0];
mid=[0.6,0.6,0.6];
ini=[1,0,0];

R=[ini(1):(mid(1)-ini(1))/NN:mid(1) mid(1):(fin(1)-mid(1))/NN:fin(1)];
G=[ini(2):(mid(2)-ini(2))/NN:mid(2) mid(2):(fin(2)-mid(2))/NN:fin(2)];
B=[ini(3):(mid(3)-ini(3))/NN:mid(3) mid(3):(fin(3)-mid(3))/NN:fin(3)];
col=[R; G; B];


colormap(col')

%% Export figure
% printPDF('2019-08-09_TmVO4_multipoles')
export_fig '2019-10-29_TmVO4_full_electronic_density_cb' -transparent -r150
% export figure in png format with transparent background and resolution 150 dpi

    

