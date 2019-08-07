
%Hamiltonian = 2*muBohr*SigmaZ+ 2e-5 * expSx*SigmaX;

Trange= 0.01:0.01:2.5;
Hrange= 0.01:0.1:5;
[TGrid,HGrid] = meshgrid(Trange,Hrange);


SxVals=[];
SzVals=[];
EVals=[];
for T = Trange
    T
    SxValsH=[];
    SzValsH=[];
    EValsH=[];
    for H = Hrange
        expSx = 0;
        newSx = 1;
        while abs(expSx-newSx) > 1e-6;
            expSx = newSx;
            Hamiltonian = 15e-5*expSx*SigmaX+2*muBohr*SigmaZ*H;
            [ newSx] = ThermalExpectation(T, Hamiltonian, SigmaX);
        end
        SxValsH=[SxValsH newSx];
        SzValsH=[SzValsH ThermalExpectation(T, Hamiltonian, SigmaZ)];
        EValsH=[EValsH ThermalExpectation(T, Hamiltonian, Hamiltonian)];
    end
    SxVals= [SxVals SxValsH']; 
    SzVals= [SzVals SzValsH']; 
    EVals= [EVals EValsH']; 
end
%%
imagesc(Trange,Hrange,SxVals)
xlabel('Temperature')
ylabel('Field')
set(gca,'YDir','normal')
%%
contour(Trange,Hrange,SxVals)
xlabel('Temperature')
ylabel('Field')
%%
imagesc(Hrange,Trange(2:end),diff(SzVals,1,2))
colorbar

set(gca,'YDir','normal')
xlabel('Temperature')
ylabel('Field')

%%
imagesc(Hrange,Trange(2:end),EVals)
colorbar
set(gca,'YDir','normal')
xlabel('Temperature')
ylabel('Field')
title('<E>')
%%
dE=diff(EVals,1,2);
imagesc(Hrange,Trange(2:end),-(dE))
colorbar
set(gca,'YDir','normal')
xlabel('Temperature')
ylabel('Field')
title('d<E>/dT')
%%
surf(HGrid,TGrid,-EVals,'EdgeColor','none')
ylabel('Temperature')
xlabel('Field')
zlabel('<E>')
%%
surf(HGrid(:,2:end),TGrid(:,2:end),-diff(EVals,1,2),'EdgeColor','none')
ylabel('Temperature')
xlabel('Field')
zlabel('Cp')

%%
mask=TGrid>0.35 & HGrid>0.03 & TGrid<3;
surf(HGrid(:,2:end),TGrid(:,2:end),TGrid(:,2:end).*diff(SzVals,1,2)./(-diff(EVals,1,2) + 1e-2).*mask(:,1:(end-1)),'EdgeColor','none')
ylabel('Temperature')
xlabel('Field')
zlabel('Cp')

%%
imagesc(Hrange,Trange(2:end),TGrid(:,2:end).*diff(SzVals,1,2)./(-diff(EVals,1,2) + 1e-2).*mask(:,1:(end-1)))
ylabel('Temperature')
xlabel('Field')
zlabel('Cp')

set(gca,'YDir','normal')
%%
plot(-dE(:,50))