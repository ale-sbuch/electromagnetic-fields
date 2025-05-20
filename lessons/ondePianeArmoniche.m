%% Onde piane armoniche

clear all; close all; clc

epsr=-1i;
sqreps=-1i*sqrt(-epsr);
k=2*pi*sqreps; % lambda=1
u=[0;0;1]; 
E0=[1;0;0];
eta=1/sqreps;
H0=cross(u,E0)/eta;
z=linspace(0,2,73);
x=zeros(size(z)); y=x;
r=[x;y;z];
Efas=E0*exp(-1i*k*(u(1)*x+u(2)*y+u(3)*z));
Hfas=H0*exp(-1i*k*(u(1)*x+u(2)*y+u(3)*z));
Nt=36*4;
wt=linspace(0,2*pi,Nt+1);
figure(1)
for nt=1:Nt
    Etd=real(Efas*exp(1i*wt(nt)));
    Htd=real(Hfas*exp(1i*wt(nt)));
    quiver3(x,y,z,Etd(1,:),Etd(2,:),Etd(3,:),0,'b')
    axis([-1 1 -1 1 0 2])
    hold on
    quiver3(x,y,z,Htd(1,:),Htd(2,:),Htd(3,:),0,'r')
    xlabel('x'); ylabel('y'); zlabel('z')
    hold off
    Mov(nt)=getframe;
end
movie(Mov,10,2*Nt)
