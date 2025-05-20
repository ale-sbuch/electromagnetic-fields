%% Trasformatore a quarto d'onda con ZL reale

clear all; close all; clc;

ZL = 100 % impedenza del carico
Z0 = 50 % impedenza caratteristica del carico

zL = ZL/Z0 % impedenza di carico normalizzata rispetto a Z0

GammaL = (zL-1)/(zL+1) % coefficiente di zL
GammaL_abs = abs(GammaL) % modulo di GammaL
GammaL_angle = angle(GammaL)*180/pi % fase di GammaL

PLa = 1-(abs(GammaL))^2 % potenza attiva al carico 88.89%

xx = linspace(0,2*pi,1001);
plot(exp(1i*xx),'k')
axis equal
hold
plot([-1 1],[0 0], 'k')
plot(1/2+1/2*exp(1i*xx),'r')
text(0.05, -0.05,'r=1', 'Color', 'r')
plot(GammaL,0,'om')
plot(GammaL,0,'.m')
text(real(GammaL)-0.1, imag(GammaL)+0.05,'z_L')

Zt = sqrt(Z0*ZL) % impedenza del trasformatore

zLt = ZL/Zt % impedenza di carico normalizzata rispetto a Zt
GammaLt = (zLt-1)/(zLt+1) % coefficiente di zLt

plot(GammaLt,0,'om')
plot(GammaLt,0,'.m')
text(real(GammaLt)-0.1, imag(GammaLt)+0.05,'z_{Lt}')

G = GammaLt*exp(-1i*linspace(0,pi,30));
plot(G,'.m')

zt = Zt/ZL % impedenza normalizzata del trasformatore
Gammat = (zt-1)/(zt+1) % coefficiente di zt
