%% Trasformatore a quarto d'onda con ZL complesso

clear all; close all; clc;

ZL = 100+100i % impedenza del carico
Z0 = 50 % impedenza caratteristica del carico

zL = ZL/Z0 % impedenza di carico normalizzata rispetto a Z0

GammaL = (zL-1)/(zL+1) % coefficiente di zL
GammaL_abs = abs(GammaL) % modulo di GammaL
GammaL_angle = angle(GammaL)*180/pi % fase di GammaL

PLa = 1-(abs(GammaL))^2 % potenza attiva al carico 39.02%

VSWR = (1+abs(GammaL))/(1-abs(GammaL)) % S = 8.1270:1

xx = linspace(0,2*pi,1001);
plot(exp(1i*xx),'k')
axis equal
hold
plot([-1 1],[0 0], 'k')
plot(1/2+1/2*exp(1i*xx),'r')
text(0.05, -0.05,'r=1', 'Color', 'r')
plot(abs(GammaL)*sin(xx), abs(GammaL)*cos(xx),'--k')
plot(GammaL,'om')
plot(GammaL,'.m')
text(real(GammaL)-0.1, imag(GammaL)+0.05,'z_L')

d_lambda = 0.25-((angle(GammaL))/(4*pi)) % 0.2087*lambda

z1 = (1+abs(GammaL))/(1-abs(GammaL)) % impedenza a d1
Gamma1 = (z1-1)/(z1+1) % coefficiente di z1

plot(Gamma1,0,'om')
plot(Gamma1,0,'.m')
text(real(Gamma1)-0.1, imag(Gamma1)+0.05,'z_1')

d1_lambda = (angle(GammaL/Gamma1))/(4*pi) % 0.0413*lambda

zLt = sqrt(z1) % impedenza normalizzata zL rispetto a Zt
GammaLt = (zLt-1)/(zLt+1) % coefficiente di zLt

plot(GammaLt,0,'om')
plot(GammaLt,0,'.m')
text(real(GammaLt)-0.1, imag(GammaLt)+0.05,'z_{Lt}')

G = GammaLt*exp(-1i*linspace(0,pi,100));
plot(G,'.m')

zt = 1/zLt % impedenza normalizzata del trasformatore
Gammat = (zt-1)/(zt+1) % coefficiente di zt

plot(Gammat,0,'om')
plot(Gammat,0,'.m')
text(real(Gammat)-0.1, imag(Gammat)+0.05,'z_t')

Z1 = z1*Z0; % impedenza Z1 rispetto a Z0
Zt = sqrt(Z0*Z1) % impedenza del trasformatore
