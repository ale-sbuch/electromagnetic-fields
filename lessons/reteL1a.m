%% Rete di adattamento a L tipo 1

clear all; close all; clc;

f = 1e9 % frequenza 1GHz

ZL = 30-30i % impedenza del carico
Z0 = 50 % impedenza caratteristica del carico

YL = 1/ZL % ammettenza del carico
Y0 = 1/Z0 % ammettenza caratteristica del carico

zL = ZL/Z0 % impedenza normalizzata
rL = real(zL) % resistenza normalizzata
xL = imag(zL) % reattanza normalizzata

yL = 1/zL % ammettenza normalizzata
gL = real(yL) % conduttanza normalizzata
bL = imag(yL) % suscettanza normalizzata

GammaL = (zL-1)/(zL+1) % coefficiente di zL
GammaL_abs = abs(GammaL) % modulo di GammaL
GammaL_angle = angle(GammaL)*180/pi % fase di GammaL

PLa = 1-(abs(GammaL))^2 % potenza attiva al carico 82.19%

xx = linspace(0,2*pi,1001);
plot(exp(1i*xx),'k')
axis equal
hold
plot([-1 1],[0 0], 'k')
plot(1/2+1/2*exp(1i*xx),'r')
text(0.05, -0.05,'r=1', 'Color', 'r')
plot(GammaL,'om')
plot(GammaL,'.m')
text(real(GammaL)-0.1, imag(GammaL)+0.05,'z_L')

b = -bL+sqrt(gL*(1-gL)) % segno + -> suscettanza negativa
LB = -Z0/(2*pi*f*b) % elemento jB -> induttanza 17.275 nH

Z1 = 1/(YL+1i*(b/Z0)) % impedenza dopo parallelo tra jB e ZL
z1 = Z1/Z0 % impedenza normalizzata di Z1
Gamma1 = (z1-1)/(z1+1) % coefficiente di z1

plot(Gamma1,'om')
plot(Gamma1,'.m')
text(real(Gamma1)-0.1, imag(Gamma1)+0.05,'z_1')

PLa = 1-(abs(Gamma1))^2 % potenza attiva al carico 95.24%

x = +sqrt(gL*(1-gL))/gL % segno + -> reattanza positiva
LX = (x*Z0)/(2*pi*f) % elemento jX -> induttanza 3.5588 nH

Zin = (1i*(x*Z0))+Z1 % impedenza dopo adattamento Zin = Z0
z_in = Zin/Z0; % impedenza normalizzata di Zin
Gamma_in = (z_in-1)/(z_in+1) % coefficiente di z_in

plot(Gamma_in,'om')
plot(Gamma_in,'.m')
text(real(Gamma_in)-0.1,imag(Gamma_in)+0.05,'z_{in}')

PLa = 1-(abs(Gamma_in))^2 % potenza attiva al carico 100%

yy = yL+1i*b*linspace(0,1,10);
G = (1-yy)./(1+yy);
plot(G,'.m')
zz = 1/(yL+1i*b)+1i*x*linspace(0,1,10);
G = (zz-1)./(1+zz);
plot(G,'.m')
