%% Rete di adattamento a L tipo 2

clear all; close all; clc;

f = 2e9 % frequenza 2GHz

ZL = 20+30i % impedenza del carico
Z0 = 50 % impedenza caratteristica del carico

YL = 1/ZL % ammettenza del carico
Y0 = 1/Z0 % ammettenza caratteristica del carico

zL = ZL/Z0 % impedenza normalizzata
rL = real(zL) % resistenza normalizzata
xL = imag(zL) % reattanza normalizzata

yL = 1/zL % ammettenza normalizzata
gL = real(yL) % conduttanza normalizzata
bL = imag(yL) % suscettanza normalizzata

GammaL = (yL-1)/(yL+1) % coefficiente di yL
GammaL_abs = abs(GammaL) % modulo di GammaL
GammaL_angle = angle(GammaL)*180/pi % fase di GammaL

PLa = 1-(abs(GammaL))^2 % potenza attiva al carico 68.97%

xx = linspace(0,2*pi,1001);
plot(exp(1i*xx),'k')
axis equal
hold
plot([-1 1],[0 0], 'k')
plot(1/2+1/2*exp(1i*xx),'b')
text(0.05, -0.05,'g=1', 'Color', 'b')
plot(GammaL,'om')
plot(GammaL,'.m')
text(real(GammaL)-0.1, imag(GammaL)+0.05,'y_L')

x = -xL+sqrt(rL*(1-rL)) % segno + -> reattanza negativa
CX = -1/(2*pi*f*(x*Z0)) % elemento jX -> capacità 14.455 pF

Y2 = 1/(ZL+(1i*(x*Z0))) % ammettenza dopo serie tra jX e ZL
y2 = Y2/Y0 % ammettenza normalizzata di Y2
Gamma2 = (y2-1)/(y2+1) % coefficiente di y2

plot(Gamma2,'om')
plot(Gamma2,'.m')
text(real(Gamma2)+0.05, imag(Gamma2)+0.05,'y_2')

PLa = 1-(abs(Gamma2))^2 % potenza attiva al carico 72.73%

b = +sqrt(rL*(1-rL))/rL % segno + -> suscettanza positiva
CB = b/(2*pi*f*Z0) % elemento jB -> capacità 1.9492 pF

Yin = (1i*(b/Z0))+Y2 % ammettenza dopo adattamento = Y0
y_in = Yin/Y0; % ammettenza normalizzata di Yin
Gamma_in = (y_in-1)/(y_in+1) % coefficiente di y_in

plot(Gamma_in,'om')
plot(Gamma_in,'.m')
text(real(Gamma_in)-0.1,imag(Gamma_in)+0.05,'y_{in}')

PLa = 1-(abs(Gamma_in))^2 % potenza attiva al carico 100%

yy = zL+1i*x*linspace(0,1,5);
G = (1-yy)./(1+yy);
plot(G,'.m')
zz = 1/(zL+1i*x)+1i*b*linspace(0,1,20);
G = (zz-1)./(1+zz);
plot(G,'.m')
