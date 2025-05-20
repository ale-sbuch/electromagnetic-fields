%% Rete di adattamento mediante stub parallelo

clear all; close all; clc;

f = 3e9 % frequenza 3GHz

ZL = 100-200i % impedenza del carico
Z0 = 50 % impedenza caratteristica del carico

eps_r = 2.2; % permittività relativa del dielettrico

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

PLa = 1-(abs(GammaL))^2 % potenza attiva al carico 32.00%

VSWR = (1+abs(GammaL))/(1-abs(GammaL)) % S = 10.4039:1

c = 3e8; % velocità della luce nello spazio libero
lambda = c/f/sqrt(eps_r) % 0.0674 m = 67.4 mm

figure;
title('Rete di adattamento mediante stub parallelo')
xx = linspace(0,2*pi,1001);
plot(exp(1i*xx),'k')
axis equal
hold
plot([-1 1],[0 0], 'k')
plot(1/2+1/2*exp(1i*xx),'b')
text(0.05, -0.05,'g=1', 'Color', 'b')
plot(2/3+1/3*exp(i*xx),'b')
text(0.4, -0.05,'g=2', 'Color', 'b')
plot(abs(GammaL)*sin(xx), abs(GammaL)*cos(xx),'--k')
plot(GammaL,'om')
plot(GammaL,'.m')
text(real(GammaL)-0.1, imag(GammaL)+0.05,'y_L')

b1 = (2*abs(GammaL))/(sqrt(1-(abs(GammaL))^2)) % segno +
y1 = 1+1i*b1 % impendenza interseca il cerchio g=1
Gamma1 = (y1-1)/(y1+1) % coefficiente di y1

plot(Gamma1,'om')
plot(Gamma1,'.m')
text(real(Gamma1)+0.05, imag(Gamma1)+0.05,'y_1')

d_lambda = (angle(GammaL/Gamma1))/(4*pi) % 0.1704
d = d_lambda*lambda % 0.0115 m = 11.5 mm

ys = -1i*b1 % ammettenza di ingresso normalizzata
Gamma_s = (ys-1)/(ys+1) % coefficiente di ys

plot(Gamma_s,'om')
plot(Gamma_s,'.m')
text(real(Gamma_s)-0.1, imag(Gamma_s)+0.05,'y_s')

l_cc_lambda = atan(1/b1)/(2*pi) % 0.0526 corto circuito
l_cc = l_cc_lambda*lambda % 0.0035 m = 3.5 mm

l_ca_lambda = 0.5+(atan(-b1)/(2*pi)) % 0.3026 circuito aperto
l_ca = l_ca_lambda*lambda % 0.0204 m = 20.4 mm
