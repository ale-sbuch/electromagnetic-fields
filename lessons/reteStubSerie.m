%% Rete di adattamento mediante stub serie

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

GammaL = (zL-1)/(zL+1) % coefficiente di zL
GammaL_abs = abs(GammaL) % modulo di GammaL
GammaL_angle = angle(GammaL)*180/pi % fase di GammaL

PLa = 1-(abs(GammaL))^2 % potenza attiva al carico 32.00%

VSWR = (1+abs(GammaL))/(1-abs(GammaL)) % S = 10.4039:1

c = 3e8; % velocità della luce nello spazio libero
lambda = c/f/sqrt(eps_r) % 0.0674 m = 67.4 mm

figure;
xx = linspace(0,2*pi,1001);
plot(exp(1i*xx),'k')
axis equal
hold
plot([-1 1],[0 0], 'k')
plot(1/2+1/2*exp(1i*xx),'r')
text(0.05, -0.05,'r=1', 'Color', 'r')
plot(2/3+1/3*exp(i*xx),'r')
text(0.4, -0.05,'r=2', 'Color', 'r')
plot(abs(GammaL)*sin(xx), abs(GammaL)*cos(xx),'--k')
plot(GammaL,'om')
plot(GammaL,'.m')
text(real(GammaL)-0.1, imag(GammaL)+0.05,'z_L')

x1 = -(2*abs(GammaL))/(sqrt(1-(abs(GammaL))^2)) % segno -
z1 = 1+1i*x1 % impendenza interseca il cerchio r=1
Gamma1 = (z1-1)/(z1+1) % coefficiente di z1

plot(Gamma1,'om')
plot(Gamma1,'.m')
text(real(Gamma1)-0.1, imag(Gamma1)+0.05,'z_1')

d_lambda = (angle(GammaL/Gamma1))/(4*pi) % 0.0161
d = d_lambda*lambda % 0.0011 m = 1.1 mm

zs = -1i*x1 % impedenza di ingresso normalizzata
Gamma_s = (zs-1)/(zs+1) % coefficiente di zs

plot(Gamma_s,'om')
plot(Gamma_s,'.m')
text(real(Gamma_s)-0.05, imag(Gamma_s)-0.05,'z_s')

l_cc_lambda = atan(-x1)/(2*pi) % 0.1974 corto circuito
l_cc = l_cc_lambda*lambda % 0.0133 m = 13.3 mm

l_ca_lambda = 0.5+(atan(1/x1)/(2*pi)) % 0.4474 circuito aperto
l_ca = l_ca_lambda*lambda % 0.0302 m = 30.2 mm
