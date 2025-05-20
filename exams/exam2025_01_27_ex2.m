%% Esame 27-01-2025 esercizio 2

clear all; close all; clc;

R = 15; % resistenza 15ohm
C = 10*1e-12; % condensatore 10pF
di = 0.64*1e-3; % diametro interno 0.64mm
de = 2.85*1e-3; % diametro esterno 2.85mm
eps_r = 1.43; % permittività relativa del cavo
Pinc = 100*1e-3; % potenza incidente 100mW

c = 3*1e8; % velocità della luce dello spazio libero

a = di/2; % raggio interno
b = de/2; % raggio esterno
Z0 = (60/sqrt(eps_r))*log(b/a) % impedenza caratteristica

%% Frequenza 1 GHz

f1 = 1e9;
w1 = 2*pi*f1;

ZL1 = R + (-1i/(w1*C))
GammaL1 = (ZL1-Z0)/(ZL1+Z0)
GammaL1_abs = abs(GammaL1)
GammaL1_angle = angle(GammaL1)*180/pi
SWR1 = (1+abs(GammaL1))/(1-abs(GammaL1))
PL1 = Pinc*(1-(abs(GammaL1))^2+2*imag(GammaL1)*i)

%% Frequenza 2 GHz

f2 = 2e9;
w2 = 2*pi*f2;

ZL2 = R + (-1i/(w2*C))
GammaL2 = (ZL2-Z0)/(ZL2+Z0)
GammaL2_abs = abs(GammaL2)
GammaL2_angle = angle(GammaL2)*180/pi
SWR2 = (1+abs(GammaL2))/(1-abs(GammaL2))
PL2 = Pinc*(1-(abs(GammaL2))^2+2*imag(GammaL2)*i)

%% Rete di adattamento mediante stub serie

zL = ZL1/Z0
lambda = (c/f1)*(1/sqrt(eps_r))

figure;
xx = linspace(0,2*pi,1001);
plot(exp(1i*xx),'k')
axis equal
hold
plot([-1 1],[0 0], 'k')
plot(1/2+1/2*exp(1i*xx),'r')
text(0.05, -0.05,'r=1', 'Color', 'r')
plot(abs(GammaL1)*sin(xx), abs(GammaL1)*cos(xx),'--k')
plot(GammaL1,'om')
plot(GammaL1,'.m')
text(real(GammaL1)+0.05, imag(GammaL1)+0.05,'z_L')

x1 = (2*abs(GammaL1))/(sqrt(1-abs(GammaL1)^2));
z1 = 1+i*x1
Gammaz1 = (z1-1)/(z1+1)
Gammaz1_abs = abs(Gammaz1)
Gammaz1_angle = angle(Gammaz1)*180/pi

plot(Gammaz1,'om')
plot(Gammaz1,'.m')
text(real(Gammaz1)+0.05, imag(Gammaz1)+0.05,'z_1')

d_lambda = (angle(GammaL1/Gammaz1))/(pi*4)
d = d_lambda*lambda

zs = -1i*x1
Gamma_s = (zs-1)/(zs+1)

plot(Gamma_s,'om')
plot(Gamma_s,'.m')
text(real(Gamma_s)-0.05, imag(Gamma_s)-0.05,'z_s')

l_cc_lambda = 0.25+(atan(1/x1)/(2*pi))
l_cc = l_cc_lambda*lambda

%% Rete di adattamento mediante trasformatore lambda/4

zL = ZL2/Z0
lambda = (c/f2)*(1/sqrt(eps_r))

figure;
xx = linspace(0,2*pi,1001);
plot(exp(1i*xx),'k')
axis equal
hold
plot([-1 1],[0 0], 'k')
plot(1/2+1/2*exp(1i*xx),'r')
text(0.05, -0.05,'r=1', 'Color', 'r')
plot(abs(GammaL2)*sin(xx), abs(GammaL2)*cos(xx),'--k')
plot(GammaL2,'om')
plot(GammaL2,'.m')
text(real(GammaL2)-0.1, imag(GammaL2)+0.05,'z_L')

d_lambda = 0.25+((angle(GammaL2))/(4*pi))
d = d_lambda*lambda

z1 = (1+abs(GammaL2))/(1-abs(GammaL2))
Gamma1 = (z1-1)/(z1+1)

plot(Gamma1,0,'om')
plot(Gamma1,0,'.m')
text(real(Gamma1)-0.1, imag(Gamma1)+0.05,'z_1')

Z1 = z1*Z0
Zt = sqrt(Z0*Z1)
l = lambda/4
