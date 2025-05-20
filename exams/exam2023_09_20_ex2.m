%% Esame 20-09-2023 esercizio 2

clear all; close all; clc;

f0 = 1.2e9; % frequenza 1.2GHz
ZL = 5+25i; % impedenza del carico
Z0 = 50; % impedenza caratteristica della linea
eps_r = 3; % permittività relativa del dielettrico
Pinc_dBm = 10; % potenza incidente verso il carico

%% Impedenza normalizzata e ammettenza normalizzata

zL = ZL/Z0 % impedenza normalizzata
rL = real(zL); % resistenza normalizzata
xL = imag(zL); % reattanza normalizzata

YL = 1/ZL; % ammettenza del carico
Y0 = 1/Z0; % ammettenza caratteristica della linea

yL = YL/Y0 % ammettenza normalizzata
gL = real(yL); % conduttanza normalizzata
bL = imag(yL); % suscettanza normalizzata

%% Rapporto d'onda stazionaria in tensione

GammaL = (zL-1)/(zL+1) % coefficiente di riflessione a zL
VSWR = (1+abs(GammaL))/(1-abs(GammaL)) % rapporto 12.5201:1

%% Potenza complessa consegnata al carico

Pinc_W = (10^(Pinc_dBm/10))/1000 % 0.0100W = 10mW
PL = Pinc_W*(1-(abs(GammaL))^2+2*imag(GammaL)*i)
% PL = 0.0027 W + j0.0137 VAR = 2.7mW + j13.7mVAR

%% Rete di adattamento mediante stub

% Lunghezza d'onda della linea
c = 3e8; % velocità della luce nello spazio libero
lambda = c/(f0*sqrt(eps_r)) % 0.1443m = 144.3mm

% Si posiziona zL nella carta di Smith e si traccia
% il cerchio VSWR con centro (0,0) e raggio |GammaL|

% Sezione utile dove inserire lo stub è sul cerchio r=1
x1 = (2*abs(GammaL)/(sqrt(1-(abs(GammaL))^2)));
z1 = 1+x1*i % impedenza normalizzata sul cerchio r=1

% Si verifica sulla carta di Smith il valore di z1

% Distanza dal carico a cui inserire lo stub
Gamma1 = (z1-1)/(z1+1) % coefficiente di riflessione a z1
d_lambda = (angle(GammaL/Gamma1))/(pi*4) % 0.1319*lambda
d = (angle(GammaL/Gamma1))/(pi*4)*lambda % 0.0190m = 19.0mm

% Si verifica sulla carta di Smith il valore di d_lambda
% usando la ghiera WTG, calcolare la differenza tra zL e z1

zs = -x1*i % impedenza dello stub in serie alla linea

% Si usa una impedenza per raggiungere il centro (0,0)

% Lunghezza dello stub terminato in circuito aperto
l_lambda = (1/(2*pi))*atan(1/x1) % 0.0474*lambda
l = (lambda/(2*pi))*atan(1/x1) % 0.0068m = 6.8mm

% Si usa circuito aperto dato che zs è nella parte negativa

%% Potenza reattiva all'ingresso dello stub

Ps = -Pinc_W*x1*i % -j0.0326VAR = -j32.6mVAR
