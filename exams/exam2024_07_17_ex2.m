%% Esame 17-07-2024 esercizio 2

clear all; close all; clc;

f0 = 750e6; % frequenza 750MHz
ZL = 125-25i; % impedenza del carico
Z0 = 50; % impedenza caratteristica della linea
eps_r = 2; % permittività relativa del dielettrico
Pinc_dBm = 33; % potenza incidente verso il carico

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
VSWR = (1+abs(GammaL))/(1-abs(GammaL)) % rapporto 2.6180:1

%% Potenza complessa consegnata al carico

Pinc_W = (10^(Pinc_dBm/10))/1000 % 1.9953W
PL = Pinc_W*(1-(abs(GammaL))^2+2*imag(GammaL)*i)
% PL = 1.5962 W - j0.3192 VAR

%% Rete di adattamento mediante stub

% Lunghezza d'onda della linea
c = 3e8; % velocità della luce nello spazio libero
lambda = c/(f0*sqrt(eps_r)) % 0.2828m = 282.8mm

% Si posiziona zL nella carta di Smith e si traccia
% il cerchio VSWR con centro (0,0) e raggio |GammaL|

% Sezione utile dove inserire lo stub è sul cerchio r=1
x1 = (-(2*abs(GammaL))/(sqrt(1-(abs(GammaL))^2)));
z1 = 1+x1*i % impedenza normalizzata sul cerchio r=1

% Si verifica sulla carta di Smith il valore di z1

% Distanza dal carico a cui inserire lo stub
Gamma1 = (z1-1)/(z1+1) % coefficiente di riflessione a z1
d_lambda = (angle(GammaL/Gamma1))/(pi*4) % 0.0737*lambda
d = (angle(GammaL/Gamma1))/(pi*4)*lambda % 0.0209m = 20.9mm

% Si verifica sulla carta di Smith il valore di d_lambda
% usando la ghiera WTG, calcolare la differenza tra zL e z1

zs = -x1*i % impedenza dello stub in serie alla linea

% Si usa una impedenza per raggiungere il centro (0,0)

% Lunghezza dello stub terminato in corto circuito
l_lambda = (1/(2*pi))*atan(-x1) % 0.1250*lambda
l = (lambda/(2*pi))*atan(-x1) % 0.0354m = 35.4mm

% Si usa corto circuito dato che zs è nella parte positiva

%% Modulo del fasore di tensione all'ingresso dello stub

Zs = zs*Z0 % modulo dell'impedenza
I = sqrt((2*Pinc_W)/Z0) % modulo del fasore di corrente
V = abs(I*Zs) % modulo del fasore di tensione
