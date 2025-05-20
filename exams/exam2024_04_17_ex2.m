%% Esame 17-04-2024 esercizio 2

clear all; close all; clc;

f0 = 750e6; % frequenza 750MHz
ZL = 25-100i; % impedenza del carico
Z0 = 50; % impedenza caratteristica della linea
eps_r = 2; % permittività relativa del dielettrico
Pinc_dBm = 17; % potenza incidente verso il carico

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
VSWR = (1+abs(GammaL))/(1-abs(GammaL)) % rapporto 10.4039:1

%% Potenza complessa consegnata al carico

Pinc_W = (10^(Pinc_dBm/10))/1000 % 0.0501mW = 50.1W
PL = Pinc_W*(1-(abs(GammaL))^2+2*imag(GammaL)*i)
% PL = 0.0160 W - j0.0642 VAR = 16.0mW - j64.2mVAR

%% Rete di adattamento mediante stub

% Lunghezza d'onda della linea
c = 3e8; % velocità della luce nello spazio libero
lambda = c/(f0*sqrt(eps_r)) % 0.2828m = 282.8mm

% Si posiziona zL nella carta di Smith e si traccia
% il cerchio VSWR con centro (0,0) e raggio |GammaL|

% Sezione utile dove inserire lo stub è sul cerchio g=1
b1 = (2*abs(GammaL))/(sqrt(1-(abs(GammaL))^2))
y1 = 1+b1*i % impedenza normalizzata sul cerchio r=1

% Si verifica sulla carta di Smith il valore di z1=1/y1

% Distanza dal carico a cui inserire lo stub
Gamma1 = (-b1*i)/(2+b1*i) % coefficiente di riflessione a y1
d_lambda = (angle(GammaL/Gamma1))/(pi*4) % 0.1314*lambda
d = (angle(GammaL/Gamma1))/(pi*4)*lambda % 0.0372m = 37.2mm

% Si verifica sulla carta di Smith il valore di d_lambda
% usando la ghiera WTG, calcolare la differenza tra zL e z1

ys = -b1*i % suscettanza dello stub in parallelo alla linea

% Si usa una suscettanza per raggiungere il centro (0,0)

% Lunghezza dello stub terminato in corto circuito
l_lambda = (1/(2*pi))*atan(1/b1) % 0.0525*lambda
l = (lambda/(2*pi))*atan(1/b1) % 0.0149m = 14.9mm

% Si usa corto circuito dato che zs=1/ys è nella parte positiva

%% Potenza complessa all'ingresso dello stub

Ps = Pinc_W*conj(ys) % 0W + j0.1461VAR = j146.1mVAR
