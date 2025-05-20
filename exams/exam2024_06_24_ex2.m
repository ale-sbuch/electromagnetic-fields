%% Esame 24-06-2024 esercizio 2

clear all; close all; clc;

%% Dati del problema

f = 800e6; % frequenza 800MHz
ZL = 25+100i; % impedenza del carico
Z0 = 50; % impedenza caratteristica del carico
eps_r = 2.2; % permittività relativa del dielettrico
Pinc_dBm = 23; % potenza incidente

%% Impedenza e ammettenza

YL = 1/ZL % ammettenza del carico
Y0 = 1/Z0 % ammettenza caratteristica del carico

zL = ZL/Z0 % impedenza normalizzata
rL = real(zL) % resistenza normalizzata
xL = imag(zL) % reattanza normalizzata

yL = 1/zL % ammettenza normalizzata
gL = real(yL) % conduttanza normalizzata
bL = imag(yL) % suscettanza normalizzata

%% Rapporto d'onda stazionaria in tensione

GammaL = (zL-1)/(zL+1) % coefficiente di riflessione

abs(GammaL)
angle(GammaL) % radianti
angle(GammaL)*180/pi % gradi



VSWR = (1+abs(GammaL))/(1-abs(GammaL)) % rapporto 10.4:1

%% Potenza consegnata al carico

Pinc_W = (10^(Pinc_dBm/10))/1000 % potenza incidente 200mW
PL = Pinc_W*(1-(abs(GammaL))^2+2*imag(GammaL)*i)
% PL = 0.0638 W + j0.2554 VAR = 63.8mW + j255.4 mVAR
