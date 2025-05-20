%% Esame 27-01-2025 esercizio 1

clear all; close all; clc

f0=15*1e9; % frequenza 15 GHz
w0=2*pi*f0; % frequenza angolare

thetainc=pi/3; % angolo di incidenza

eps0=8.854*1e-12; % permittività nello spazio libero
mu0=4*pi*1e-7; % permeabilità nello spazio libero

% mezzo 1
epsr1=1; mur1=1;
eps1=epsr1*eps0;
mu1=mur1*mu0;

% mezzo 2
epsr2=3; mur2=1;
eps2=epsr2*eps0;
mu2=mur2*mu0;

% numero d'onda
k1=w0*sqrt(eps1*mu1)
k2=w0*sqrt(eps2*mu2)

% impedenza intrinseca
eta1=sqrt(mu1/eps1)
eta2=sqrt(mu2/eps2)

% angolo di trasmissione
thetatra=asin(sqrt((eps1*mu1)/(eps2*mu2))*sin(thetainc))
thetatra_degree=thetatra*180/pi % cos(thetatra) = sqrt(3)/2

% coefficiente di riflessione parallelo e perpendicolare
Gammapar=(eta1*cos(thetainc)-eta2*cos(thetatra))./(eta1*cos(thetainc)+eta2*cos(thetatra))
Gammaprp=(eta2./cos(thetatra)-eta1./cos(thetainc))./(eta2./cos(thetatra)+eta1./cos(thetainc))

% % coefficiente di trasmissione parallelo e perpendicolare
Taupar=2*eta2*cos(thetainc)./(eta1*cos(thetainc)+eta2*cos(thetatra))
Tauprp=2*eta2./cos(thetatra)./(eta2./cos(thetatra)+eta1./cos(thetainc))
