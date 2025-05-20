%% Calcola i coefficienti di riflessione e trasmissione
%% tra 2 mezzi al variare dell'angolo di incidenza

clear all; close all; clc

%% Mezzo 1 più denso del mezzo 2
epsr1=3; mur1=1;
epsr2=1; mur2=1;

%% Mezzo 1 meno denso del mezzo 2
% epsr1=1; mur1=1;
% epsr2=3; mur2=1;

%% Mezzo 1 meno denso del mezzo 2 con piccole perdite
% epsr1=1; mur1=1;
% epsr2=3-0.03i; mur2=1;

%% Mezzo 1 meno denso del mezzo 2 con perdite
% epsr1=1; mur1=1;
% epsr2=3-0.3i; mur2=1;

%% Mezzo 1 è aria e mezzo 2 è un conduttore (rame) a 1GHz
% epsr1=1; mur1=1;
% sigma=5.8e7; f=1e9; omega=2*pi*f; eps0=8.854e-12;
% epsr2=1-1i*sigma/omega/eps0; mur2=1;

%% Plots

% impedenze normalizzate
eta1=sqrt(mur1/epsr1);
eta2=sqrt(mur2/epsr2);

% angolo di incidenza
thetainc=linspace(0,pi/2,1001);

% angolo di trasmissione dalla legge di Snell
thetatra=conj(asin(sqrt(epsr1*mur1)/sqrt(epsr2*mur2)*sin(thetainc)));

% coefficiente di riflessione parallelo e perpendicolare
Gammapar=(eta1*cos(thetainc)-eta2*cos(thetatra))./(eta1*cos(thetainc)+eta2*cos(thetatra));
Gammaprp=(eta2./cos(thetatra)-eta1./cos(thetainc))./(eta2./cos(thetatra)+eta1./cos(thetainc));

% coefficiente di trasmissione parallelo e perpendicolare
Taupar=2*eta2*cos(thetainc)./(eta1*cos(thetainc)+eta2*cos(thetatra));
Tauprp=2*eta2./cos(thetatra)./(eta2./cos(thetatra)+eta1./cos(thetainc));

figure(1)
plot(thetainc,abs(Gammapar),thetainc,abs(Gammaprp),thetainc,abs(Taupar),thetainc,abs(Tauprp),'linewidth',2)
xlabel('\theta^{inc} (rad)'), ylabel('|\Gamma|_{||} or |\Gamma|_{\perp} or |\tau|_{||} or |\tau|_{\perp}')
legend('\Gamma_{||}','\Gamma_{\perp}','\tau_{||}','\tau_{\perp}')

figure(2)
plot(thetainc,real(Gammapar),thetainc,real(Gammaprp),thetainc,imag(Gammapar),thetainc,imag(Gammaprp),'linewidth',2)
xlabel('\theta^{inc} (rad)'), ylabel('\Gamma')
legend('\Ree\{\Gamma_{||}\}','\Ree\{\Gamma_{\perp}\}','\Imm\{\Gamma_{||}\}','\Imm\{\Gamma_{\perp}\}')

figure(3)
plot(thetainc,real(Taupar),thetainc,real(Tauprp),thetainc,imag(Taupar),thetainc,imag(Tauprp),'linewidth',2)
xlabel('\theta^{inc} (rad)'), ylabel('\tau')
legend('\Ree\{\tau_{||}\}','\Ree\{\tau_{\perp}\}','\Imm\{\tau_{||}\}','\Imm\{\tau_{\perp}\}')
