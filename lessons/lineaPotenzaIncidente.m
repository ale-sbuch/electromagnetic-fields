%% Tensione, corrente e potenza incidente di una linea

clear all; close all; clc

beta=2*pi;
z=linspace(-2,0,501);

Z0=50;
ZL=30+70i;

GammaL=(ZL-Z0)/(ZL+Z0);
Gamma=GammaL*exp(2i*beta*z);

Vp=1;
Pinc=1/2*Vp^2/Z0;
V=Vp*exp(-1i*beta*z).*(1+Gamma);
I=Vp/Z0*exp(-1i*beta*z).*(1-Gamma);
P=1/2*V.*conj(I);

plot(z,abs(V),z,abs(I)*Z0,z,real(P/Pinc),z,imag(P/Pinc))
legend('|V|','|I|','\Re\{P/P^{inc}\}','\Im\{P/P^{inc}\}')

wt=linspace(0,2*pi,300);

figure(2)
for nt=1:length(wt)-1
    v=real(V*exp(1i*wt(nt)));
    i=real(I*exp(1i*wt(nt)));

    plot(z,v,z,Z0*i)
    axis([-2 0 -2 2])
    Mov(nt)=getframe;
end

movie(Mov,5)
