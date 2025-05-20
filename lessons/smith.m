% smith chart
close all; clear all
%YL=3.33e-3;
% YL=3.2e-3-2.4e-3i
% ZL=1/YL;
%Z0=75;
Z0=50;
ZL=30;
GammaL=(ZL-Z0)/(ZL+Z0);
f=0.5e9;
%epsr=1;
%epsr=3;% 2.4169;
Pinc=250e-3; % Pinc=10^(30/10)*1e-3; 
c=299792458;
eta0=376.7303;% d=1e-3; D=d*cosh(pi*sqrt(epsr)*Z0/eta0);
%vp=c/sqrt(epsr);
vp=2e8;
%k=2*pi*f/c*sqrt(epsr);
k=2*pi*f/vp;
lambda=2*pi/k;
lambda_mm=lambda*1e3
xx=linspace(0,2*pi,1001);
figure(1)
plot(exp(1i*xx),'k')
hold on
axis equal tight
plot(1/2+exp(1i*xx)/2,'k')
plot(1+1i-exp(1i*xx/4),'k')
plot(1-1i-exp(-1i*xx/4),'k')
plot([-1 1],[0 0],'k')
plot(real(GammaL),imag(GammaL),'or')
plot(real(GammaL),imag(GammaL),'.r')
text(real(GammaL),imag(GammaL),'Z_L')
VSWR=(1+abs(GammaL))/(1-abs(GammaL))
PL=(1-abs(GammaL)^2+2i*imag(GammaL))*Pinc

zL=ZL/Z0;
rL=real(zL);
xL=imag(zL);
yL=1/zL;
gL=real(yL);
bL=imag(yL);

% %% Lumped elements L-network
% % type A     --+---[jx]--+
% %             [jb]      [zL]
% %            --+---------+
% if rL<1
%     xA1=-xL+sqrt(rL*(1-rL));
%     xA2=-xL-sqrt(rL*(1-rL));
%     
%     bA1=sqrt(1/rL-1);
%     bA2=-sqrt(1/rL-1);
%     
%     if xA1>0
%         LxA1=Z0*xA1/(2*pi*f)
%     else
%         CxA1=-1/(2*pi*f*Z0*xA1)
%     end
%     if bA1>0
%         CbA1=bA1/(2*pi*f*Z0)
%     else
%         LbA1=-Z0/(2*pi*f*bA1)
%     end
%     
%    if xA2>0
%         LxA2=Z0*xA2/(2*pi*f)
%     else
%         CxA2=-1/(2*pi*f*Z0*xA2)
%     end
%     if bA2>0
%         CbA2=bA2/(2*pi*f*Z0)
%     else
%         LbA2=-Z0/(2*pi*f*bA2)
%     end
% end
% 
% 
% %% type B
% %            --[jx]--+-----+
% %                   [jb]  [zL]
% %            --------+-----+
% if gL<1
%     bB1=-bL+sqrt(gL*(1-gL));
%     bB2=-bL-sqrt(gL*(1-gL));
%     
%     xB1=sqrt(1/gL-1);
%     xB2=-sqrt(1/gL-1);
%     
%     if xB1>0
%         LxB1=Z0*xB1/(2*pi*f)
%     else
%         CxB1=-1/(2*pi*f*Z0*xB1)
%     end
%     if bB1>0
%         CbB1=bB1/(2*pi*f*Z0)
%     else
%         LbB1=-Z0/(2*pi*f*bB1)
%     end
%     if xB2>0
%         LxB2=Z0*xB2/(2*pi*f)
%     else
%         CxB2=-1/(2*pi*f*Z0*xB2)
%     end
%     if bB2>0
%         CbB2=bB2/(2*pi*f*Z0)
%     else
%         LbB2=-Z0/(2*pi*f*bB2)
%     end
%     
% end
% ff=linspace(0.999*f,1.001*f,1001);
% jomega=2i*pi*ff;
% 
% % ZZL=Z0*(rL+jomega*xL/(2*pi*f)); YYL=1./ZZL; % carico serie R+jwL
% YYL=1/Z0*(gL-bL*(2*pi*f)./(jomega)); ZZL=1./YYL; % carico parallelo G+1/jwL
% 
%  ZinA1=1./(jomega*CbA1+1./(1./(jomega*CxA1)+ZZL));
%  ZinA2=1./(1./(jomega*LbA2)+1./(1./(jomega*CxA2)+ZZL));
%  
%  ZinB1=jomega*LxB1+1./(jomega*CbB1+YYL);
%  ZinB2=1./(jomega*CxB2)+1./(jomega*CbB2+YYL);
% 
% GammaA1=(ZinA1-Z0)./(ZinA1+Z0);
% GammaA2=(ZinA2-Z0)./(ZinA2+Z0);
% GammaB1=(ZinB1-Z0)./(ZinB1+Z0);
% GammaB2=(ZinB2-Z0)./(ZinB2+Z0);
% 
% 
% 
% dGdwA1=(GammaA1(end)-GammaA1(1))/(ff(end)-ff(1))*f;
% abs(dGdwA1)
% %1/2*abs(abs(bA1)-(1-1i*bA1).^2*(abs(xA1)+abs(xL)))
% 1/2*abs(abs(bA1)-(1-1i*bA1).^2*(abs(xA1)-abs(bL)/yL^2))
% 
% dGdwA2=(GammaA2(end)-GammaA2(1))/(ff(end)-ff(1))*f;
% abs(dGdwA2)
% %1/2*abs(abs(bA2)-(1-1i*bA2).^2*(abs(xA2)+abs(xL)))
% 1/2*abs(abs(bA2)-(1-1i*bA2).^2*(abs(xA2)-abs(bL)/yL^2))
% 
% dGdwB1=(GammaB1(end)-GammaB1(1))/(ff(end)-ff(1))*f;
% abs(dGdwB1)
% %1/2*abs(abs(xB1)-(1-1i*xB1).^2*(abs(bB1)-abs(xL)/zL^2))
% 1/2*abs(abs(xB1)-(1-1i*xB1).^2*(abs(bB1)+abs(bL)))
% 
% dGdwB2=(GammaB2(end)-GammaB2(1))/(ff(end)-ff(1))*f;
% abs(dGdwB2)
% %1/2*abs(abs(xB2)-(1-1i*xB2).^2*(abs(bB2)-abs(xL)/zL^2))
% 1/2*abs(abs(xB2)-(1-1i*xB2).^2*(abs(bB2)+abs(bL)))
% 
% figure(10)
% plot(ff,abs(GammaA1),ff,abs(GammaA2),ff,abs(GammaB1),ff,abs(GammaB2))
% legend('\Gamma_{A1}','\Gamma_{A2}','\Gamma_{B1}','\Gamma_{B2}')
% 
% figure(11)
% plot(real(GammaA1),imag(GammaA1),real(GammaA2),imag(GammaA2),real(GammaB1),imag(GammaB1),real(GammaB2),imag(GammaB2))
% legend('\Gamma_{A1}','\Gamma_{A2}','\Gamma_{B1}','\Gamma_{B2}')
% 

ll=linspace(0,lambda/4,10000);
Gamma=GammaL*exp(-2i*k*ll);
z=(1+Gamma)./(1-Gamma);
y=1./z;
% q=find(abs(imag(y))==min(abs(imag(y)))); % z reale per trasf lambda/4
q=find(abs(1-real(y))==min(abs(1-real(y)))); % y=1+jb stub parallelo
% q=find(abs(1-real(z))==min(abs(1-real(z)))); % z=1+jx stub serie

plot(Gamma(1:q),'b','linewidth',2)
plot(Gamma(q),'og')
plot(Gamma(q),'.g')
text(real(Gamma(q)),imag(Gamma(q)),'Z_1')

line_length_lambda=ll(q)/lambda
line_length_mm=ll(q)*1e3
text(0,0,'d')
z1=z(q)
Z1=z(q)*Z0

%% stub parallelo 
y1=1/z1;
ys=-1i*imag(y1);
zs=1/ys;
plot(-1/2-exp(1i*xx)/2,':k') % g=1
Gbcst=(1-(ys+linspace(0,100,1000)))./(1+(ys+linspace(0,100,1000))); % b=bstub
plot(Gbcst,':k')
plot(conj(Gbcst),':k')
text(0,0,'g=1')
text(0,0,'b=b_1')
text(0,0,'b=b_s')

%% stub serie
% zs=-1i*imag(z1);
% Zs=zs*Z0
% 
% Gxcst=((zs+linspace(0,100,1000))-1)./((zs+linspace(0,100,1000))+1); % x=xstub
% plot(Gxcst,':k')
% plot(conj(Gxcst),':k')

Gammas=(zs-1)/(zs+1);

%phs=linspace(pi,angle(Gammas)); % sc stub
 phs=linspace(0,angle(Gammas)); % oc stub
lstub_mm=(phs(1)-phs(end))*lambda/(4*pi)*1e3
lstub_lambda=(phs(1)-phs(end))/(4*pi)
% 
plot(1,0,'og'), plot(1,0,'.g'), text(1,0,'o.c.') 
% plot(-1,0,'og'), plot(-1,0,'.g'), text(-1,0,'s.c.') 
% % 
plot(exp(1i*phs),'r','linewidth',2)
plot(Gammas,'ob')
plot(Gammas,'.b')
text(real(Gammas),imag(Gammas),'Z_s')
text(0,0,'d')
text(0,0,'l_s')
 
%% trasformatore lambda/4
% Zt=sqrt(z(q))*Z0;
% ll=linspace(0,lambda/4,1000);
% text(0,0,'l_t=\lambda/4')
% text(0,0,'l_t=\lambda/4')
% Gamma1I=(Z1-Zt)/(Z1+Zt)*exp(-2i*k*ll);
% 
% plot(Gamma1I,':','linewidth',2);
% plot(Gamma1I(1),'og')
% plot(Gamma1I(1),'.g')
% text(real(Gamma1I(1)),imag(Gamma1I(1)),'Z^\prime_1')
% plot(Gamma1I(end),'oc')
% plot(Gamma1I(end),'.c')
% text(real(Gamma1I(end)),imag(Gamma1I(end)),'Z^\prime_2')
% 
% Z2=(1+Gamma1I(end))./(1-Gamma1I(end))*Zt;
% Gamma2I=(Z2-Z0)/(Z2+Z0);
% plot(Gamma2I,'oc')
% plot(Gamma2I,'.c')
% text(real(Gamma2I),imag(Gamma2I),'Z_2=Z_0')
% 
% 
% zz1=(1+Gamma1I)./(1-Gamma1I);
% ZZ1=Zt*zz1;
% zz1=ZZ1/Z0;
% Gamma1=(ZZ1-Z0)./(ZZ1+Z0);
% plot(Gamma1,'linewidth',2);

axis off

PL=Pinc*(1+1i*imag(zL)/real(zL))
Ps=zs*Pinc % potenza reattiva stub serie
%Ps=conj(ys)*Pinc % potenza reattiva stub parallelo
P1=z1*Pinc
%Z2=(1+Gamma1(end))/(1-Gamma1(end))*Zt;
% Gammam=1i*(Z2-Zt)/(Z2+Zt);
% Pm=Pinc*(1+2i*imag(Gammam)/(1-abs(Gammam)^2))
