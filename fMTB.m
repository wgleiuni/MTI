function [Ei]=fMTB(eta1,eta2,kappa1,kappa2,ks,ke)
% Function version of MTB, used to find whether there are three band
% situation

%% K-space line
kstep=100;
kx=linspace(ks(1),ke(1),kstep)';
ky=linspace(ks(2),ke(2),kstep)';
K=[kx,ky];
%% Parameters
R1=[sqrt(3);1];
R2=[-sqrt(3);1];
R3=[0;-2];
a1=[1;sqrt(3)];
a2=[-1;sqrt(3)];

%% Normalize vectors R1, R2, R3
R1=R1/norm(R1);
R2=R2/norm(R2);
R3=R3/norm(R3);
a1=a1/norm(a1);
a2=a2/norm(a2);

%% Gamma matrix
Gamma1a=getGamma(R1,eta1);
Gamma1b=getGamma(R2,eta1);
Gamma1c=getGamma(R3,eta1);

Gamma2a=getGamma(R1,eta2);
Gamma2b=getGamma(R2,eta2);
Gamma2c=getGamma(R3,eta2);

%% Calculate band structure
% dim=1 for 1D figure, dim=2 for 2D BZ

% kx1=linspace(0,1,2*kstep);
% ky1=linspace(0,sqrt(3),2*kstep);
% kx2=linspace(1,3/2,kstep);
% ky2=linspace(sqrt(3),sqrt(3)/2,kstep);
% kx3=linspace(3/2,0,2*kstep);
% ky3=linspace(sqrt(3)/2,0,2*kstep);
% Kx=[kx1,kx2,kx3];
% Ky=[ky1,ky2,ky3];
% K=[Kx',Ky']/(3/2/pi);

numK=length(K);
Ei=zeros(numK,6);

for i=1:length(K)
    k=K(i,:);
    GammaAB=-kappa1*(Gamma1c+exp(1i*k*a1)*Gamma1a+exp(1i*k*a2)*Gamma1b);
    GammaBC=-kappa2*(Gamma2a+exp(-1i*k*(a1-a2))*Gamma2b+exp(-1i*k*a1)*Gamma2c);
    diag1=3*kappa1*(1-eta1/2)*eye(2);
    diag2=3*kappa2*(1-eta2/2)*eye(2);
    Gamma=[diag1,GammaAB,zeros(2);GammaAB',diag1+diag2,GammaBC;zeros(2),GammaBC',diag2];
    Ei(i,:)=sort(eig(Gamma));
end
end