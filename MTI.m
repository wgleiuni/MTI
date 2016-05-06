% Mechanical Topological Insulator, graphene and three-band material
% First to reproduce Srep18107

%% Parameters
eta=.1;
kappa=1/(1-eta/2);
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
Gamma1=getGamma(R1,eta);
Gamma2=getGamma(R2,eta);
Gamma3=getGamma(R3,eta);

%% Calculate band structure
kx1=linspace(0,1,2000);
ky1=linspace(0,sqrt(3),2000);
kx2=linspace(1,3/2,1000);
ky2=linspace(sqrt(3),sqrt(3)/2,1000);
kx3=linspace(3/2,0,2000);
ky3=linspace(sqrt(3)/2,0,2000);
Kx=[kx1,kx2,kx3];
Ky=[ky1,ky2,ky3];
K=[Kx',Ky']/(3/2/pi);

numK=length(Kx);
Ei=zeros(numK,4);

for i=1:numK
    k=K(i,:);
    GammaAB1=-kappa*(Gamma3+exp(-1i*k*a1)*Gamma1+exp(-1i*k*a2)*Gamma2);
    GammaAB2=-kappa*(Gamma3+exp(1i*k*a1)*Gamma1+exp(1i*k*a2)*Gamma2);
    diag=3*kappa*(1-eta/2)*eye(2);
    Gamma=[diag,GammaAB1;GammaAB2,diag];
    Ei(i,:)=eig(Gamma);
end

%% Figure
figure
plot(Ei)