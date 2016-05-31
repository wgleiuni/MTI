% Mechanical Topological Insulator, graphene and three-band material
% First to reproduce Srep18107

%% Parameters
eta=2/3;
kappa=1/(1-eta/2);
BS=1;
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
    GammaAB1=-kappa*(Gamma3+exp(1i*k*a1)*Gamma1+exp(1i*k*a2)*Gamma2);
    GammaAB2=-kappa*(Gamma3+exp(-1i*k*a1)*Gamma1+exp(-1i*k*a2)*Gamma2);
    diag=3*kappa*(1-eta/2)*eye(2);
    Gamma=[diag,GammaAB1;GammaAB2,diag];
    Ei(i,:)=eig(Gamma);
end

%% Figure

figure
if BS~=1
    set(gcf,'position',[2000,400,570,422],'color','w')
    plot(sqrt(abs(Ei)))
else
    set(gcf,'position',[2000,400,570,844],'color','w')
    subplot(2,1,1);
    plot(sqrt(abs(Ei)))
    title(['$\eta=$',num2str(eta),' $\kappa=$',num2str(1)],'interpreter','latex')
end

%% Ribbon band structure
if BS==1
    n=60;
    H=zeros((4*n+4)*2);
    H=H+3*(1-eta/2)*eye(length(H));
    kpara=linspace(-pi,pi,500);
    BEi=zeros(length(kpara),length(H));
    for i=1:length(kpara)
        tH=H;
        k=kpara(i);
        temp=-(Gamma1+Gamma2*exp(-1i*k));
        tempB=zeros(8);
        tempB(1:2,1:2)=-Gamma3;
        tempB(3:4,3:4)=-(Gamma2+Gamma1*exp(1i*k));
        tempB(5:6,5:6)=-Gamma3;
        tempB(7:8,7:8)=temp;
        tH(3:4,1:2)=temp;
        tH(1:2,3:4)=conj(temp);
        tH(end-1:end,end-3:end-2)=-(Gamma2+Gamma1*exp(1i*k));
        tH(end-3:end-2,end-1:end)=-(Gamma2+Gamma1*exp(-1i*k));
        tH(end-5:end-4,end-3:end-2)=-Gamma3;
        tH(end-3:end-2,end-5:end-4)=-Gamma3;
        tH(5:end-4,3:end-6)=tH(5:end-4,3:end-6)+kron(eye(n),tempB);
        tH(3:end-6,5:end-4)=tH(3:end-6,5:end-4)+kron(eye(n),conj(tempB));
        BEi(i,:)=eig(tH*kappa);
    end
    
    BEi(abs(imag(BEi))>1e-6 | real(BEi)<0)=NaN;
    BEi=sqrt(real(BEi));
    %tBEi=sqrt(real(BEi));
    %% Figure
    %figure
    subplot(2,1,2)
    
    plot(kpara/pi,BEi,'bo','Markersize',1)
    ylim([1.3 2.1]);
end