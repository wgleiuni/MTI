% Mechanical Topological Insulator, graphene and three-band material
% Three band material

%% Parameters
eta1=0;
eta2=0;
kappa1=1;
kappa2=1;
BS=1;
m1=2;
m2=2;
m3=2;
R1=[sqrt(3);1];
R2=[-sqrt(3);1];
R3=[0;-2];
a1=[1;sqrt(3)];
a2=[-1;sqrt(3)];
dim=1;

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
if dim==1
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
    Ei=zeros(numK,6);
    
    for i=1:numK
        k=K(i,:);
        GammaAB=-kappa1*(Gamma1c+exp(1i*k*a1)*Gamma1a+exp(1i*k*a2)*Gamma1b);
        GammaBC=-kappa2*(Gamma2a+exp(-1i*k*(a1-a2))*Gamma2b+exp(-1i*k*a1)*Gamma2c);
        diag1=3*kappa1*(1-eta1/2)*eye(2);
        diag2=3*kappa2*(1-eta2/2)*eye(2);
        Gamma=[diag1/m1,GammaAB/m1,zeros(2);GammaAB'/m2,(diag1+diag2)/m2,GammaBC/m2;zeros(2),GammaBC'/m3,diag2/m3];
        Ei(i,:)=eig(Gamma);
    end
    
    %% Figure
    figure
    if BS~=1
        set(gcf,'position',[2000,400,570,422],'color','w')
        tEi=Ei;
        tEi(tEi<0)=NaN;
        plot(sqrt((tEi)))
        title(['$\eta_1=$',num2str(eta1),' $\eta_2=$',num2str(eta2),' $\kappa_1=$',num2str(kappa1),' $\kappa_2=$',num2str(kappa2)],'interpreter','latex')
    elseif BS==1
        set(gcf,'position',[2000,400,570,844],'color','w')
        subplot(2,1,1);
        tEi=Ei;
        tEi(tEi<0)=NaN;
        plot(sqrt((tEi)))
        title(['$\eta_1=$',num2str(eta1),' $\eta_2=$',num2str(eta2),' $\kappa_1=$',num2str(kappa1),' $\kappa_2=$',num2str(kappa2)],'interpreter','latex')
    end
elseif dim==2
    kx=linspace(-2*pi,2*pi,101);
    ky=kx;
    [Kx,Ky]=meshgrid(kx,ky);
    Ei=zeros(length(kx)*length(ky),6);
    for i=1:length(kx)*length(ky)
        k=[Kx(i),Ky(i)];
        GammaAB=-kappa1*(Gamma1c+exp(1i*k*a1)*Gamma1a+exp(1i*k*a2)*Gamma1b);
        GammaBC=-kappa2*(Gamma2a+exp(-1i*k*(a1-a2))*Gamma2b+exp(-1i*k*a1)*Gamma2c);
        diag1=3*kappa1*(1-eta1/2)*eye(2);
        diag2=3*kappa2*(1-eta2/2)*eye(2);
        Gamma=[diag1,GammaAB,zeros(2);GammaAB',diag1+diag2,GammaBC;zeros(2),GammaBC',diag2];
        Ei(i,:)=eig(Gamma);
    end
    tEi=Ei;
    tEi(tEi<0)=NaN;
    %% Figure
    figure
    set(gcf,'position',[2000,400,570,422],'color','w')
    hold on
    for i=1:6
        surf(Kx,Ky,reshape(sqrt(tEi(:,i)),length(kx),length(ky)),'linestyle','none')
    end
    view(45,45)
    title(['$\eta_1=$',num2str(eta1),' $\eta_2=$',num2str(eta2),' $\kappa_1=$',num2str(kappa1),' $\kappa_2=$',num2str(kappa2)],'interpreter','latex')
end

%% Ribbon band structure
if BS==1
    n=40;
    H=zeros((6*n+6)*2);
    temp=zeros(6);
    temp(1:2,1:2)=3*kappa2*(1-eta2/2)*eye(2);
    temp(3:4,3:4)=3*kappa1*(1-eta1/2)*eye(2)+3*kappa2*(1-eta2/2)*eye(2);
    temp(5:6,5:6)=3*kappa1*(1-eta1/2)*eye(2);
    H=H+kron(eye(2*(n+1)),temp);
    temp=zeros(2*3);
    temp(3:4,3:4)=-kappa2*Gamma2c;
    temp(5:6,5:6)=-kappa1*Gamma1c;
    temp=kron(eye(2*n+1),temp);
    H(5:end-2,1:end-6)=H(5:end-2,1:end-6)+temp;
    H(1:end-6,5:end-2)=H(1:end-6,5:end-2)+conj(temp);
    
    kpara=linspace(-pi,pi,500);
    BEi=zeros(length(kpara),length(H));
    for i=1:length(kpara)
        tH=H;
        k=kpara(i);
        T1=Gamma1a+exp(1i*k)*Gamma1b;
        T2=Gamma1b+Gamma1a*exp(-1i*k);
        Tb1=Gamma2a+exp(-1i*k)*Gamma2b;
        Tb2=exp(1i*k)*Gamma2a+Gamma2b;
        temp=zeros(2*6);
        temp(1:2,1:2)=-kappa2*Tb2;
        temp(3:4,3:4)=-kappa1*conj(T1);
        temp(5:6,5:6)=0*eye(2);
        temp(7:8,7:8)=-kappa2*conj(Tb1);
        temp(9:10,9:10)=-kappa1*conj(T2);
        temp=kron(eye(n+1),temp);
        temp=temp(1:end-2,1:end-2);
        tH(3:end,1:end-2)=tH(3:end,1:end-2)+temp;
        tH(1:end-2,3:end)=tH(1:end-2,3:end)+conj(temp);
        BEi(i,:)=eig(tH);
    end
    BEi(abs(imag(BEi))>1e-6 | real(BEi)<0)=NaN;
    BEi=sqrt(real(BEi));
    %% Figure
    %figure
    subplot(2,1,2)
    
    plot(kpara/pi,BEi,'bo','Markersize',1)
    %ylim([1.3 2.1]);
end