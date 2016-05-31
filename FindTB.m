% To find three-band parameters
tic
method=2;

if method==1
    eta1=1;
    eta2=1;
    kappa1=linspace(0,2,41);
    kappa1=kappa1(2:end);
    kappa2=linspace(0,2,41);
    kappa2=kappa2(2:end);
    
    result=zeros(length(kappa1)*length(kappa2),5);
    err=0.001;
    ks=[2/3*pi,2*sqrt(3)/3*pi];
    ke=[4*pi/3,0];
    
    idx=1;
    
    for k=1:length(kappa1)
        for l=1:length(kappa2)
            result(idx,1)=eta1;
            result(idx,2)=eta2;
            result(idx,3)=kappa1(k);
            result(idx,4)=kappa2(l);
            temp=fMTB(eta1,eta2,kappa1(k),kappa2(l),ks,ke);
            if (nnz(abs(temp(:,2)-temp(:,5))<err)>0)
                result(idx,3)=4;
            elseif (nnz(abs(temp(:,2)-temp(:,4))<err)>0)
                result(idx,3)=3;
            elseif (nnz(abs(temp(:,3)-temp(:,5))<err)>0)
                result(idx,3)=3;
            end
            idx=idx+1;
        end
        
    end
elseif method==2
    eta1=linspace(0,2,41);
    eta1=eta1(2:end);
    eta2=linspace(0,2,41);
    eta2=eta2(2:end);
    kappa1=linspace(0,2,41);
    kappa1=kappa1(2:end);
    kappa2=linspace(0,2,41);
    kappa2=kappa2(2:end);
    
    result=zeros(length(kappa1)*length(kappa2)*length(eta1)*length(eta2),5);
    err=0.001;
    ks=[2/3*pi,2*sqrt(3)/3*pi];
    ke=[4*pi/3,0];
    
    idx=1;
    for i=1:length(eta1)
        for j=1:length(eta2)
            for k=1:length(kappa1)
                for l=1:length(kappa2)
                    result(idx,1)=eta1(i);
                    result(idx,2)=eta2(j);
                    result(idx,3)=kappa1(k);
                    result(idx,4)=kappa2(l);
                    temp=fMTB(eta1(i),eta2(j),kappa1(k),kappa2(l),ks,ke);
                    if (nnz(abs(temp(:,2)-temp(:,5))<err)>0)
                        result(idx,3)=4;
                    elseif (nnz(abs(temp(:,2)-temp(:,4))<err)>0)
                        result(idx,3)=3;
                    elseif (nnz(abs(temp(:,3)-temp(:,5))<err)>0)
                        result(idx,3)=3;
                    end
                    idx=idx+1;
                end
            end
        end
        disp(i)
    end
end
%%
if (nnz(result(:,5))>0)
    disp(result(result(:,3)>0,:))
end
toc