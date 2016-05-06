function [out]=getGamma(R,eta)
R=R/norm(R);
out=zeros(2);
out(1,1)=(1-eta)+eta*R(1)*R(1);
out(1,2)=eta*R(1)*R(2);
out(2,1)=eta*R(2)*R(1);
out(2,2)=(1-eta)+eta*R(2)*R(2);
end