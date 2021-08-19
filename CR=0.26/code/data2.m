function[Z,X,delta,V,eta]=data2(beta1,beta2,beta3,gamma0,gamma1,gamma2,gamma3,n,B)
%lambda_0(t)=5+t;
Z=binornd(1,0.5,n,1);
Phi=[1,0.2;0.2,1];
nu_ksi=[0;0];
ksi=mvnrnd(nu_ksi,Phi,n);
P=exp(gamma0+gamma1.*Z+gamma2.*ksi(:,1)+gamma3.*ksi(:,2))./(1+exp(gamma0+gamma1.*Z+gamma2.*ksi(:,1)+gamma3.*ksi(:,2)));
eta=binornd(1,P);
C=unifrnd(0,1,n,1);
U=unifrnd(0,1,n,1);
temp1=-log(U)./exp(beta1.*Z+beta2.*ksi(:,1)+beta3.*ksi(:,2));
T=(-10+sqrt(100+8*temp1))./2;
X=(1-eta).*C+eta.*(min(T,C));
delta=eta.*(T<=C);
nu_epsilon=zeros(6,1);
Phi_epsilon=eye(6).*0.3;
epsilon=mvnrnd(nu_epsilon,Phi_epsilon,n);
V=zeros(n,6);
for i=1:n
    V(i,:)=(B*(ksi(i,:))'+(epsilon(i,:))')';
end