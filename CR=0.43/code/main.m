gamma0=-0.5;
gamma1=2;
gamma2=1;
gamma3=1;
beta1=1;
beta2=0.2;
beta3=0.5;
B=zeros(6,2);
B(1:3,1)=0.8;
B(4:6,2)=0.8;
BB=1;
n=1000;
for i=1:BB
    
    [eta,Z,X,delta,V]=data2(beta1,beta2,beta3,gamma0,gamma1,gamma2,gamma3,n,B);
    
    id=(1:n)';
    V1=V(:,1);
    V2=V(:,2);
    V3=V(:,3);
    V4=V(:,4);
    V5=V(:,5);
    V6=V(:,6);
    columns ={'id','delta','X','V1','V2','V3','V4','V5','V6','Z'};
    data=table(id,delta,X,V1,V2,V3,V4,V5,V6,Z,'VariableNames',columns);
    filename=sprintf('Coxlatent%d.csv',i);
    writetable(data,filename);
    
   
end