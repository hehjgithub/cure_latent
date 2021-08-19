%macro nlm(data);
proc nlmixed data=&data qpoints=10 hess;
parms 
b11=0.1  b21=0.1  b31=0.1   b42=0.1  b52=0.1  b62=0.1
sigma1=0.1 sigma2=0.1 sigma3=0.1 sigma4=0.1 sigma5=0.1 sigma6=0.1
sigma12=0.1 beta1=0.1 beta2=0.1 beta3=0.1 gamma0=0.1 gamma1=1.1 gamma2=0.1 gamma3=0.1
h1=0.5 h2=0.5 h3=0.5 h4=0.5 h5=0.5 ;
bounds sigma1 sigma2 sigma3 sigma4 sigma5 sigma6   h1 h2 h3 h4 h5 >0;

temp1=-1/2*(V1-b11*nu1)**2/sigma1-1/2*(V2-b21*nu1)**2/sigma2
      -1/2*(V3-b31*nu1)**2/sigma3-1/2*(V4-b42*nu2)**2/sigma4
	  -1/2*(V5-b52*nu2)**2/sigma5-1/2*(V6-b62*nu2)**2/sigma6
	  -1/2*(log(sigma1)+log(sigma2)+log(sigma3)+log(sigma4)+log(sigma5)+log(sigma6));
loglik1=temp1;
temp2=gamma0+gamma1*Z+gamma2*nu1+gamma3*nu2;
temp3=beta1*Z+beta2*nu1+beta3*nu2;
temp4=h1*L1+h2*L2+h3*L3+h4*L4+h5*L5;
if delta=0 then do;
temp5=1/(1+exp(temp2))+exp(temp2)/(1+exp(temp2))*exp(-temp4*exp(temp3));
loglik=loglik1+log(temp5);
end;
if delta=1 then do;
temp6=p1*h1+p2*h2+p3*h3+p4*h4+p5*h5;
temp7=exp(temp2)/(1+exp(temp2))*exp(-temp4*exp(temp3))*temp6*exp(temp3);
loglik=loglik1+log(temp7);
end;



model X~ general(loglik);
random nu1 nu2 ~normal([0,0],[1,sigma12,1]) subject=id;
ods output ParameterEstimates=two.test1 hessian=two.cov;
run;
%mend();
%macro estimate();
libname two 'E:\article\AH-cure-latent\simulation\code\CR=0.43\temp result';
%do jj=1%to 1;
proc import out=two.Coxlatent&jj 
datafile="E:\article\AH-cure-latent\simulation\code\CR=0.43\data\Coxlatent&jj..csv" 
dbms=csv
replace;
run;


data two.one;
set two.Coxlatent&jj;
if delta=1;
run;

proc capability data=two.one noprint;
var X; 
output out=two.quant_death pctlpts=0 20 40  60  80 100 pctlpre=qd; 
run;

data two.quant_death;
set two.quant_death;
aa=1;
run;


data two.Coxlatent&jj;
set two.Coxlatent&jj;
aa=1;
run;

data two.Coxlatent&jj;
merge two.Coxlatent&jj two.quant_death;
by aa;
run;

data two.Coxlatent&jj;
set two.Coxlatent&jj;
array p {5} p1-p5;
array L {5} L1-L5;
array q{6} qd0 qd20 qd40 qd60 qd80 qd100 ;

do ii=1 to 5;
	p{ii}=0;
	L{ii}=0;
end;

if X>=q{1} then do;
P{1}=1;
end;
do ii=2 to 5;
if X>=q{ii} then do;
p{ii-1}=0;
p{ii}=1;
end;
end;
if X>q{6} then do;
p{5}=0;
end;

do ii=1 to 5;
if X>=q{ii} then do;
if X>=q{ii+1} then do;
L{ii}=q{ii+1}-q{ii};
end;
else do;
L{ii}=X-q{ii};
end;
end;
end;
run;






%nlm(two.Coxlatent&jj) 

data two.one&jj;
set _null_;
run;
data two.one&jj;
set two.one&jj two.test1;
KEEP parameter Estimate;
run;

data two.cov&&jj;
set  two.cov;
run;




%end;
%mend();
%estimate()

%macro toge();
data two.estall;  
set _null_;
run;

data two.hesall;  
set _null_;
run;


%do ii=1%to 1;
data two.estall;   
set two.estall two.one&ii;   
run;

data two.hesall;  
set two.hesall two.cov&ii;  
run;

%end;

proc export data=two.estall
outfile="E:\article\AH-cure-latent\simulation\code\CR=0.43\result\estall.csv"
dbms=csv
label
replace;
run;

proc export data=two.hesall
outfile="E:\article\AH-cure-latent\simulation\code\CR=0.43\result\hesall.csv"
dbms=csv
label
replace;
run;

%mend();
%toge()
