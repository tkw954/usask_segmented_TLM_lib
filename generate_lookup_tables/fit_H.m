function [W1, W2, n1, n2, n3, E, x,y,x_hat,y_hat]=fit_H(lambdaA,beta,N_n1,N_n2,N_X)


if nargin<3
    N_n1=6;
end

if nargin<4
    N_n2=6;
end

if nargin<5
    N_X=200;
end



L=1000;
rA=12e-3;
rB=lambdaA*rA;


% Density
rho = 870; %kg/m^3
% Bulk Modulus
K=1.5e9;
c_bar=sqrt(K/rho);
nu = beta/((L/(c_bar*rA^2))*((lambdaA^2+lambdaA+1)^2/(9*lambdaA^3)));




if isempty(which('MOCinit'))
    error('MOC_solver  library must be added to path')
end

N_cycles=10;
%% H A

[pA, ~, t_HA, T, ZcA]=QA_unit_MOC_anechoic(rA,rB,L, K, rho, nu,N_X,N_cycles);
dt=t_HA(2)-t_HA(1);


Y_imp=diff(pA/ZcA);%impulse response

B_HA=Y_imp;
B_HA(1)=B_HA(1)-1;
A_HA=Y_imp;
A_HA(1)=A_HA(1)+1;

u_imp=[1 zeros(1,numel(Y_imp)-1)];

HA=filter(B_HA,A_HA,u_imp);

HA=filter([1/2 1/2],1,HA);

t=(0:(numel(HA)-1))*dt;
%% part 1

idx1=find((t>(0*T))&(t<(2*T)));
k_skip1=1;
idx1(1:k_skip1)=[];
k_endskip1=1;
idx1((end-k_endskip1+1):end)=[];
HA1=HA(idx1);
t1=t(idx1);


y1=HA1/dt*T;

x1=t1/T/2;


n1=1./linspace(0,1*0.3,N_n1);
n1(1)=0;
phi1=exp(-repmat(n1',1,numel(x1)).*repmat(x1,N_n1,1));

W1=y1/phi1;
y_hat1=W1*phi1;





E1=sqrt(mean((y_hat1-y1).^2));

E_rel1=E1/(max(y1)-min(y1));

HA_DC1=2*(W1(1)+sum(W1(2:end)./n1(2:end).*(1-exp(-n1(2:end)))));

%% part 2

idx2=find((t>(2*T))&(t<(4*T)));
k_skip2=k_skip1;
idx2(1:k_skip2)=[];
k_endskip2=k_endskip1;
idx2((end-k_endskip2+1):end)=[];
HA2=HA(idx2);
t2=t(idx2);

y2=HA2/dt*T;

x2=(t2-2*T)/T/2;



n2=1./linspace(0,1*0.3,N_n2);
n2(1)=0;
phi2=exp(-repmat(n2',1,numel(x2)).*repmat(x2,N_n2,1));
W2=y2/phi2;
y_hat2=W2*phi2;



E2=sqrt(mean((y_hat2-y2).^2));
E_rel2=E2/(max(y2)-min(y2));

HA_DC2=2*(W2(1)+sum(W2(2:end)./n2(2:end).*(1-exp(-n2(2:end)))));



%% Part 3
idx3=find(t>(4*T));
k_skip3=k_skip1;
idx3(1:k_skip3)=[];
k_endskip3=k_endskip1;
idx3((end-k_endskip3+1):end)=[];
HA3=HA(idx3);
t3=t(idx3);

y3=HA3/dt*T;

x3=(t3-4*T)/T/2;

Y=log(abs(y3));

X=[x3;ones(size(x3))];

W3_tmp=Y/X;

W3=[ W3_tmp(1) sign(y3(1))*exp(W3_tmp(2))];

y_hat3=W3(2)*exp(W3(1)*x3);



E3=sqrt(mean((y_hat3-y3).^2));
E_rel3=E3/(max(y3)-min(y3));

n3=-W3(1);

E=[E1 E2 E3;E_rel1 E_rel2 E_rel3];


x=t/T/2;
y=HA/dt*T;
 
x_hat=[x1 x2+1 x3+2];
y_hat=[y_hat1 y_hat2 y_hat3];



