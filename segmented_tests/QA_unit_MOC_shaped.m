function [pA, qB, t, T]=QA_unit_MOC_shaped(r_LU,x_LU, K, rho, nu,N_x,N_cycles)


if nargin<6
     N_x=200;
end

if nargin<7
    N_cycles=1;
end


if isempty(which('MOCinit'))
    error('MOC_solver  library must be added to path')
end

x_LU=x_LU-x_LU(1);
L=x_LU(end);


%% MOC params


N_t=round((N_x-1)*N_cycles)+1;

p_IC=0;
q_IC=0;

p_BC=[nan 0];

q_BC=repmat([1 nan],N_t);

r=@(x) interp1(x_LU,r_LU,x);

c=@(x) sqrt(K/rho)*ones(size(x));

%% friction
%steady only
% n=0;
% m=0;

%Trikha
% n=[26.4 200 8000];
% m=[1.0 8.1 40.0];

%johnston 2006
beta_f=2;
m1=1.4064;
m2=2.5200;
n1=33.104;
% r_approx=(r(0)+r(L))/2;
% dt_approx=1.4e-3;
%k=ceil((2*log(r_approx)-log(n1*nu*dt_approx))/(2*log(beta_f)))
k=4;

m=nan(k,1);
n=nan(k,1);
m(1)=m1;
m(2)=m2;
n(1)=n1;

for i=3:k
    m(i)=beta_f*m(i-1);
end

for i=2:k
    n(i)=beta_f^2*n(i-1);
end


%% solve
[ x,t,Zc,c_bar  ] = MOCinit( N_x,N_t, L, c, rho, r  );

T=(x_LU(end)-x_LU(1))/c_bar;


[ p, q ] =  MOCsolverF(x, t, p_IC, q_IC, p_BC, q_BC, Zc, r, nu, n, m  );

pA=p(:,1)';

qB=q(:,end)';

