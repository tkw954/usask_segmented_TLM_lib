function [pA, qB, t, T, ZcA, beta, gamma, lambda_A]=QA_unit_MOC_anechoic(rA,rB,L, K, rho, nu,N_x,N_cycles)


if nargin<7
     N_x=200;
end

if nargin<8
    N_cycles=1;
end


if isempty(which('MOCinit'))
    error('MOC_solver  library must be added to path')
end



c_bar=sqrt(K/rho);% Local speed of sound (m/s) assuming constant
lambda_A = rB/rA; % Taper ratio (dimensionless), Equation (11)
%lambda_2 = r1/r2; % Taper ratio (dimensionless), Equation (11)
%beta = ((nu*L)/(c_bar*r1^2))*((lambda_1^2+lambda_1+1)^2/(9*lambda_1^3)); % Dissipation number (dimensionless), Equation (13)
ZcA=c_bar*rho/(pi*rA^2);
%Zc_2=c_bar*rho/(pi*rB^2);
T=L/c_bar; % Wave propagation time (s)
%nu=gamma/(((r1/abs(dr_dx))/(c_bar*r1^2))*((lambda_1^2+lambda_1+1)^2/(9*lambda_1^3)));


%    nu=gamma/((r1/dr_dx)/(c_bar*r1^2));

dr_dx=(rB-rA)/L;
gamma=nu*((rA/dr_dx)/(c_bar*rA^2));


beta = ((nu*L)/(c_bar*rA^2))*((lambda_A^2+lambda_A+1)^2/(9*lambda_A^3));

%% MOC params
N_t=round((N_x-1)*N_cycles)+1;

p_IC=0;
q_IC=0;

p_BC=[nan 0];

q_BC=repmat([1 nan],N_t,1);

r=@(x) rA+(rB-rA)/L*x;

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
r_approx=(r(0)+r(L))/2;
dt_approx=1.4e-3;
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
[ x,t,Zc ] = MOCinit( N_x,N_t, L, c, rho, r  );


tic
[ p, q, ] =  MOCsolverF_anechoic(x, t, p_IC, q_IC, p_BC, q_BC, Zc, r, nu, n, m  );
dt=toc;


%%check impuse

pA=p(:,1)';

qB=q(:,end)';



