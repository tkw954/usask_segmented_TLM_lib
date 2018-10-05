function [HA,HB,GA,GB,dt,T,beta,gamma,ZcA,ZcB] = calc_FIR_MOC_anechoic(rA,rB,L, K, rho, nu,N_X,N_cycles)
%CALC_FIR_MOC Summary of this function goes here
%   Detailed explanation goes here
if nargin<7
    N_X=200;
end

if nargin<8
    N_cycles=3.5;
end

if isempty(which('MOCinit'))
    error('MOC_solver  library must be added to path')
end


lambda_A=rB/rA;
lambda_B=rA/rB;

%% H A

[pA, qB, t_HA, T, ZcA, beta, gamma]=QA_unit_MOC_anechoic(rA,rB,L, K, rho, nu,N_X,N_cycles);
dt=t_HA(2)-t_HA(1);


Y_imp=diff(pA/ZcA);%impulse response

B_HA=Y_imp;
B_HA(1)=B_HA(1)-1;
A_HA=Y_imp;
A_HA(1)=A_HA(1)+1;

u_imp=[1 zeros(1,numel(Y_imp)-1)];

HA=filter(B_HA,A_HA,u_imp);


HA=filter([1/2 1/2],1,HA);%filter result to take into account only every second point is directly effected by impulse



%% H B

[pB, qA, ~, ~, ZcB]=QA_unit_MOC_anechoic(rB,rA,L, K, rho, nu,N_X,N_cycles);

Y_imp=diff(pB/ZcB);%impulse response

B_HB=Y_imp;
B_HB(1)=B_HB(1)-1;
A_HB=Y_imp;
A_HB(1)=A_HB(1)+1;

u_imp=[1 zeros(1,numel(Y_imp)-1)];

HB=filter(B_HB,A_HB,u_imp);


    HB=filter([1/2 1/2],1,HB);


%% G A
qBqA_imp=diff(qB);

%(1-H_A)
B_A=-HA;
B_A(1)=1-HA(1);


GA=1/lambda_A^2*filter(B_A,1,qBqA_imp);

GA=filter([1/2 1/2],1,GA);


%% G B
qAqB_imp=diff(qA);

%(1-H_B)
B_B=-HB;
B_B(1)=1-HB(1);


GB=1/lambda_B^2*filter(B_B,1,qAqB_imp);

GB=filter([1/2 1/2],1,GB);


end

