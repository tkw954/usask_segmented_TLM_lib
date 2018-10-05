lambdaA=1;
%beta=1e-2;
L=100;


r_LU=12e-3*[1 0.5 0.5 1];
x_LU=L*[0 0.4 0.6 1];


% Density
rho = 870; %kg/m^3
% Bulk Modulus
%K = 1324.356^2*rho*2; %Pa
K=1.5e9;
nu=10e-6;
c_bar=sqrt(K/rho);
%beta = ((nu*L)/(c_bar*r1^2))*((lambda_1^2+lambda_1+1)^2/(9*lambda_1^3))
%nu = beta/((L/(c_bar*rA^2))*((lambdaA^2+lambdaA+1)^2/(9*lambdaA^3)));

T=L/c_bar;

N_X=200;



if isempty(which('MOCinit'))
    error('MOC_solver  library must be added to path')
end

N_cycles=12;

%% open outlet
tic;
[pA, qB, t, T]=QA_unit_MOC_shaped(r_LU,x_LU, K, rho, nu,N_X,N_cycles);
dt_MOC_QA=toc;


tic;
simout=sim('test_segmented_QA');
dt_TLM_QA=toc;
pA_sim=simout.yout{1}.Values.Data;
qB_sim=-simout.yout{2}.Values.Data;
t_sim=simout.tout;

E_pA=sqrt(mean((pA(:)-pA_sim(:)).^2));
E_qB=sqrt(mean((qB(:)-qB_sim(:)).^2));

%% closed outlet
tic;
[ pB, qA, t, T]=PA_unit_MOC_shaped(r_LU,x_LU, K, rho, nu,N_X,N_cycles);
dt_MOC_PA=toc;


tic;
simout=sim('test_segmented_PA');
dt_TLM_PA=toc;
qA_sim=simout.yout{1}.Values.Data;
pB_sim=simout.yout{2}.Values.Data;
t_sim=simout.tout;

E_qA=sqrt(mean((qA(:)-qA_sim(:)).^2));
E_pB=sqrt(mean((pB(:)-pB_sim(:)).^2));


Zc_A=rho*c_bar/(pi/4*r_LU(1)^2);

figure(1)
plot(t/T,pA/Zc_A)
hold all
plot(t_sim/T, pA_sim/Zc_A);
hold off
legend('MOC','TLM')
xlabel('t/T')
ylabel('p_A/Zc_A')

figure(2)
plot(t/T,qB)
hold all
plot(t_sim/T, qB_sim);
hold off
legend('MOC','TLM')
xlabel('t/T')
ylabel('q_B')

figure(3)
plot(t/T,qA*Zc_A)
hold all
plot(t_sim/T, qA_sim*Zc_A);
hold off
legend('MOC','TLM')
xlabel('t/T')
ylabel('q_A Z_{cA}')

figure(4)
plot(t/T,pB)
hold all
plot(t_sim/T, pB_sim);
hold off
legend('MOC','TLM')
xlabel('t/T')
ylabel('p_B')
