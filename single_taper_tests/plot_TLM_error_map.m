tmp=load('G_lookup.mat');

lambda_mat=tmp.lambda;
beta_mat=tmp.beta;

N_lambda=numel(lambda_mat);
N_beta=numel(beta_mat);

E_pA=nan(N_lambda,N_beta);
E_qA=nan(N_lambda,N_beta);
E_pB=nan(N_lambda,N_beta);
E_qB=nan(N_lambda,N_beta);
ZcA_save=nan(N_lambda,N_beta);

dt_MOC_PA=nan(N_lambda,N_beta);
dt_TLM_PA=nan(N_lambda,N_beta);
dt_MOC_QA=nan(N_lambda,N_beta);
dt_TLM_QA=nan(N_lambda,N_beta);

gamma=nan(N_lambda,N_beta);

for i=1:N_lambda
    fprintf('%d/%d\n',i,N_lambda)
    for j=1:N_beta
        
        lambdaA=lambda_mat(i);
        beta=beta_mat(j);
        L=100;
        rA=12e-3;
        rB=lambdaA*rA;
        
        
        % Density
        rho = 870; %kg/m^3
        % Bulk Modulus
        %K = 1324.356^2*rho*2; %Pa
        K=1.5e9;
        %nu=1000e-6;
        c_bar=sqrt(K/rho);
        %beta = ((nu*L)/(c_bar*r1^2))*((lambda_1^2+lambda_1+1)^2/(9*lambda_1^3))
        nu = beta/((L/(c_bar*rA^2))*((lambdaA^2+lambdaA+1)^2/(9*lambdaA^3)));
        
        
        
        N_X=200;
              
        N_cycles=50;
        
        %% open outlet
        tic;
        [pA, qB, t, T, ZcA, ~, gamma(i,j), ~]=QA_unit_MOC(rA,rB,L, K, rho, nu,N_X,N_cycles);
        dt_MOC_QA(i,j)=toc;
        
        
        tic;
        simout=sim('QA_test');
        dt_TLM_QA(i,j)=toc;
        pA_sim=simout.yout{1}.Values.Data;
        qB_sim=-simout.yout{2}.Values.Data;
        t_sim=simout.tout;
        
        E_pA(i,j)=sqrt(mean((pA(:)-pA_sim(:)).^2));
        E_qB(i,j)=sqrt(mean((qB(:)-qB_sim(:)).^2));
        
        
        %% blocked outlet
        
        tic;
        [ pB, qA, t2]=PA_unit_MOC(rA,rB,L, K, rho, nu,N_X,N_cycles);
        dt_MOC_PA(i,j)=toc;
        
        
        tic;
        simout=sim('PA_test');
        dt_TLM_PA(i,j)=toc;
        qA_sim=simout.yout{1}.Values.Data;
        pB_sim=simout.yout{2}.Values.Data;
        t_sim2=simout.tout;
        
        E_qA(i,j)=sqrt(mean((qA(:)-qA_sim(:)).^2));
        E_pB(i,j)=sqrt(mean((pB(:)-pB_sim(:)).^2));
        
        ZcA_save(i,j)=ZcA;
    end
end

figure(1)
[c,h] = contour((lambda_mat-1)./(lambda_mat+1),log10(beta_mat),log10(E_qA'.*ZcA_save'),'k');
clabel(c,h)
xlabel('\phi')
ylabel('log10(\beta)')
title('log_{10}(E_{qA} Z_{cA})')
grid on


figure(2)
[c,h] = contour((lambda_mat-1)./(lambda_mat+1),log10(beta_mat),log10(E_pA'./ZcA_save'),'k');
clabel(c,h)
xlabel('\phi')
ylabel('log10(\beta)')
title('log_{10}(E_{pA}/Z_{cA})')
grid on

figure(3)
[c,h] = contour((lambda_mat-1)./(lambda_mat+1),log10(beta_mat),log10(E_qB'),'k');
clabel(c,h)
xlabel('\phi')
ylabel('log10(\beta)')
title('log_{10}(E_{qB})')
grid on


figure(4)
[c,h] = contour((lambda_mat-1)./(lambda_mat+1),log10(beta_mat),log10(E_pB'),'k');
clabel(c,h)
xlabel('\phi')
ylabel('log10(\beta)')
title('log_{10}(E_{pB})')
grid on

figure(5)
[c,h] = contour((lambda_mat-1)./(lambda_mat+1),log10(beta_mat),(dt_TLM_PA'),'k');
clabel(c,h)
xlabel('\phi')
ylabel('log10(\beta)')
title('dt_{PA}')
grid on

