% Range of dissipation number (beta) (dimensionless)
N_beta=60;
beta = logspace(-5,0,N_beta);

% Range of taper ratio (lambda) (dimensionless)
N_lambda=51;%must be odd
%lambda_A = linspace(1,0.1,N_lambda);
phi_A=linspace(-0.5,0.5,N_lambda);
lambda_A=(1+phi_A)./(1-phi_A);

N_n1=6;
N_n2=6;

tstart=tic;
for i=1:N_beta
    fprintf('%d/%d\n',i,N_beta)
    for j=1:N_lambda
        [W1(:,:,i,j), W2(:,:,i,j), n1(:,:,i,j), n2(:,:,i,j), n3(i,j), E(:,:,i,j)]=fit_H(lambda_A(j),beta(i),N_n1,N_n2);
    end
end
dt=toc(tstart);

%n1 and n2 don't change
n1=n1(:,:,1,1);
n2=n2(:,:,1,1);

figure(1)

[cc,h] = contour(phi_A,log10(beta),squeeze(log10(E(2,1,:,:))));
clabel(cc,h)
xlabel('\phi')
ylabel('log_{10}(\beta)')
title('Relative Error Part 1')

figure(2)

[cc,h] = contour(phi_A,log10(beta),squeeze(log10(E(2,2,:,:))));
clabel(cc,h)
xlabel('\phi')
ylabel('log_{10}(\beta)')
title('Relative Error Part 2')


figure(3)
set(gcf,'defaultAxesColorOrder',jet(N_beta))
colormap(jet)
plot(phi_A,squeeze(W1(1,1,:,:)))
xlabel('\phi')
ylabel('W1(1)')
caxis([log10(beta(1)) log10(beta(end))])
y=colorbar;
ylabel(y,'log_{10}(\beta)')


figure(4)
set(gcf,'defaultAxesColorOrder',jet(N_lambda))
colormap(jet)
plot(log10(beta),squeeze(W1(1,1,:,:))')
xlabel('log_{10}(\beta)')
ylabel('W1(1)')
caxis([phi_A(1) phi_A(end)])
y=colorbar;
ylabel(y,'\phi')

figure(5)
set(gcf,'defaultAxesColorOrder',jet(N_beta))
colormap(jet)
N_r=ceil(sqrt(N_n1));
N_c=ceil(N_n1/N_r);
for i=1:N_n1
    subplot(N_r,N_c,i)
plot(phi_A,squeeze(W1(1,i,:,:)))
xlabel('\phi')
ylabel(sprintf('W1(%d)',i))
end

figure(6)
set(gcf,'defaultAxesColorOrder',jet(N_beta))
colormap(jet)
N_r=ceil(sqrt(N_n2));
N_c=ceil(N_n2/N_r);
for i=1:N_n2
    subplot(N_r,N_c,i)
plot(log10(beta),squeeze(W2(1,i,:,:))')
xlabel('log_{10}(\beta)')
ylabel(sprintf('W2(%d)',i))
end

figure(7)
set(gcf,'defaultAxesColorOrder',jet(N_beta))
colormap(jet)


plot(log10(beta),n3')
xlabel('log_{10}(\beta)')
ylabel('n_3')


%% save results
fname=['H_lookup_' datestr(now,30) '.mat'];
save(fname,'beta','lambda_A','phi_A','W1','W2','n1','n2','n3','dt')