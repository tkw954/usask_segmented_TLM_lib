

tmp=load('H_lookup.mat');
beta=tmp.beta;
N_beta=numel(beta);
lambda_A=tmp.lambda_A;
lambda_A=lambda_A(lambda_A<=1);
lambda_A=sort(lambda_A,'descend');
N_lambda=numel(lambda_A);


tstart=tic;

k=6;
mG_A=nan(N_lambda,k,N_beta);
mG_B=nan(N_lambda,k,N_beta);
dt_elapsed=nan(N_lambda,N_beta);
epsilon=nan(N_lambda,N_beta);
% Loop through the range of dissipation number
for i=1:N_beta
    fprintf('beta %d/%d\n',i,N_beta)
    for j=1:N_lambda
        tic;
        [mG_A(j,:,i),mG_B(j,:,i),n,...
    epsilon(j,i), t_GA,GA_step,GB_step,G_hatA_step,G_hatB_step] ...
    = fit_G(lambda_A(j), beta(i));


      

    end
end

dt=toc(tstart);

gamma=repmat(beta,N_lambda,1)./repmat((sign(lambda_A-1).*(lambda_A-1).*(lambda_A.^2+lambda_A+1).^2./(9*lambda_A.^3))',1,N_beta);




% mE_lookup = permute(mE,[1 3 2]);
% mG_lookup = permute(mG,[1 3 2]);
% tau_lookup = tau;

idx=1;
figure(1)
set(gcf,'defaultaxescolororder',jet(N_beta))
plot((lambda_A-1)./(lambda_A+1),squeeze(mG_A(:,idx,:)))
hold all
plot((1./lambda_A-1)./(1./lambda_A+1),squeeze(mG_B(:,idx,:)))
hold off
xlabel('(\lambda-1)/(\lambda+1)=(r_2-r_1)/(r_2+r_1)')
ylabel(sprintf('param(%d)',idx(1)))
grid on



figure(2)
[cc,hc]=contour(lambda_A,log10(beta),log10(epsilon'));
clabel(cc,hc);
hold all
[cc,hc]=contour(lambda_A,log10(beta),log10(gamma'),'k:');
clabel(cc,hc);
caxis([min(log10(epsilon(:))) max(log10(epsilon(:)))])
hold off
xlabel('\lambda_1')
ylabel('log_{10}(\beta)')

figure(3)
Nrow=ceil(sqrt(2*k));
Ncol=ceil(2*k/Nrow);
for i=1:k
    subplot(Nrow,Ncol,i)
    pcolor(lambda_A,log10(beta),squeeze(mG_A(:,i,:))')
    title(sprintf('mG_A(%d)',i))
    colorbar
    shading interp
        subplot(Nrow,Ncol,i+k)
    pcolor(lambda_A,log10(beta),squeeze(mG_B(:,i,:))')
    title(sprintf('mG_B(%d)',i))
    colorbar
    shading interp
end

figure(4)
set(gcf,'defaultaxescolororder',jet(N_beta))
Nrow=ceil(sqrt(k));
Ncol=ceil(k/Nrow);
for i=1:k
    subplot(Nrow,Ncol,i)

    
    plot((lambda_A-1)./(lambda_A+1),squeeze(mG_A(:,i,:)))
    hold all
    plot((1./lambda_A-1)./(1./lambda_A+1),squeeze(mG_B(:,i,:)))
    hold off
    xlabel('(\lambda-1)/(\lambda+1)=(r_2-r_1)/(r_2+r_1)')
    ylabel(sprintf('mG(i)',i))
    grid on
end


%% save results (this only works for N_params=12)
lambda=[fliplr(lambda_A) 1./lambda_A(2:end)];
phi=(lambda-1)./(lambda+1);

mG=[flipud(mG_A); mG_B(2:end,:,:)];

fname=['G_lookup_' datestr(now,30) '.mat'];
save(fname,'beta','lambda','phi','mG','n','dt')