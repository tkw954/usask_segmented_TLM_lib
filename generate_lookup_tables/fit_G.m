function [mG_A,mG_B,n,...
    epsilon, t_G,GA_step,GB_step,G_hatA_step,G_hatB_step] ...
    = fit_G(lambda_A, beta)

%[ mE, mG, n ] = fit_G(lambda_A, beta)
% Inputs: lambda_A  = Taper ratio (dimensionless) r2/r1
%         beta    = dissipation number (dimensionless)
%
% Outputs: Optimized weighting factors
%

%
% Reference:
% T Wiens. Using a Segmented Transmission Line Model to 
% Simulate Laminar Pipeline Dynamics in Compound Tapered Tubes

if isempty(which('MOCinit'))
    error('MOC_solver  library must be added to path')
end

%assume some values, the optimized results are not sensitive to these
%values, only lambda_A and beta
rA = 20e-3; %(m) radius A
rB = rA*lambda_A; %(m) radius B
nu=100e-6;%(m^2/s) kinematic viscosity
K=1.5e9;%(Pa) bulk modulus
rho=890;%(kg/m^3) density
c=sqrt(K/rho);%(m/s) sonic speed

lambda_B=rA/rB;

% Assign length, L (m), based off of given dissipation number
L = (beta*c*rA^2/nu)*((9*lambda_A^3)/((lambda_A^2+lambda_A+1)^2));


N_X=200;%number of points in MOC solutions
N_cycles=6;

[HA,HB,GA,GB,dt,T,beta,gamma,ZcA,ZcB] = calc_FIR_MOC_anechoic(rA,rB,L, K, rho, nu,N_X,N_cycles);
N_tH=numel(HA);
N_tG=numel(GA);
t_H=dt*(0:(N_tH-1));
t_G=dt*(0:(N_tG-1))-T;

R_ZcA=24*beta/(1+lambda_A+lambda_A^2);
HA_DC=(R_ZcA+1/lambda_A^2-1)./(R_ZcA+1/lambda_A^2+1);

R_ZcB=24*beta/(1+lambda_B+lambda_B^2);
HB_DC=(R_ZcB+1/lambda_B^2-1)./(R_ZcB+1/lambda_B^2+1);

idx=2:(numel(GA)-1);



%idx_G=(t_H>=((1*T)-eps));
idx_G=find((t_G>0));
idx_G(1)=[];

GA_step=cumsum(GA);
GB_step=cumsum(GB);

GA=GA(idx_G);
GB=GB(idx_G);
GA_step=GA_step(idx_G);
GB_step=GB_step(idx_G);
t_G=t_G(idx_G);

%N_tG=numel(GA_step);
%t_G=dt*(0:(N_tG-1));

y=1-GA_step/(1-HB_DC);
%y=y(idx);
x=t_G/T;
%x=x(idx);

N_n=6;
%n=logspace(0,2,N_p);
n=1./linspace(0,1,N_n+1);
n(1)=[];
%phi=exp(-repmat(n',1,numel(x)).*repmat(x,N_n,1));

phi=exp(-repmat(n',1,numel(x)).*repmat(x,N_n,1));



mG_A=y/phi;
y_hat=mG_A*phi;

figure(100)
plot(x,[y;y_hat])
hold all
plot(xlim,0*[1 1],'k--')
hold off

G_hatA_step=(1-y_hat)*(1-HB_DC);


y=1-GB_step/(1-HA_DC);
%y=y(idx);

mG_B=y/phi;
y_hat=mG_B*phi;

G_hatB_step=(1-y_hat)*(1-HA_DC);



epsilon=sqrt(mean((GA_step-G_hatA_step).^2)+mean((GB_step-GB_step).^2));
