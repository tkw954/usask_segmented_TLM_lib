function [W1,W2,W3,n1,n2,n3,H_DC] = interpolate_H(lambda,beta)
%INTERPOLATE_H Summary of this function goes here
%   Detailed explanation goes here
lu_data=load('H_lookup.mat');




phi=(lambda-1)./(lambda+1);

if any(phi<min(lu_data.phi_A))||any(phi>max(lu_data.phi_A))
    warning('lambda out of range')
end

if any(beta<min(lu_data.beta))||any(beta>max(lu_data.beta))
    warning('beta out of range')
end

n1=lu_data.n1;
W1=nan(size(n1));

for i=1:numel(W1)
W1(i)=interp2(lu_data.phi_A,log10(lu_data.beta),squeeze(lu_data.W1(1,i,:,:)),phi,log10(beta),'spline');
end


idx=n1~=0;
H_DC1=2*(W1(~idx)+sum(W1(idx)./n1(idx).*(1-exp(-n1(idx)))));

n2=lu_data.n2;
W2=nan(size(n2));

for i=1:numel(W2)
W2(i)=interp2(lu_data.phi_A,log10(lu_data.beta),squeeze(lu_data.W2(1,i,:,:)),phi,log10(beta),'spline');
end

idx=n2~=0;
H_DC2=2*(W2(~idx)+sum(W2(idx)./n2(idx).*(1-exp(-n2(idx)))));

%n3=interp2(lu_data.phi_A,log10(lu_data.beta),squeeze(lu_data.W3(1,1,:,:)),phi,log10(beta),'spline');
n3=interp2(lu_data.phi_A,log10(lu_data.beta),lu_data.n3,phi,log10(beta),'spline');


%% Calculate part 3 based on DC values

R_Zc=24*beta/(1+lambda+lambda^2);
H_DC=(R_Zc+1/lambda^2-1)./(R_Zc+1/lambda^2+1);

H_DC3=H_DC-H_DC1-H_DC2;
W3=H_DC3*n3/2;

end

