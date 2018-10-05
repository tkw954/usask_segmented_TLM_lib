function [mG, n] = interpolate_G(lambda,beta)
%INTERPOLATE_H Summary of this function goes here
%   Detailed explanation goes here
lu_data=load('G_lookup.mat');

phi=(lambda-1)./(lambda+1);

if any(phi<min(lu_data.phi))||any(phi>max(lu_data.phi))
    warning('G:lambda out of range')
end

if any(beta<min(lu_data.beta))||any(beta>max(lu_data.beta))
    warning('G: beta out of range')
end

% k=size(lu_data.mG,2); %number of weighting function terms
% n=nan(1,k);%weighting function coefficient
% n(1)=0.3/(1+3*beta);% Equation (19)
% for i=2:k
%     n(i)=n(i-1)*3;% Equation (19)
% end

n=lu_data.n;
k=numel(n);

mG=nan(1,k);
for i=1:k
    mG(i)=interp2(log10(lu_data.beta),lu_data.phi',squeeze(lu_data.mG(:,i,:)),log10(beta),phi,'spline');
end


end

