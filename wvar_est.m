function [wvar,CI_low,CI_up,df]=wvar_est(W,j,jj,Lj,Ljj,Nj,Mjj,N,M,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wavelet variance and confidence interval estimation for specific level j,j'
%
%Input: W = wavelet-wavelet,scaling-wavelet or wavelet-scaling
%           coefficient matrix
%    j,j' = levels j, j'
% Lj, Ljj = length of level j-th, j'-th wavelet filters
% Nj, Mjj = number of points whose calculations is not affected by boundary
%    N, M = image size (N rows, M columns)
%       p = confidence level for (1-p) confidence interval calculation
%Output: wvar = wavelet variance for level j,j'
%      CI_low = lower confidence limit
%       CI_up = upper confidence limit
%          df = degrees of freedom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_cut=W(Lj:N,Ljj:M).^2;
wvar=mean(W_cut(:));

temp=zeros(2*Nj,2*Mjj);
temp(1:Nj,1:Mjj)=W(Lj:N,Ljj:M);
sW = ifft2(fft2(temp).*conj(fft2(temp)))/Nj/Mjj;
sigma_W=sum(sum(sW(:).^2));

if Nj*Mjj>=128   
    df=(2*Nj*Mjj*wvar^2)/sigma_W;
else
    df=max((Nj*Mjj)/(2^j*2^jj),1);
end

CI_low=df*wvar/chi2inv(1-(p/2),df);
CI_up=df*wvar/chi2inv(p/2,df);