function [wvar_median,CI_low_median,CI_up_median,df_median]=...
    wvar_est_median(W,j,jj,Lj,Ljj,Nj,Mjj,N,M,C,K,nu_k,nu_kk,Vzero,...
    SumSqVzero,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Median estimator for wavelet variance and confidence interval for specific
%level j,j'
%
%Input:     W = wavelet-wavelet,scaling-wavelet or wavelet-scaling
%               coefficient matrix
%        j,j' = levels j, j'
%     Lj, Ljj = length of level j-th, j'-th wavelet filters
%     Nj, Mjj = number of points whose calculations is not affected by boundary
%        N, M = image size (N rows, M columns)
%           C = constant used for median and its CI estimation
%           K = number of Slepian tapers
% nu_k, nu_kk = \nu_{k,:}, \nu_{k',:}
%       Vzero = V_{k,k'}(0)
%  SumSqVzero = sum of squares of Vzero over all k,k' 
%           p = confidence level for (1-p) confidence interval calculation
%Output: wvar_median = wavelet variance for level j,j' (median estimate)
%      CI_low_median = lower confidence limit (median estimate)
%       CI_up_median = upper confidence limit (median estimate)
%          df_median = degrees of freedom (median estimate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=1:K; kk=1:K;
u=Lj:N; v=Ljj:M;

temp=log(W(u,v).^2);
muzero=median(temp(:));

Jzero=zeros(K,K);
Jzero(k,kk)=nu_k(k,u)*sign(log(W(u,v).^2)-muzero)*nu_kk(kk,v)';

mu=sum(sum(Jzero.*Vzero))/SumSqVzero;

Szero=1/(K^2)*sum(sum((Jzero-mu*Vzero).^2));

temp=W(u,v).^2;
wvar_median=median(temp(:))*exp(-Szero/(2*Nj*Mjj*C))/(norminv(0.75,0,1)^2);

if Nj*Mjj>=128   
    df_median=(2*Nj*Mjj*C)/Szero;
else
    df_median=max((Nj*Mjj)/(2^j*2^jj),1);
end

CI_low_median=df_median*wvar_median/chi2inv(1-(p/2),df_median);
CI_up_median=df_median*wvar_median/chi2inv(p/2,df_median);
