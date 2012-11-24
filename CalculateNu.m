function [nu_k,nu_kk,Vzero,SumSqVzero] = CalculateNu(Lj,Ljj,Nj,Mjj,N,M,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculations needed for median wavelet variance estimator and its CIs
%
%Input:      Lj, Ljj = length of level j-th, j'-th wavelet filters
%            Nj, Mjj = number of points whose calculations is not affected
%                      by the boundary
%               N, M = image size (N rows, M columns)
%                  K = number of Slepian tapers
%Output: nu_k, nu_kk = \nu_{k,:}, \nu_{k',:}
%              Vzero = V_{k,k'}(0)
%         SumSqVzero = sum of squares of Vzero over all k,k' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nu_k=zeros(K,N);
nu_kk=zeros(K,M);

Vzero=zeros(K,K);

BW=1/2*7;

if Nj>2*BW
    nu_k(:,Lj:N)=dpss(Nj,BW,K)';
else
    nu_k=NaN(K,N);
end

if Mjj>2*BW
    nu_kk(:,Ljj:M)=dpss(Mjj,BW,K)';
else
    nu_kk=NaN(K,M);
end

u=Lj:N; v=Ljj:M;
k=1:K; kk=1:K;
Vzero(k,kk)=nu_k(k,u)*ones(length(u),length(v))*nu_kk(kk,v)';
SumSqVzero=sum(sum(Vzero.^2));