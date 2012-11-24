function [beta_wls,y] = wvar_slope(J1,J0,wvar,df)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate slope for diagonal wavelet variances by linear weighted least
%squares (wls) according to section 9.5 in Percival and Walden (2000)
%
%Input: J1   = scale from where slope should start
%       wvar = wavelet variance
%              (i.e. wavelet-wavelet,wavelet-scaling or scaling-wavelet
%                    variance)
%       df   = degrees of freedom for wvar
%
%Output: beta_wls = wls estimator beta
%             y   = linear wls slope from scale J1 to J0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=ones(J0-J1+1,2);
Sigma_e=zeros(J0-J1+1);
Y=zeros(J0-J1+1,1);
for i=1:J0-J1+1
    j=i+J1-1;
    A(i,2)=log2(2^(j-1));
    Sigma_e(i,i)=psi(1,df(j,j)/2);
    Y(i)=log2(wvar(j,j))-psi(0,df(j,j)/2)+log(df(j,j)/2);
end
beta_wls=lscov(A,Y,Sigma_e);
y=zeros(size(Y));
for i=1:J0-J1+1
    j=i+J1-1;
    y(i)=beta_wls(1)+beta_wls(2)*log2(2^(j-1));
end