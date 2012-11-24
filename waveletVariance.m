function [wwvar,swvar,wsvar,df_ww,CI_low_ww,CI_up_ww,df_sw,CI_low_sw,...
    CI_up_sw,df_ws,CI_low_ws,CI_up_ws,wwvar_median,swvar_median,...
    wsvar_median,df_ww_median,CI_low_ww_median,CI_up_ww_median,...
    df_sw_median,CI_low_sw_median,CI_up_sw_median,df_ws_median,...
    CI_low_ws_median,CI_up_ws_median] = ...
    waveletVariance(im,max_j_scale,max_jj_scale,h,g,L,p,mediancalc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2-D pyramid algorithm for calculating wavelet variances
%
%Input:           im = image,
%        max_j_scale = maximum j scale according to image size
%       max_jj_scale = maximum j' scale according to image size
%                  h = MODWT level 1 wavelet filters
%                  g = MODWT level 1 scaling filters
%                  L = filter length
%                  p = significance level for (1-p) confidence intervals
%         mediancalc = calculate also median estimates (true, false)
%
%Output:       wwvar = wavelet-wavelet variance matrix (for levels j,j')
%              swvar = scaling-wavelet variance matrix (for levels j,j')
%              wsvar = wavelet-scaling variance matrix (for levels j,j')
%  df_ww,df_ws,df_sw = degrees of freedom
% CI_low_ww,CI_low_sw,CI_low_ws = lower confidence intervals
%    CI_up_ww,CI_up_sw,CI_up_ws = upper confidence intervals
%         ..._median = same as above, but for the median type estimator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M]=size(im);

W_temp=zeros(size(im));
Z_temp=zeros(size(im));
W=W_temp; U=W_temp; V=W_temp; Z=W_temp;
W_j=W; U_j=U; V_j=V; Z_j=Z;
W_jj_col=W; U_jj_col=U; V_jj_col=V; Z_jj_col=Z;
W_jj_row=W; U_jj_row=U; V_jj_row=V; Z_jj_row=Z;


%% mean case initialization
wwvar=zeros(max_j_scale,max_jj_scale); swvar=wwvar; wsvar=wwvar;
df_ww     = zeros(max_j_scale,max_jj_scale);
CI_low_ww = zeros(max_j_scale,max_jj_scale);
CI_up_ww  = zeros(max_j_scale,max_jj_scale);
df_sw     = zeros(max_j_scale,max_jj_scale);
CI_low_sw = zeros(max_j_scale,max_jj_scale);
CI_up_sw  = zeros(max_j_scale,max_jj_scale);
df_ws     = zeros(max_j_scale,max_jj_scale);
CI_low_ws = zeros(max_j_scale,max_jj_scale);
CI_up_ws  = zeros(max_j_scale,max_jj_scale);

%% median case initialization
C=4*(normpdf(norminv(0.75,0,1),0,1)*norminv(0.75,0,1))^2; %constant used for median and its CI estimation
K=5; %number of Slepian tapers used for median and its CI estimation
wwvar_median=zeros(max_j_scale,max_jj_scale); swvar_median=wwvar_median; wsvar_median=wwvar_median;
CI_low_ww_median = zeros(max_j_scale,max_jj_scale);
CI_up_ww_median  = zeros(max_j_scale,max_jj_scale);
CI_low_sw_median = zeros(max_j_scale,max_jj_scale);
CI_up_sw_median  = zeros(max_j_scale,max_jj_scale);
CI_low_ws_median = zeros(max_j_scale,max_jj_scale);
CI_up_ws_median  = zeros(max_j_scale,max_jj_scale);
df_ww_median     = zeros(max_j_scale,max_jj_scale);
df_sw_median     = zeros(max_j_scale,max_jj_scale);
df_ws_median     = zeros(max_j_scale,max_jj_scale);

%% 2-D pyramid algorithm

for j=1:max_j_scale
    
    if j==1
        input=im;
    else
        input=Z_new;
    end
    
    for col=1:M
        [W_temp(:,col),Z_temp(:,col)]=modwt(input(:,col),h,g,j);
    end
    
    for row=1:N
        input=W_temp;
        [W_j(row,:),V_j(row,:)]=modwt(input(row,:),h,g,j);
        input=Z_temp;
        [U_j(row,:),Z_j(row,:)]=modwt(input(row,:),h,g,j);
    end
    
    for jj=j:max(max_j_scale,max_jj_scale)
        
        Lj=(2^j-1)*(L-1)+1;
        Ljj=(2^jj-1)*(L-1)+1;
        Nj=N-Lj+1;
        Mjj=M-Ljj+1;
             
        if j==jj
            if jj<=max_jj_scale
                Z_new=Z_j;
                W=W_j;
                U=U_j;
                V=V_j;
                %% mean (diagonal expansion)
                [wwvar(j,jj),CI_low_ww(j,jj),CI_up_ww(j,jj),df_ww(j,jj)]=wvar_est(W,j,jj,Lj,Ljj,Nj,Mjj,N,M,p);
                [swvar(j,jj),CI_low_sw(j,jj),CI_up_sw(j,jj),df_sw(j,jj)]=wvar_est(U,j,jj,Lj,Ljj,Nj,Mjj,N,M,p);
                [wsvar(j,jj),CI_low_ws(j,jj),CI_up_ws(j,jj),df_ws(j,jj)]=wvar_est(V,j,jj,Lj,Ljj,Nj,Mjj,N,M,p);
                
                %% median (diagonal expansion)
                if mediancalc
                    [nu_k,nu_kk,Vzero,SumSqVzero] = CalculateNu(Lj,Ljj,Nj,Mjj,N,M,K);
                    [wwvar_median(j,jj),CI_low_ww_median(j,jj),CI_up_ww_median(j,jj),df_ww_median(j,jj)]=...
                        wvar_est_median(W,j,jj,Lj,Ljj,Nj,Mjj,N,M,C,K,nu_k,nu_kk,Vzero,SumSqVzero,p);
                    [swvar_median(j,jj),CI_low_sw_median(j,jj),CI_up_sw_median(j,jj),df_sw_median(j,jj)]=...
                        wvar_est_median(U,j,jj,Lj,Ljj,Nj,Mjj,N,M,C,K,nu_k,nu_kk,Vzero,SumSqVzero,p);
                    [wsvar_median(j,jj),CI_low_ws_median(j,jj),CI_up_ws_median(j,jj),df_ws_median(j,jj)]=...
                        wvar_est_median(V,j,jj,Lj,Ljj,Nj,Mjj,N,M,C,K,nu_k,nu_kk,Vzero,SumSqVzero,p);  
                end
                %%
            end
        else            
            if jj==j+1
                U_temp_col=U_j;
                Z_temp_col=Z_j;
                V_temp_row=V_j;
                Z_temp_row=Z_j;  
            else
                U_temp_col=U_jj_col;
                Z_temp_col=Z_jj_col;
                V_temp_row=V_jj_row;
                Z_temp_row=Z_jj_row;
            end

            if ((j<=max_jj_scale && max_j_scale>max_jj_scale) || (jj<=max_j_scale && max_j_scale<=max_jj_scale))
                for col=1:M
                    input=U_temp_col;
                    [W_jj_col(:,col),U_jj_col(:,col)]=modwt(input(:,col),h,g,jj);
                    input=Z_temp_col;
                    [V_jj_col(:,col),Z_jj_col(:,col)]=modwt(input(:,col),h,g,jj);  
                end
                Njj_col=N-Ljj+1;
                Mj_col=M-Lj+1;
                %% mean (vertical expansion)
                [wwvar(jj,j),CI_low_ww(jj,j),CI_up_ww(jj,j),df_ww(jj,j)]=wvar_est(W_jj_col,jj,j,Ljj,Lj,Njj_col,Mj_col,N,M,p);
                [swvar(jj,j),CI_low_sw(jj,j),CI_up_sw(jj,j),df_sw(jj,j)]=wvar_est(U_jj_col,jj,j,Ljj,Lj,Njj_col,Mj_col,N,M,p);
                [wsvar(jj,j),CI_low_ws(jj,j),CI_up_ws(jj,j),df_ws(jj,j)]=wvar_est(V_jj_col,jj,j,Ljj,Lj,Njj_col,Mj_col,N,M,p);   
                
                %% median (vertical expansion)
                if mediancalc
                    [nu_k_col,nu_kk_col,Vzero_col,SumSqVzero_col] = CalculateNu(Ljj,Lj,Njj_col,Mj_col,N,M,K);
                    [wwvar_median(jj,j),CI_low_ww_median(jj,j),CI_up_ww_median(jj,j),df_ww_median(jj,j)]=...
                        wvar_est_median(W_jj_col,jj,j,Ljj,Lj,Njj_col,Mj_col,N,M,C,K,nu_k_col,nu_kk_col,Vzero_col',SumSqVzero_col,p);
                    [swvar_median(jj,j),CI_low_sw_median(jj,j),CI_up_sw_median(jj,j),df_sw_median(jj,j)]=...
                        wvar_est_median(U_jj_col,jj,j,Ljj,Lj,Njj_col,Mj_col,N,M,C,K,nu_k_col,nu_kk_col,Vzero_col',SumSqVzero_col,p);
                    [wsvar_median(jj,j),CI_low_ws_median(jj,j),CI_up_ws_median(jj,j),df_ws_median(jj,j)]=...
                        wvar_est_median(V_jj_col,jj,j,Ljj,Lj,Njj_col,Mj_col,N,M,C,K,nu_k_col,nu_kk_col,Vzero_col',SumSqVzero_col,p);
                end
                %%
            end
            
            if (jj<=max_jj_scale)
                for row=1:N
                    input=V_temp_row;
                    [W_jj_row(row,:),V_jj_row(row,:)]=modwt(input(row,:),h,g,jj);
                    input=Z_temp_row;
                    [U_jj_row(row,:),Z_jj_row(row,:)]=modwt(input(row,:),h,g,jj);   
                end

                %% mean (horizontal expansion)
                [wwvar(j,jj),CI_low_ww(j,jj),CI_up_ww(j,jj),df_ww(j,jj)]=wvar_est(W_jj_row,j,jj,Lj,Ljj,Nj,Mjj,N,M,p);
                [swvar(j,jj),CI_low_sw(j,jj),CI_up_sw(j,jj),df_sw(j,jj)]=wvar_est(U_jj_row,j,jj,Lj,Ljj,Nj,Mjj,N,M,p);
                [wsvar(j,jj),CI_low_ws(j,jj),CI_up_ws(j,jj),df_ws(j,jj)]=wvar_est(V_jj_row,j,jj,Lj,Ljj,Nj,Mjj,N,M,p);
                
                %% median (horizontal expansion)
                if mediancalc
                    [nu_k,nu_kk,Vzero,SumSqVzero] = CalculateNu(Lj,Ljj,Nj,Mjj,N,M,K);
                    [wwvar_median(j,jj),CI_low_ww_median(j,jj),CI_up_ww_median(j,jj),df_ww_median(j,jj)]=...
                        wvar_est_median(W_jj_row,j,jj,Lj,Ljj,Nj,Mjj,N,M,C,K,nu_k,nu_kk,Vzero,SumSqVzero,p);
                    [swvar_median(j,jj),CI_low_sw_median(j,jj),CI_up_sw_median(j,jj),df_sw_median(j,jj)]=...
                        wvar_est_median(U_jj_row,j,jj,Lj,Ljj,Nj,Mjj,N,M,C,K,nu_k,nu_kk,Vzero,SumSqVzero,p);
                    [wsvar_median(j,jj),CI_low_ws_median(j,jj),CI_up_ws_median(j,jj),df_ws_median(j,jj)]=...
                        wvar_est_median(V_jj_row,j,jj,Lj,Ljj,Nj,Mjj,N,M,C,K,nu_k,nu_kk,Vzero,SumSqVzero,p);     
                end
            end
        end
    end
end
