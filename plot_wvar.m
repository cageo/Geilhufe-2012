function plot_wvar(filename,mediancalc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Different plots of wavelet variances:
% - boxplot of diagonal scales for ww,sw and ws variance for mean and
%   median including confidence intervals
% - boxplot of all scale combinations for sw (column-wise) and ws
%   (row-wise) variance for mean and median including confidence intervals
% - plots of sw and ws variance matrix (i.e. all scale combinations)
%   for mean and median
%
%Input:  filename = name of the .mat-file with the wavelet variance results
%      mediancalc = plot median cases (true / false)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize mean estimations
load([filename,'.mat'])

% calculate min and max axis values
if mediancalc
    min_CI_sw =min(min(swvar(:)),min(min(swvar_median(:)),min(min(CI_low_sw(:)),min(CI_low_sw_median(:)))));
    max_CI_sw =max(max(swvar(:)),max(max(swvar_median(:)),max(max(CI_up_sw(:)),max(CI_up_sw_median(:)))));
    min_CI_ws =min(min(wsvar(:)),min(min(wsvar_median(:)),min(min(CI_low_ws(:)),min(CI_low_ws_median(:)))));
    max_CI_ws =max(max(wsvar(:)),max(max(wsvar_median(:)),max(max(CI_up_ws(:)),max(CI_up_ws_median(:)))));
    min_CI_ww =min(min(wwvar(:)),min(min(wwvar_median(:)),min(min(CI_low_ww(:)),min(CI_low_ww_median(:)))));
    max_CI_ww =max(max(wwvar(:)),max(max(wwvar_median(:)),max(max(CI_up_ww(:)),max(CI_up_ww_median(:)))));
else
    min_CI_sw =min(min(swvar(:)),min(CI_low_sw(:)));
    max_CI_sw =max(max(swvar(:)),max(CI_up_sw(:)));
    min_CI_ws =min(min(wsvar(:)),min(CI_low_ws(:)));
    max_CI_ws =max(max(wsvar(:)),max(CI_up_ws(:)));
    min_CI_ww =min(min(wwvar(:)),min(CI_low_ww(:)));
    max_CI_ww =max(max(wwvar(:)),max(CI_up_ww(:)));
end
min_value_sw_ws=log2(min(min_CI_sw,min_CI_ws));
max_value_sw_ws=log2(max(max_CI_sw,max_CI_ws));  
min_CI_ww = log2(min_CI_ww);
max_CI_ww = log2(max_CI_ww);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% boxplots for diagonal scales
%wavelet-wavelet variance
wwvar_mean_diag=log2(diag(wwvar));
low_ww_mean_diag=wwvar_mean_diag-log2(diag(CI_low_ww));
up_ww_mean_diag=log2(diag(CI_up_ww))-wwvar_mean_diag;

%calculate wls slope from scale level J1 to J0 for ww
if(size(df_ww,1)>3)
    J1 = 3;
else
    J1 = 1;
end
J0 = min(size(df_ww));

[beta_wls,y] = wvar_slope(J1,J0,wwvar,df_ww);

%boxplot diag wwvar including slope
fullscreen = get(0,'ScreenSize');
figure('Position',[0 0 3/4*fullscreen(3) 3/4*fullscreen(4)])
subplot(1,3,1), hold on, errorbar(1:size(wwvar_mean_diag,1),wwvar_mean_diag,low_ww_mean_diag,up_ww_mean_diag,'o');
                if mediancalc    
                    wwvar_median_diag=log2(diag(wwvar_median));
                    low_ww_median_diag=wwvar_median_diag-log2(diag(CI_low_ww_median));
                    up_ww_median_diag=log2(diag(CI_up_ww_median))-wwvar_median_diag;
                    errorbar(1:size(wwvar_median_diag,1),wwvar_median_diag,low_ww_median_diag,up_ww_median_diag,'xr');
                end
                plot(J1:J0,y,'k-','LineWidth',1);
                axis([0 size(wwvar_mean_diag,1)+1 (1-sign(min_CI_ww)*0.02)*min_CI_ww (1+sign(max_CI_ww)*0.02)*max_CI_ww]);
                xlabel('Level j=j'''); ylabel(['WWVAR on log2 scale, \beta = ',num2str(beta_wls(2))]); 
                axis square

%scaling-wavelet variance
swvar_mean_diag=log2(diag(swvar));  
low_sw_mean_diag=swvar_mean_diag-log2(diag(CI_low_sw));
up_sw_mean_diag=log2(diag(CI_up_sw))-swvar_mean_diag;

%calculate slope from scale level J1 to J0 for sw
if(size(df_sw,1)>3)
    J1 = 3;
else
    J1 = 1;
end
J0 = min(size(df_sw));
[beta_wls,y] = wvar_slope(J1,J0,swvar,df_sw);

%boxplot diag swvar including slope
subplot(1,3,2), hold on, errorbar(1:size(swvar_mean_diag,1),swvar_mean_diag,low_sw_mean_diag,up_sw_mean_diag,'o');
                if mediancalc
                    swvar_median_diag=log2(diag(swvar_median));  
                    low_sw_median_diag=swvar_median_diag-log2(diag(CI_low_sw_median));
                    up_sw_median_diag=log2(diag(CI_up_sw_median))-swvar_median_diag;
                    errorbar(1:size(swvar_median_diag,1),swvar_median_diag,low_sw_median_diag,up_sw_median_diag,'xr');
                end
                plot(J1:J0,y,'k-','LineWidth',1);
                axis([0 size(swvar_mean_diag,1)+1 (1-sign(min_value_sw_ws)*0.02)*min_value_sw_ws (1+sign(max_value_sw_ws)*0.02)*max_value_sw_ws]);
                xlabel('Level j=j'''); ylabel(['SWVAR on log2 scale, \beta = ',num2str(beta_wls(2))]); 
                axis square

%wavelet-scaling variance                
wsvar_mean_diag=log2(diag(wsvar));  
low_ws_mean_diag=wsvar_mean_diag-log2(diag(CI_low_ws));
up_ws_mean_diag=log2(diag(CI_up_ws))-wsvar_mean_diag;

%calculate slope from scale level J1 to J0 for ws
if(size(df_ws,1)>3)
    J1 = 3;
else
    J1 = 1;
end
J0 = min(size(df_ws));
[beta_wls,y] = wvar_slope(J1,J0,wsvar,df_ws);

%boxplot diag wsvar including slope
subplot(1,3,3), hold on, errorbar(1:size(wsvar_mean_diag,1),wsvar_mean_diag,low_ws_mean_diag,up_ws_mean_diag,'o');
                if mediancalc
                    wsvar_median_diag=log2(diag(wsvar_median));  
                    low_ws_median_diag=wsvar_median_diag-log2(diag(CI_low_ws_median));
                    up_ws_median_diag=log2(diag(CI_up_ws_median))-wsvar_median_diag;
                    errorbar(1:size(wsvar_median_diag,1),wsvar_median_diag,low_ws_median_diag,up_ws_median_diag,'xr');
                end
                plot(J1:J0,y,'k-','LineWidth',1);
                axis([0 size(wsvar_mean_diag,1)+1 (1-sign(min_value_sw_ws)*0.02)*min_value_sw_ws (1+sign(max_value_sw_ws)*0.02)*max_value_sw_ws]);
                xlabel('Level j=j'''); ylabel(['WSVAR on log2 scale, \beta = ',num2str(beta_wls(2))]); 
                axis square                               
hold off
set(gcf,'NextPlot','add');
axes;
h = title(['Wavelet variances for ',filename],'FontSize',14);
set(gca,'Visible','off');
set(h,'Visible','on'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% boxplots for all scale-combinations for swvar, wsvar (column- and rowwise)
%scaling-wavelet mean
swvar_mean=log2(swvar(:));
low_sw_mean=swvar_mean-log2(CI_low_sw(:));
up_sw_mean=log2(CI_up_sw(:))-swvar_mean;

%wavelet-scaling (transposed) mean
wsvar_mean_transposed=log2(wsvar)'; wsvar_mean_transposed=wsvar_mean_transposed(:);
temp=CI_low_ws'; temp2=CI_up_ws';
low_ws_mean_transposed=wsvar_mean_transposed-log2(temp(:));
up_ws_mean_transposed=log2(temp2(:))-wsvar_mean_transposed;

if mediancalc
    %scaling_wavelet median  
    swvar_median_=log2(swvar_median(:));
    low_sw_median=swvar_median_-log2(CI_low_sw_median(:));
    up_sw_median=log2(CI_up_sw_median(:))-swvar_median_;
    
    %wavelet-scaling (transposed) median
    wsvar_median_transposed=log2(wsvar_median)'; wsvar_median_transposed=wsvar_median_transposed(:);
    temp=CI_low_ws_median'; temp2=CI_up_ws_median';
    low_ws_median_transposed=wsvar_median_transposed-log2(temp(:));
    up_ws_median_transposed=log2(temp2(:))-wsvar_median_transposed;
end

%boxplot swvar column-wise
fullscreen = get(0,'ScreenSize');
figure('Position',[0 0 3/4*fullscreen(3) 3/4*fullscreen(4)])
subplot(1,2,1), hold on, errorbar(1:size(swvar_mean,1),swvar_mean,low_sw_mean,up_sw_mean,'o');
                if mediancalc
                    errorbar(1:size(swvar_median_,1),swvar_median_,low_sw_median,up_sw_median,'xr'),
                end
                xlim([0 size(swvar_mean,1)+1]); xlabel('Level combination j,j'' column-wise')
                ylim([(1-sign(min_value_sw_ws)*0.02)*min_value_sw_ws (1+sign(max_value_sw_ws)*0.02)*max_value_sw_ws]);
                ylabel('SWVAR on log2 scale');
                axis square
  
                %marking of diagonal values
                temp=sqrt(size(swvar_mean,1));
                for i=0:temp-1
                    j=i*(temp+1)+1;
                    plot([j j],[min_value_sw_ws max_value_sw_ws],'k:','LineWidth',3)
                end
                
%boxplot wsvar row-wise
subplot(1,2,2), hold on, errorbar(1:size(wsvar_mean_transposed,1),wsvar_mean_transposed,low_ws_mean_transposed,up_ws_mean_transposed,'o');
                if mediancalc
                    errorbar(1:size(wsvar_median_transposed,1),wsvar_median_transposed,low_ws_median_transposed,up_ws_median_transposed,'xr');
                end
                xlim([0 size(wsvar_mean_transposed,1)+1]); xlabel('Levels combination j,j'' row-wise')
                ylim([(1-sign(min_value_sw_ws)*0.02)*min_value_sw_ws (1+sign(max_value_sw_ws)*0.02)*max_value_sw_ws]);
                ylabel('WSVAR on log2 scale');
                axis square
                  
                %marking of diagonal values
                temp=sqrt(size(wsvar_mean_transposed,1));
                for i=0:temp-1
                    j=i*(temp+1)+1;
                    plot([j j],[min_value_sw_ws max_value_sw_ws],'k:','LineWidth',3)
                end 
hold off
set(gcf,'NextPlot','add');
axes;
h = title(['All combinations of scaling-wavelet and transposed wavelet-scaling variances for ',filename],'FontSize',14);
set(gca,'Visible','off');
set(h,'Visible','on'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% sw/ws variance plots mean_median

fullscreen = get(0,'ScreenSize');
figure('Position',[0 0 3/4*fullscreen(3) 3/4*fullscreen(4)])
if mediancalc
    subplot(2,2,1), imagesc(log2(swvar)), colorbar, caxis([min_value_sw_ws,max_value_sw_ws]), xlabel('Level j'), ylabel('Level j'''), set(gca,'YDir','normal'), axis square, title(['SWVAR (mean) on log2 scale for ',filename]);
    subplot(2,2,3), imagesc(log2(wsvar')), colorbar, caxis([min_value_sw_ws,max_value_sw_ws]), xlabel('Level j'''), ylabel('Level j'), set(gca,'YDir','normal'), axis square, title(['WSVAR (mean) transposed on log2 scale for ',filename]);
    subplot(2,2,2), imagesc(log2(swvar_median)), colorbar, caxis([min_value_sw_ws,max_value_sw_ws]), xlabel('Level j'), ylabel('Level j'''), set(gca,'YDir','normal'), axis square, title(['SWVAR (median) on log2 scale for ',filename]);
    subplot(2,2,4), imagesc(log2(wsvar_median')), colorbar, caxis([min_value_sw_ws,max_value_sw_ws]), xlabel('Level j'''), ylabel('Level j'), set(gca,'YDir','normal'), axis square, title(['WSVAR (median) transposed on log2 scale for ',filename]);
else
    subplot(1,2,1), imagesc(log2(swvar)), colorbar, caxis([min_value_sw_ws,max_value_sw_ws]), xlabel('Level j'), ylabel('Level j'''), set(gca,'YDir','normal'), axis square, title(['SWVAR (mean) on log2 scale for ',filename]);
    subplot(1,2,2), imagesc(log2(wsvar')), colorbar, caxis([min_value_sw_ws,max_value_sw_ws]), xlabel('Level j'''), ylabel('Level j'), set(gca,'YDir','normal'), axis square, title(['WSVAR (mean) transposed on log2 scale for ',filename]);
end