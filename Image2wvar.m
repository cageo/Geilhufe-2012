function Image2wvar(im,wfilter,alpha,mediancalc,plot_results,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%Mainfunction%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate wavelet variance estimates and their CIs, which will be saved in
%a separate *.mat-file with user specified filename. Example plots for the
%results can be displayed.
%
%Input:        im = grayscale image (matrix)
%         wfilter = type of wavelet filter ('haar', 'd4', 'la8', 'la16'),
%                   (default: 'haar')                    
%           alpha = significance level for (1-alpha) confidence intervals,
%                   (default: '0.05')
%      mediancalc = calculate also median estimates (true, false)
%                   (default: false)
%    plot_results = plot the results (true, false) (default: true)
%        filename = filename for saving the results (default: 'wvar_image')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set default values if not specified differently by user
if nargin == 1
    wfilter = 'haar';
    alpha=0.05;
    mediancalc = false;
    plot_results = true;
    filename = 'wvar_image';
elseif nargin == 2
    alpha=0.05;
    mediancalc = false;
    plot_results = true;
    filename = 'wvar_image';
elseif nargin == 3
    mediancalc = false;
    plot_results = true;
    filename = 'wvar_image';
elseif nargin == 4
    plot_results = true;
    filename = 'wvar_image';
elseif nargin == 5
    filename = 'wvar_image';
end

%assure that image is of precision double
im=double(im);

%load the desired basic wavelet filters for modwt
[h,g,L]=myfilter(wfilter); h=h./sqrt(2); g=g./sqrt(2);

%determine the maximum scale in each direction according to image size and
%filter type
switch wfilter
  case 'haar' 
      max_j_scale=floor(log2(size(im,1)));
      max_jj_scale=floor(log2(size(im,2)));
  case 'd4' 
      max_j_scale=floor(log2((size(im,1)-1)/(L-1)+1)); 
      max_jj_scale=floor(log2((size(im,2)-1)/(L-1)+1));
  case 'la8' 
      max_j_scale=floor(log2((size(im,1)-1)/(L-1)+1));
      max_jj_scale=floor(log2((size(im,2)-1)/(L-1)+1));
  case 'la16' 
      max_j_scale=floor(log2((size(im,1)-1)/(L-1)+1));
      max_jj_scale=floor(log2((size(im,2)-1)/(L-1)+1));
end


if mediancalc
    [wwvar,swvar,wsvar,df_ww,CI_low_ww,CI_up_ww,df_sw,CI_low_sw,CI_up_sw,df_ws,CI_low_ws,CI_up_ws,...
        wwvar_median,swvar_median,wsvar_median,df_ww_median,CI_low_ww_median,CI_up_ww_median,...
        df_sw_median,CI_low_sw_median,CI_up_sw_median,df_ws_median,CI_low_ws_median,CI_up_ws_median] = ...
        waveletVariance(im,max_j_scale,max_jj_scale,h,g,L,alpha,true);
else
    [wwvar,swvar,wsvar,df_ww,CI_low_ww,CI_up_ww,df_sw,CI_low_sw,CI_up_sw,df_ws,CI_low_ws,CI_up_ws,~]=...
        waveletVariance(im,max_j_scale,max_jj_scale,h,g,L,alpha,false);
end

%save all results
save(filename,'-regexp','^file|^wwvar|^swvar|^wsvar|^df|^CI|^im')

%plot results
if plot_results
    plot_wvar(filename,mediancalc);
end