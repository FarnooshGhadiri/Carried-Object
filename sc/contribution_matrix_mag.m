function [contri_mtx] = contribution_matrix_mag(r,theta,mag,ori_edge,bin_r, nb_bin_theta, nb_ori,...
    blur_r,blur_t, blur_method)
%[contri_mtx] = contribution_matrix_mag(ori_edge,bin_r, nb_bin_theta,...
%     blur_o, blur_r,blur_t, blur_method)
%   this is an inner function and all the paremeters are supposed to be
%   passed in correctly.
% 
% Liming Wang, Jan, 2008
%

if (~exist('blur_r','var'))
    blur_r  = 0.2;
end

if (~exist('blur_t','var'))
    blur_t  = 0.8;
end

if (~exist('blur_method','var'))
    blur_method='gaussian';
end

if(strcmp(blur_method,'gaussian'))
    blur_method=1;
else % using cosine function to blur
    blur_method=2;
end

%
% for each of the orientation
% compute each (x,y)'s contribution to each bin
%   this contribution is stored in a matrix of
%       #(x,y)    x    #bin
%   each entry indicated the pixel (x,y)'s weight's in this bin
% 

nb_pix  = length(theta);
nb_r    = length(bin_r)-1;
nb_bins = nb_r*nb_bin_theta*nb_ori;

si = [];
sj = [];
val = [];

% prepare the blur function for each bin
rad_bins_len    = bin_r(2:end) - bin_r(1:end-1);
rad_bins_len2   = rad_bins_len/2;
rad_blur_part   = rad_bins_len*blur_r;
rad_bins        = bin_r;
rad_bins_mid    = (rad_bins(1:end-1) + rad_bins(2:end))/2;

rad_bins_lambda = rad_bins_len2+rad_blur_part;

ori_unit        = 2*pi/nb_bin_theta;
ori_unit2       = ori_unit/2;
ori_bins        = (0:nb_bin_theta)*ori_unit;
ori_bins_mid    = (ori_bins(1:end-1) + ori_bins(2:end))/2;
ori_blur_part   = ori_unit*blur_t;

ori_bins_lambda = ori_unit2+ori_blur_part;


if(blur_method==2)  % using consine blur
    rad_blur_part_scale   = (pi/2)./(rad_blur_part);
    ori_blur_part_scale   = (pi/2)/(ori_blur_part);
end
% % % % % % % % % % % % %  this is athree-level iteration % % % % % % % % %
% % % % % iterate on edge_ori                             % % % % % % % % %
% % % % % % % % % iterate on radial bin                   % % % % % % % % %
% % % % % % % % % % % % % iterate on angular bin          % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

nb_ori  = length(ori_edge);

% bin order: 1st dim is angular, 2nd order is radial, 3rd order is
% orientation, for each point
for ori_id  = 1:nb_ori
    rt_idx  = ori_edge(ori_id).idx;
    if(isempty(rt_idx))
        continue;
    end    
    mr  = r(rt_idx);
    mth = theta(rt_idx);
    
    ori_st_id=(ori_id-1)*nb_r*nb_bin_theta;    
    for r_bin=1:nb_r
        r_idx   = find(mr>(rad_bins(r_bin)-rad_blur_part(r_bin)) &...            
            mr< (rad_bins(r_bin+1)+rad_blur_part(r_bin)));
        if(isempty(r_idx))
            continue;
        end
        mr1     = mr(r_idx);
        mth1    = mth(r_idx);
        for t_bin=1:nb_bin_theta % To find head location we need to process this part
            if(t_bin==1)
                t_idx   = find((mth1>=ori_bins(t_bin)&mth1<ori_bins(t_bin+1)+ori_blur_part)...
                    |mth1>ori_bins(end)-ori_blur_part);
            elseif(t_bin==nb_bin_theta)
                t_idx   = find(mth1>(ori_bins(t_bin)-ori_blur_part) | mth1<ori_blur_part);
            else
                t_idx   = find(mth1>(ori_bins(t_bin)-ori_blur_part) &...
                    mth1<(ori_bins(t_bin+1)+ori_blur_part));
            end
            if(isempty(t_idx))
                continue;
            end
            mr2     = mr1(t_idx);
            mth2    = mth1(t_idx);
            
            r_dis   = abs(mr2 - rad_bins_mid(r_bin));
            t_dis   = abs(angle_diff(mth2,ori_bins_mid(t_bin)));
            
            if(blur_method==1)  % using gaussian blur function
                contri  = Gauss(r_dis,rad_bins_lambda(r_bin)).*Gauss(t_dis,ori_bins_lambda);
                si = [si; rt_idx(r_idx(t_idx))];
                sj = [sj; (ori_st_id+(r_bin-1)*nb_bin_theta+t_bin)*ones(length(t_idx), 1)];
                val = [val; contri];
            else
                r_dis_within_idx=find(r_dis<=rad_bins_len2(r_bin));
                t_dis_within_idx=find(t_dis<=ori_unit2);
                r_dis   = r_dis-rad_bins_len2(r_bin);
                t_dis   = t_dis-ori_unit2;
                r_dis(r_dis_within_idx) = 0;
                t_dis(t_dis_within_idx)  = 0;               
                contri  = cos(r_dis*rad_blur_part_scale(r_bin)).*cos(t_dis*ori_blur_part_scale);
                si = [si; rt_idx(r_idx(t_idx))];
                sj = [sj; (ori_st_id+(r_bin-1)*nb_bin_theta+t_bin)*ones(length(t_idx), 1)];
                val = [val; contri];
            end
        end
    end
end



% % % % normalize each point's weights to sum up to 1 % % % % %
contri_mtx = sparse(si, sj, val, nb_pix, nb_bins);


contri_1    = full(sum(contri_mtx,2));
contri_1(find(contri_1<eps)) = 1;

contri_mtx  = spmtimesd(contri_mtx,mag./contri_1,[]);


function y=Gauss(x,sigma)

% y=Gauss(x,sigma)
%
% calculate the Gaussian function of x
x=x(:);

y=exp(-x.^2/(2*sigma^2))/(sqrt(2*pi)*sigma);
