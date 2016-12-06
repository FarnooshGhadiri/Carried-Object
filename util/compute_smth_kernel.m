function vote_filter=compute_smth_kernel(scale_height, asp_ratio)
%voter_filter=compute_smth_kernel()
% compute kernels
% 
% 
% Liming Wang, Jan 2008
% 


% these ratios are based on pedestrian height
% hypo_cluster_rad_ratio  = 0.1;
% ellipse_b_ratio         = 0.04;
% gauss_ellipse_b_ratio   = 0.06;

% ellipse_b_ratio         = 0.02;
% gauss_ellipse_b_ratio   = 0.03;
ellipse_b_ratio         = 0.03;
gauss_ellipse_b_ratio   = 0.05;



% Disc and smooth tool specification
% Ellipse disc (circular disc can be viewed as a special case)

disc_type = 0; % 0: Equally Weighted (all equals 1); 1: Gaussian Weighted
disc_sigma_ratio = 0.5; % Gaussian sigma (standard deviation) over a or b. Available only when disc is Gaussian Weighted.
% Smooth tool
% -- EllipseGaussian (Normal 2D Gaussian can be viewed as a special case)
smooth_sigma_ratio = 1/3; % Gaussian sigma (standard deviation) over a or b.

vote_filter = cell(length(scale_height),1);
for ns=1:length(scale_height)
    disc_b = round(scale_height(ns)*ellipse_b_ratio);
    disc_a = round(disc_b*asp_ratio(ns));
    smooth_b = round(scale_height(ns)*gauss_ellipse_b_ratio);
    smooth_a = round(smooth_b*asp_ratio(ns));
    % Generate disc
    if disc_type == 0 % Equally Weighted
        disc = double(genEllipse(2*disc_a+1,2*disc_b+1));
    else % Gaussian Weighted
        [x,y] = meshgrid(-disc_b:disc_b,-disc_a:disc_a);
        X = [x(:),y(:)];
        C = [disc_b*disc_sigma_ratio,0;0,disc_a*disc_sigma_ratio].^2;
        disc = gaussian(X, [0;0], C);
        disc = reshape(disc, 2*disc_a+1, 2*disc_b+1); % Do not need to normalize.
    end

    % Generate smoothing filter
    [x,y] = meshgrid(-smooth_b:smooth_b,-smooth_a:smooth_a);
    X = [x(:),y(:)];
    C = [smooth_b*smooth_sigma_ratio,0;0,smooth_a*smooth_sigma_ratio].^2;
    smooth_filter = gaussian(X, [0;0], C);
    smooth_filter = reshape(smooth_filter, 2*smooth_a+1, 2*smooth_b+1); % Do not need to normalize.
    % Generate final filter
    vote_filter{ns} = conv2(disc,smooth_filter);
end
