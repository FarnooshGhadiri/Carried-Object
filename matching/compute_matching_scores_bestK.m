function [scoresK,scoresK_id] = compute_matching_scores_bestK(fea,codes,...
    codes_weight, K, sc_threshold, mask_fcn)
%[scoresK,scoresK_id]%=compute_matching_scores_bestK(fea,codes,para)
%compute matching score between test features and codebook features
%
% 
%INPUT:
%   fea:        test features, [#feature points x #dim of feature]
%   codes:      codes to be matched, [#code  x #dim of feature]
%   codes_weight:used to weight the distance function
%   K:          compute best K matches
%   sc_threshold: to prune those zero features
%   mask_fcn:   whether to use mask function or not, 1:use, 0: not use
%OUTPUT:
%   scoresK:    scores for best K matches for each feature points,
%                   [#feature points  x K]
%   scoresK_id: indicate which code is matched, [#feature pts  x  K]
%
%   Liming Wang, Jan 2008
%

if(~exist('mask_fcn','var'))
    mask_fcn=1;
end

if(mask_fcn)
    [scoresK,scoresK_id] = mex_compute_mscores_K_chi2(fea,codes,codes_weight,K,sc_threshold);
else
    fea_sum = sum(fea,2);
    fea	= spmtimesd(sparse(fea),1./fea_sum,[]);
    fea = full(fea);
    [scoresK,scoresK_id] = hist_cost_chi2(fea,codes,K);
end

