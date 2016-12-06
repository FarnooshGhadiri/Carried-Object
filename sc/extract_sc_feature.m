function [fea,fea_sum]  = extract_sc_feature(edge_map,theta_map,testpos,para_sc)
%[fea,fea_sum]  = extract_sc_feature(edge_map,theta_map,testpos,para_sc)
%   extract shape context feature from edge map for point set 'testpos'
%
%INPUT:
%   EDGE_MAP:   edge detection results
%   THETA_MAP:  orientation at each point on edge map
%   TESTPOS:    feature points
%   PARA_SC:    parameter for feature extraction
%
%OUTPUT:
%   FEA:        feature set.
%               [# of feature points x dim of the feature histogram]
%   FEA_SUM:    sum of each feature, used for normalization
%               [# of feature points x 1]
%
%   Liming Wang, Jan 2008
%

[ey,ex] = find(edge_map>para_sc.edge_thresh);
eind    = sub2ind(size(edge_map), ey, ex);

ori     = theta_map(eind);
mag     = edge_map(eind);

fea     = compute_sc(testpos(:,1),testpos(:,2),...
    ex,ey,ori,mag,...
    para_sc.bin_r,para_sc.nb_bin_theta,para_sc.nb_ori,...
    para_sc.blur_r,para_sc.blur_t,para_sc.blur_o,para_sc.blur_method);

fea_sum = sum(fea,2);
