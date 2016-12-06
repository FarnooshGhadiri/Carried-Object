function [recog_result] = sc_vote_1_scale_backup(Obj_Msk,edge_map,theta_map,codebook,para,verbose)
%[recog_result] = sc_vote_1_scale(edge_map,theta_map,codebook,para,verbose)
% voting on one scale image
% INPUT:
%   EDGE_MAP:   edge detection result of one scale
%   THETA_MAP:  orientation of points on edge map
%   CODEBOOK:   model information
%   PARA:       parameter
%   VERBOSE:    control debug output level
% Liming Wang, Jan 2008
global Center_Msk
if(~exist('verbose','var'))
    verbose=0;
end

if(verbose);
    last_time=cputime;
end

para_sc = para{1};
para_fea= para{2};
para_vote=para{3};
para_fea.K=20;
voter_filter    = para_vote.voter_filter;

[imgh,imgw]     = size(edge_map);
%testpos         = sample_grid_location([imgh,imgw], para_fea.sample_step);
sample_step=9;% Similar to the training set
testpos=sample_location_in_mask(Obj_Msk,sample_step);

if(para_sc.edge_bivalue)
    edge_map    = double(edge_map>para_sc.edge_thresh);
end

[fea,fea_sum]   = extract_sc_feature(edge_map, theta_map, testpos, para_sc);

fea_idx         = find(fea_sum>para_sc.sum_total_thresh);

if(length(fea_idx)<size(fea,1))
    if(verbose>1)
        disp(sprintf('prune %d zero test features',size(fea,1)-length(fea_idx)));
    end
    fea     = feature_from_ind(fea,fea_idx);
    testpos = feature_from_ind(testpos,fea_idx);
end

if(verbose>1)
    now_time    = cputime;
    fprintf(1,'Feature extraction: %f\n', now_time-last_time);
    last_time   = now_time;
end

[scoresK, scoresK_id]   = compute_matching_scores_bestK(fea,codebook.codes, codebook.sc_weight,...
    para_fea.K, para_sc.sum_total_thresh, para_fea.mask_fcn);

if(verbose>1)
    now_time    = cputime;
    fprintf(1,'Matching: %f\n', now_time-last_time);
    last_time   = now_time;
end

valid_vote_idx  = find(scoresK>para_vote.vote_thresh);
valid_scoresK_id=scoresK_id(valid_vote_idx);
Valid_img_id=codebook.obj_id(valid_scoresK_id);
% % % % % % % % %  NView=Number of samples in each class()% % % % %% % % %
% %
%
Person_view=mode(Valid_img_id);
%Nxt_view=Person_view+1;
if Person_view==3 | Person_view==4
    Person_view=[3 4];
else
    if Person_view==7 | Person_view==8
        Person_view=[7 8];
     end

end

valid_vote_idx=find(ismember(Valid_img_id,Person_view));
    % valid_vote_idx=find(Valid_img_id==Person_view | Valid_img_id==Nxt_view | Valid_img_id==Pre_view);
    
    valid_scoresK_id=scoresK_id(valid_vote_idx); % Recalculate those score that are related to just one view
    
    [candidate_pos] = get_candidate_pos(valid_vote_idx, scoresK_id, codebook.relpos, testpos);
    
    [hypo_list, score_list, vote_map, winThreshold]	= get_hypo_center(candidate_pos, scoresK, ...
        [imgh,imgw], valid_vote_idx, para_vote.vote_offset, voter_filter);
    
    % % % % % cluter hypo centers, avoid too close hypo because of scale or
    % % % % % deformation
    [hypo_list, score_list] = clusterHypo(hypo_list, score_list, [], ...
        para_vote.elps_ab, para_vote.nb_iter);
    
    % % % %% %% % % % % % % % % Find the hypo which is close to the center box % %  %% %
    % % %
    Disx=hypo_list(:,1)-Center_Msk(:,1);
    Disy=hypo_list(:,2)-Center_Msk(:,2);
    Dis=sqrt(Disx.^2+Disy.^2);
    [m1,ind]=min(Dis);
    %idx=find(Dis(ind)<40 & abs(Disx(ind))<25 & abs(Disx(ind))<25);
    hypo_list=hypo_list(ind,:);
    score_list=score_list(ind,:);
    
    % % % % % % trace back to find voters for each hypo % % % % % % % % % %
    [vote_record,valid_hypo_idx] = trace_back_vote_record(hypo_list, ...
        candidate_pos, para_vote.vote_disc_rad, para_vote.min_vote);


while m1>40 | isempty(vote_record)
    Second_fqt=find(~ismember(Valid_img_id,Person_view));
    Person_view2=mode(Valid_img_id(Second_fqt));
        
        if Person_view2==7 | Person_view2==8
            Person_view2=[7 8];
        end
        
    Person_view=[Person_view Person_view2];
    valid_vote_idx=find(ismember(Valid_img_id,Person_view2));
    % valid_vote_idx=find(Valid_img_id==Person_view | Valid_img_id==Nxt_view | Valid_img_id==Pre_view);
    
    valid_scoresK_id=scoresK_id(valid_vote_idx); % Recalculate those score that are related to just one view
    
    [candidate_pos] = get_candidate_pos(valid_vote_idx, scoresK_id, codebook.relpos, testpos);
    
    [hypo_list, score_list, vote_map, winThreshold]	= get_hypo_center(candidate_pos, scoresK, ...
        [imgh,imgw], valid_vote_idx, para_vote.vote_offset, voter_filter);
    
    % % % % % cluter hypo centers, avoid too close hypo because of scale or
    % % % % % deformation
    [hypo_list, score_list] = clusterHypo(hypo_list, score_list, [], ...
        para_vote.elps_ab, para_vote.nb_iter);
    
    % % % %% %% % % % % % % % % Find the hypo which is close to the center box % %  %% %
    % % %
    Disx=hypo_list(:,1)-Center_Msk(:,1);
    Disy=hypo_list(:,2)-Center_Msk(:,2);
    Dis=sqrt(Disx.^2+Disy.^2);
    [m1,ind]=min(Dis);
    %idx=find(Dis(ind)<40 & abs(Disx(ind))<25 & abs(Disx(ind))<25);
    hypo_list=hypo_list(ind,:);
    score_list=score_list(ind,:);
    
    % % % % % % trace back to find voters for each hypo % % % % % % % % % %
    [vote_record,valid_hypo_idx] = trace_back_vote_record(hypo_list, ...
        candidate_pos, para_vote.vote_disc_rad, para_vote.min_vote);
end
%  hypo_list=hypo_list(idx(1,1),:);
%  score_list=score_list(idx(1,1),:);





if(verbose>1)
    now_time    = cputime;
    fprintf(1,'Voting: %f\n', now_time-last_time);
    last_time   = now_time;
end
%%%%% save results for post processing %%%%%%%%%%%
recog_result.edge       = edge_map;
recog_result.theta      = theta_map;
recog_result.imgsz      = [imgh,imgw];
recog_result.testpos    = testpos;
recog_result.features   = fea;
recog_result.scoresK    = scoresK;
recog_result.scoresK_id   = scoresK_id;
recog_result.candidate_pos= candidate_pos;
recog_result.vote_record  = vote_record;
recog_result.vote_map     = vote_map;
recog_result.hypo_list    = hypo_list;
recog_result.score_list   = score_list;
recog_result.valid_vote_idx = valid_vote_idx;
recog_result.win_thresh = winThreshold;
%%%%% save end %%%%%%%%%