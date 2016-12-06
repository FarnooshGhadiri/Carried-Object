function [recog_result] = sc_vote_1_scale(Obj_Msk,edge_map,theta_map,codebook,para,verbose)
%[recog_result] = sc_vote_1_scale(edge_map,theta_map,codebook,para,verbose)
% voting on one scale image
% INPUT:
%   Obj_Msk:   Foreground Mask
%   EDGE_MAP:   edge detection result of one scale
%   THETA_MAP:  orientation of points on edge map
%   CODEBOOK:   model information
%   PARA:       parameter
%   VERBOSE:    control debug output level

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
para_fea.K=80;
voter_filter    = para_vote.voter_filter;

[imgh,imgw]     = size(edge_map);
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
%para_vote.vote_thresh=0.7; % change the thresh to see the effect
valid_vote_idx  = find(scoresK>para_vote.vote_thresh);
valid_scoresK_id=scoresK_id(valid_vote_idx);
Valid_img_id=codebook.view(valid_scoresK_id);

 [candidate_pos] = get_candidate_pos(valid_vote_idx, scoresK_id, codebook.relpos, testpos);
    
    [hypo_list, score_list, vote_map, winThreshold]	= get_hypo_center(candidate_pos, scoresK, ...
        [imgh,imgw], valid_vote_idx, para_vote.vote_offset, voter_filter);
    
    % % % % % cluter hypo centers, avoid too close hypo because of scale or
    % % % % % deformation
    [hypo_list, score_list] = clusterHypo(hypo_list, score_list, [], ...
        para_vote.elps_ab, para_vote.nb_iter);

    % % % % % % trace back to find voters for each hypo % % % % % % % % % %
    [vote_record_T,valid_hypo_idx] = trace_back_vote_record(hypo_list, ...
        candidate_pos, para_vote.vote_disc_rad, para_vote.min_vote);
    
%      % Choose the center that is voted more by top part of person's body
     Dst_up=[];
     for i=1:length(valid_hypo_idx)
          N_row=max(testpos(:,1));
          nb_test = size(testpos,1);
          testpos_idx=mod(vote_record_T(1,i).voter_id, nb_test) + 1;
          Valid_testpos=testpos(testpos_idx,2);
          Up_part_Dst=length(find(Valid_testpos<(N_row/(3))));% distribution of vote at upper body part
          Dst_up=[Dst_up Up_part_Dst];
     end
     [Mvalue,idx]=max(Dst_up);
     Dst_up(idx)=0;
     [Mvalue2,idx2]=max(Dst_up);
     if abs(Mvalue-Mvalue2)<5
         [MaxS]=max(score_list(idx,1),score_list(idx2,1));
         idx=find(score_list(:,1)==MaxS);
                
    end
% % % % % % % % %  NView=Number of samples in each class()% % % % %% % % %
Dis_Cnter=sqrt((hypo_list(:,1)-Center_Msk(1,1)).^2+(hypo_list(:,2)-Center_Msk(1,2)).^2);

if Dis_Cnter(1)>=35 
    if length(Dis_Cnter)>1 & min(Dis_Cnter(:,1))<35
    [min_value,idx]=min(Dis_Cnter(:,1));
    else
        hypo_list(1,1)=Center_Msk(1,1);
        hypo_list(1,2)=Center_Msk(1,2);
        [vote_record_T,valid_hypo_idx] = trace_back_vote_record(hypo_list, ...
        candidate_pos, para_vote.vote_disc_rad, para_vote.min_vote);
    idx=1;
    end
else
    idx=1; %Choose the center that gets more vore
  end
    




Person_view=mode(Valid_img_id(vote_record_T(1,idx).voter_id));
 
 
% add hypo that is close to the high probable hypo

Dis_hypo=sqrt((hypo_list(idx,1)-hypo_list(:,1)).^2+((hypo_list(idx,2)-hypo_list(:,2)).^2));
Valid_Center=find(Dis_hypo<30);
vote_record.voter_id=[];
for i=1:size(Valid_Center,1)
    % % % % Valid vote on the specific view %  % % % 

    Valid_Vote_view=find(Valid_img_id(vote_record_T(1,Valid_Center(i,1)).voter_id)==Person_view);
    vote_record_T(1,Valid_Center(i,1)).voter_id=vote_record_T(1,Valid_Center(i,1)).voter_id(Valid_Vote_view);
    vote_record.voter_id=[vote_record.voter_id;vote_record_T(1,Valid_Center(i,1)).voter_id];
end
 

hypo_list=hypo_list(idx,:);
score_list=score_list(idx,:);

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
recog_result.PsnVeiw=Person_view;
%%%%% save end %%%%%%%%%