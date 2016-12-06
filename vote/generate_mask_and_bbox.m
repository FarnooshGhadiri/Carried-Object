function [recog_result_new] =generate_mask_and_bbox(recog_result, codebook,...
    max_score_thresh,para_vote, verbose)
%[recog_result_new] =generate_mask_and_bbox(recog_result, codebook,
%max_score_thresh,para_vote, verbose) 
% 
para_vote.min_height=90;
para_vote.max_height=300;
if(~exist('verbose','var') || isempty(verbose))
    verbose = 0;
end 

recog_result_new    = [];

    recog_res   = generate_mask_and_bbox_1_scale(recog_result, codebook, para_vote, max_score_thresh);
    
    % prune hypo by mask heights
    mask_heights    = recog_res.mask_heights;
       if length(recog_res.vote_record)==2
        recog_res.vote_record(1,1).voter_id=[recog_res.vote_record(1,1).voter_id ; recog_res.vote_record(1,2).voter_id];
    end
   

    
    recog_result_new= [recog_result_new,recog_res];
    if(verbose>3)
        display_hypo_rect_mask_wrapper( recog_res );
    end    

