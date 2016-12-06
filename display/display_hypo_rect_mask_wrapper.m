function display_hypo_rect_mask_wrapper(recog_result)
%UNTITLED1 Summary of this function goes here
%   Detailed explanation goes here

edge_map     = recog_result.edge;
vote_map    = recog_result.vote_map;
hypo_list   = recog_result.hypo_list;
score_list  = recog_result.score_list;
hypo_bbox   = recog_result.hypo_bbox;
vote_record = recog_result.vote_record;
valid_vote_idx= recog_result.valid_vote_idx;
testpos     = recog_result.testpos;
hypo_mask   = recog_result.hypo_mask;
display_hypo_rect_mask(edge_map, vote_map,...
    hypo_list,score_list,hypo_bbox,...
    vote_record, valid_vote_idx, testpos,hypo_mask);
