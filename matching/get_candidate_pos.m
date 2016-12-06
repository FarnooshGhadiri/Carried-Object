function [candidate_pos] = get_candidate_pos(valid_idx, scoresK_id, relpos, testpos)
%[candidate_pos] = get_candidate_pos(scoresK, scoresK_id, relpos, testpos, vote_thresh)
% get candiates positions by matching position and relative distance to
% object center
%
%INPUT:
%   valid_idx:  valid vote indices to 'scoresK'
%   scoresK_id: indices to codebook entries
%   relpos:     the relative position computed from codebook entries
%   testpos:    feature points position
%OUTPUT:
%   candidate_pos: possible winner positions
%
%   Liming Wang, Jan 2008
%

nb_test = size(testpos,1);
testpos_idx = mod(valid_idx-1, nb_test) + 1;

candidate_pos   = zeros(length(valid_idx),2);
scoresK_id      = scoresK_id(valid_idx);
relpos          = relpos(scoresK_id,:);
testpos         = testpos(testpos_idx,:);
candidate_pos   = testpos+relpos;

