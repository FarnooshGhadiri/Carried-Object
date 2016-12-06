function [vote_record,valid_hypo_idx] = trace_back_vote_record(hypo_list, ...
    candidate_pos, vote_disc_rad, min_vote)
% [vote_record,valid_hypo_idx] = trace_back_vote_record(hypo_list,
% candidate_pos,scoresK_scale, vote_disc_rad)
% 
% 
%
%
% Liming Wang, Jan 2008
%vote_disc_rad=12;
vote_record = []; % zeros(size(match_pos,1),1);
%vote_disc_rad=6;
hypo_cnt = size(hypo_list,1);

valid_hypo  = zeros(hypo_cnt,1);
vote_disc_rad=7;
hypoIdx = 0;

% Elips type 2 ((y-h)^2)/a^2+((x-k)^2)/b^2=1  (h,k) is a elips center ...b
% is a major axis
Major_axis=14;
for hypo=1:hypo_cnt
    % Search in the circular neihborhood 
   % dist2= sqrt((candidate_pos(:,1)-hypo_list(hypo,1)).^2 + (candidate_pos(:,2)-hypo_list(hypo,2)).^2);
    %  [hypo_offset,match_id]  = sort(dist2,'ascend');
    %voter_id = find(hypo_offset<=vote_disc_rad);
 dist2=(((candidate_pos(:,1)-hypo_list(hypo,1)).^2 )/vote_disc_rad^2+ ((candidate_pos(:,2)-hypo_list(hypo,2)).^2)/Major_axis^2);
  [hypo_offset,match_id]  = sort(dist2,'ascend');
 voter_id = find(hypo_offset<=1);
    if(length(voter_id)>=min_vote)
        hypoIdx = hypoIdx + 1;
        vote_record(hypoIdx).voter_id   = match_id(voter_id);
        valid_hypo(hypo)=1;
    end
end

valid_hypo_idx  = find(valid_hypo);
