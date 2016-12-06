function [cost,r_id,c_id] = myhungarian2(A, c_thresh)
%[cost,r_id,c_id] = myhungarian(A, c_thresh)
% 
% Liming Wang, Jan 2008

[m,n]   = size(A);

r_id    = mex_lap_munkres(A);

c_id    = find(r_id>0);
r_id    = r_id(c_id);

rc_ind  = sub2ind([m,n],r_id,c_id);

cost_a  = A(rc_ind);

valid_match_id  = find(cost_a<c_thresh);

r_id    = r_id(valid_match_id);
c_id    = c_id(valid_match_id);

rc_ind1 = sub2ind([m,n], r_id, c_id);
cost    = sum(A(rc_ind1));
