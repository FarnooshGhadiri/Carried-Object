function [new_voter_idx,new_score,uni_test_id,uni_code_id]=one2one_match(test_id,code_id,cost)
% [new_voter_idx,new_score,uni_test_id,uni_code_id]=one2one_match(test_id,code_id,cost)

code_id = double(code_id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   construct a matrix with #unique(test_id) x #unique(scoresK_id)    %
[uni_code_id,code_m,code_n] = unique(code_id);
[uni_test_id,test_m,test_n] = unique(test_id);
nb_voter    = length(test_m);
nb_code     = length(code_m);
A           = ones(nb_voter, nb_code)*2;
uni_id      = sub2ind([nb_voter,nb_code],test_n,code_n);
A(uni_id)   = cost; %1 - scoresK(v_id);
% one-2-one match constrains                          %
[cost,t_id,c_id]= myhungarian2(A, 2);
new_score   = length(t_id)  - cost;
%                recover vote_record.voter_id                         %
nb_voter_new= length(t_id);
% new_voter_id= zeros(nb_voter_new,1);
uni_test_id = uni_test_id(t_id);
uni_code_id = uni_code_id(c_id);

dim_t   = max(test_id);
dim_c   = max(code_id);

idx_mtx = sparse(test_id,double(code_id),[1:length(test_id)],dim_t,dim_c);

new_idx = sub2ind([dim_t,dim_c],uni_test_id,uni_code_id);

new_voter_idx	= full(idx_mtx(new_idx));
