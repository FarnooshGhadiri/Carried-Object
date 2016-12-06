function para=set_parameter(codebook,ratio)
%para=set_parameter(codebook, ratio)
% set parameters for detection
%INPUT: 
%   CODEBOOK:   this variable could be empty or not. In case of empty, you
%   are setting parameters for training phrase. Then only parameters
%   related to shape context will be returned. If not empty, then
%   parameters for shape context will derived from codebook.
%OUTPUT:
%   PARA:   Parameters for shape context, feature extraction and matching,
%   voting will be returned if codebook is not empty. Only parameters 
%   related to shape context will be returned if codebook is empty.
% 

% 

%%%%%%%%%%%%%%%%%%%%%%%(I) parameter for training%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following parameters are for training stage, these parameters should
% be set up during the training phase, which is codebook generating. these 
% parameters should be fixed once the codebook is generated. You could 
% possibly change them if you really want to, but I suggest not.
% 

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!! I will put '+' ahead of the parameter comment, indicating
% you need to pay attention to this paramter; '-' means you can ignore them
% most of the time. But still, paying attention does not mean you have to
% change it, you just need to check if your object is different from
% pedestrian people; or your object has different height !!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(~exist('codebook','var') || isempty(codebook))
    
    %++++++++++ important parameter for codebook
    para_sc.model_height    = 170;
    para_sc.model_width=80;
    
    % 1. parameters for edge detection
    
    %+++++++++ which edge detect should be used, take value of ['pb'|'flt']
    para_sc.detector    = 'MR'; % DR
    %+++++++++ whether extract shape context on edge magnitude or not, take
    % value among [0|1], 0: use magnitude, 1: use edge as binary value
    para_sc.edge_bivalue= 0;
    
    %--------- the threshold for edge map, edge magnitude above this
    % threshold will be used to compute shape context
    if(strcmp(para_sc.detector,'MR'))
        para_sc.edge_thresh = 0.17;
    else
        %'flt' will produce more clutter, so set it higher to prune some of
        %the clutters
        para_sc.edge_thresh = 0.37;         
    end
    
    
    % 2. parameters for shape context feature
    %++++++++++ radial bin for shape context
    %  here is an illustration of bin_r
    % |<--6-->|
    % |<--------15------->|
    % so the actual bin value would be [ 6, 9], there are 2 bins.
    para_sc.bin_r   = [0,6,15];
    %++++++++++ number of angular bins
    para_sc.nb_bin_theta    = 12;
    %+++++++ number of orientation to split edge
    para_sc.nb_ori  = 4;
    %+++++++ blurring factors, blur_r is on radial, blur_t is angular,
    %blur_o is edge orientation
    para_sc.blur_r	= 0.2;
    para_sc.blur_t	= 1.0;
    para_sc.blur_o  = 0.2;
    %------- blurring method, take value of ['gaussian'|'cosine']
    para_sc.blur_method = 'gaussian';
    %------- threshold used to prune nearly-zero shape context vector, will
    %be set up on the fly of generating codebook
    para_sc.sum_total_thresh   = 5;
    %------- grid sampling rate on codebook images
    para_sc.codebook_sample_step = 5;
    
    para    = para_sc;
    return;
else
    para_sc = codebook.para;
end

%%%%%%%%%%%%%%(II) parameter for feature extraction and matching%%%%%%%%%%%
% the following parameters are for feature extraction and matcing
%

if(exist('ratio','var') && ratio>1)
    error('ratio cannot be bigger than 1.0');
end

model_height    = para_sc.model_height;

%+++++++++ this ratio is used for resizing the image to achieve
%detection across scale, it should be a value less than 1.0, otherwise the
%image dimension would increase infinitely.
para_fea.ratio  = ratio;
%+++++++++ the rough aspect ratio of object height/width, for pedestrain
para_fea.asp_ratio  = 2.1; 
%--------- whether use mask function or not when do the maching, the mask
%funcion comes from model mask
para_fea.mask_fcn   = 1;
%---------  best K matches to be found for each test point
para_fea.K          = 8;

%---------  point sampling rate on test image, I want sample at least 18
%points along the height dimension
if(para_fea.asp_ratio>1)
    max_dim = model_height;
else
    max_dim = model_height/para_fea.asp_ratio;
end
para_fea.sample_step= [round(max_dim/18),round(max_dim/18)];


%%%%%%%%%%%%%%(II) parameter for voting procedure%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following parameters are for voting process, most of the parameters
% here are derived from the model height
%
%--------minimum score a voter should have, for chi2 distance, this is a
%reasonable value
para_vote.vote_thresh   = 0.6;
%--------radius used for tracing back voters
para_vote.vote_disc_rad = max_dim*0.053; %8 %15;

%--------the minimus number of votes a hypotheses should have
para_vote.min_vote      = (model_height/para_fea.sample_step(2)) *...
    (model_height/para_fea.asp_ratio/para_fea.sample_step(1))* 0.2;
%--------the number of iteration to cluster hypothesis
para_vote.nb_iter       = 2;
%------- radii used to cluster hypotheses centers
para_vote.elps_ab       = [(model_height/para_fea.asp_ratio)*0.15, model_height*0.15];
%------- the radius used to extract mask on codebook images
para_vote.maskRadius	= max(para_fea.sample_step);
%------- minimus height a 150-pixels height hypothesis should have
para_vote.min_height    = model_height*(1-0.3);
%------- maximum height a 150-pixels height hypothesis should have
para_vote.max_height    = model_height*(1+0.3);
%------- maximum voting position offset, this offset can avoid voting
%outside the image dimension
para_vote.vote_offset   = max(max(abs(codebook.relpos))) + 5;
voter_filter            = compute_smth_kernel(model_height, para_fea.asp_ratio);
%------- this filter is used to collect votes from voters
para_vote.voter_filter  = voter_filter{1};

para{1} = para_sc;
para{2} = para_fea;
para{3} = para_vote;
