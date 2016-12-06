
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %  This is a demo file to generate codebook of your own       % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/home/farnoosh/Projects/CO_Detector/'))
dataroot    = getDataRoot;
img_dir = fullfile('/home/farnoosh/Projects/Data/Train/Train&Mask2/');

codebook_file= fullfile(dataroot,'/codebook','codebookMR.mat');
 
% load default parameters
para_sc = set_parameter;

 anno_data   = load_anno_data(img_dir,para_sc);
 codebook    = load_codebook(anno_data,para_sc);
 
 save(codebook_file,'codebook');
