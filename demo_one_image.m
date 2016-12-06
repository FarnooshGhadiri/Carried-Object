%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demo file to Generate hypo for a pedestrian% % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
global Img
global Center_Msk
global PN


fullfile(pwd,'lib')
addpath(genpath('/home/farnoosh/Projects/CO_Detector/'))
%[img,Mask]=Read_Img_Msk(Img);
dataroot    = getDataRoot;

Path_test=[dataroot,'/data/Test/Test_img_PETs.mat'];
load(Path_test)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cb_file     = fullfile(dataroot,'/codebook/codebookMR.mat'); %cb_pb_height_150_bin_30.mat');
for i=1:length(annot_list)
    
    i
    
img=annot_list(1,i).img;
Obj_Msk=annot_list(1,i).mask;
Img=annot_list(1,i).name;
Center_Msk=annot_list(1,i).Center;
PN=annot_list(1,i).PersonNum;
%%% 

edge_file   = [dataroot,'/data/edge/edge_PETs/Edg_',Img,num2str(PN),'.mat'];

load(cb_file);

verbose = 3;

ratio       = 1/1.2;
para        = set_parameter(codebook,ratio);

%%% % % %% % Detect Contour using Contour detection and Image segmentaion by Michael Randolph Maire



if(~exist(edge_file,'file'))
    if(verbose>1)
        fprintf(1,'Begin edge detection...');
        
    end
   
    obj_I_im = double(img) / 255;
     [obj_edge,obj_theta] = multiscalePb(obj_I_im);
    I_edge.edge=obj_edge;
    I_edge.theta=obj_theta;
    
    if(verbose>1)
        fprintf(1,'Edge detection: %f secs\n',toc);
    end
    % save this file so that you do not need to compute edge again in
    % the later tuning
    save(edge_file,'I_edge','img','ratio');
else
    load(edge_file);
       
end
%-----------------% % % Detect contour using Learning to Detect Natural Image Boundaries...by
% David ..



% if(~exist(edge_file,'file'))
%     if(verbose>1)
%         fprintf(1,'Begin edge detection...');
%         tic;
%     end
%    
%     I_edge  = compute_edge_pyramid(img, para{1}.detector,...
%         para{3}.min_height, para{2}.ratio,Obj_Msk); 
%     if(verbose>1)
%         fprintf(1,'Edge detection: %f secs\n',toc);
%     end
%     % save this file so that you do not need to compute edge again in
%     % the later tuning
%     save(edge_file,'I_edge','img','ratio');
% else
%     load(edge_file);
%        
% end


para{2}.ratio   = ratio;

[Hypo_Prob] = sc_detector(Obj_Msk,codebook, I_edge , para, verbose);
save([dataroot,'/Hypo/Hypo_PETs/',Img,'_',num2str(PN), '.mat'],'Hypo_Prob');

end


