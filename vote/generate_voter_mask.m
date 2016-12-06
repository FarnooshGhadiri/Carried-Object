function [Temp] = generate_voter_mask(vote_record,scoresK,scoresK_id,...
        testpos, codebook, valid_vote_idx, imgsz, maskRadius)
%[retVoterMask] = generate_voter_mask(vote_record,scoresK,scoresK_id,
%   testpos, codebook, valid_vote_idx, imgsz, maskRadius);
% 1. look up the 'vote_record', know which codebook entry voted
% 2. 'testpos' tells where does the codebook entry matched
% 3. 'scoresK_id' tells which entry it is
% 4. 'codebook' tells the image,code ids information
% 5. 'annolist' tells the mask and image information
%
% to obtain the mask infrmation, we first check the 'vote_record' to find the
% positive voter, then check the 'scoresK_id' to see which codebook entry
% gives the vote, then check the 'codebook' to obtain the object and image id
% then go back to 'annolist' to get the mask information and then put the
% mask information at 'testpos'
% 
%

if(length(maskRadius)==1)
    maskRadius  = [maskRadius,maskRadius];
end

annolist    = codebook.annolist;
nb_test     = size(testpos, 1);
nhypo       = length(vote_record);
retVoterMask= zeros(imgsz(1),imgsz(2));

for hypoIdx = 1:nhypo    
    idx_code_curhypo= vote_record(hypoIdx).voter_id;
    idx_code_curhypo= valid_vote_idx(idx_code_curhypo);
    idx_fea_curhypo = mod(idx_code_curhypo-1, nb_test) + 1; % feature id
    x_pos = testpos(idx_fea_curhypo, 1);
    y_pos = testpos(idx_fea_curhypo, 2);
    codeIdx = scoresK_id(idx_code_curhypo);
    
    hypoMask= zeros(imgsz+2*[maskRadius(2),maskRadius(1)]);
    sx  = maskRadius;
    sy  = maskRadius;
    Temp=zeros(size(hypoMask));
    for codeEntry = 1:length(codeIdx);
        codeIdxCur=codeIdx(codeEntry);
        c_y = y_pos(codeEntry);
        c_x = x_pos(codeEntry);

        img_id  = codebook.img_id(codeIdxCur);
        obj_id  = codebook.view(codeIdxCur);
        loc     = codebook.location(codeIdxCur,:);
        patch_score = scoresK(idx_code_curhypo(codeEntry));
       objMask = getObjMask(annolist,img_id,obj_id,loc,maskRadius);
       objMask=objMask>0;
       
        hypoMask= appendDisc(hypoMask,c_x,c_y,sx,sy,maskRadius,patch_score,objMask);
Temp=max(Temp,hypoMask);
    end
   
    retVoterMask_T= hypoMask(sy+1:sy+imgsz(1),sx+1:sx+imgsz(2));
    winMtx_mask = retVoterMask_T>eps;
    S_point=round((1/3)*(size(winMtx_mask,1)));
    E_point=round((2/3)*(size(winMtx_mask,1)));
    winMtx_mask_temp=imcrop(winMtx_mask,[1 S_point size(winMtx_mask,2) E_point-S_point]);
   [winMtx_mask2,nb_clr] = bwlabel(winMtx_mask_temp);
    
   winMtx_mask_temp = dropSmallRegion(winMtx_mask2,nb_clr,365);
   winMtx_mask(S_point:E_point,:)=winMtx_mask_temp;
    
    retVoterMask = retVoterMask | winMtx_mask;
    Temp=Temp(sy+1:sy+imgsz(1),sx+1:sx+imgsz(2)).*retVoterMask;
    
end




function objMask = getObjMask(annolist,img_id,obj_id,loc,maskRadius)

masksz=size(annolist(img_id).mask);

mask_tmp = zeros(masksz(1:2) + 2*[maskRadius(2),maskRadius(1)]);

sx = maskRadius(1);
sy = maskRadius(2);

mask_tmp(sy+1:sy+masksz(1),sx+1:sx+masksz(2)) = annolist(img_id).mask(:,:,1);

objMask=mask_tmp(sy+loc(2)-maskRadius(2):sy+loc(2)+maskRadius(2),...
    sx+loc(1)-maskRadius(1):sx+loc(1)+maskRadius(1));


function VR = appendDisc(VR,cx,cy,sx,sy,disc_rad,add_value,DISC)
%
rx = round(sx+cx);
ry = round(sy+cy);
if(~exist('DISC','var'))
    [X,Y]=meshgrid(-disc_rad:disc_rad);
    DISC = sqrt(X.^2+Y.^2);
    DISC = (DISC<=disc_rad);
end
DISC = DISC*add_value;

VR((ry-disc_rad):(ry+disc_rad),(rx-disc_rad):(rx+disc_rad))=...
   max(VR((ry-disc_rad):(ry+disc_rad),(rx-disc_rad):(rx+disc_rad)),DISC);