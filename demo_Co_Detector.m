
% for threshhold codebook 
clear
clc
addpath lib;
Person_Height=180;
addpath(genpath('/home/farnoosh/Projects/CO_Detector/'))
dataroot    = getDataRoot;

Path_seed=[dataroot,'/Hypo/Hypo_PETs/'];
Path_edge=[dataroot,'/data/edge/edge_PETs/'];
load([dataroot,'/data/Test/Test_img_PETs.mat'])

for i=1:length(annot_list)
    hh=i
    tic
    Img=annot_list(1,i).img;
    
    Mask=annot_list(1,i).mask;
    
    [R,C]=size(Mask);
    img_name=annot_list(1,i).name;
    
    PN=annot_list(1,i).PersonNum;
    edge_file   = [dataroot,'/data/edge/edge_PETs/Edg_',img_name,num2str(PN),'.mat'];
    load([Path_seed,img_name,'_',num2str(PN), '.mat']);
    Hypo=Hypo_Prob.hypo;
    Seed_Point=Hypo_Prob.seedpoint;
    Cor_img=annot_list(1,i).Corrdinate;
    
    frame1=annot_list(1,i).Frame;
    Carried_Object=[];
     
    nvec = 17;
    [EigVal, EigVect, D] = im2evec(Img,nvec,Mask,img_name,edge_file);
    
    %% biased normalized cut demo
    [nr,nc,nb] = size(Img);
   
    %% display biased ncuts for various values of gamma
    ROIMODE    = 0; % input is ROI or a set of points
    numinput   = 10; % number of input points
    NUMDISPLAY = 1; % number of display panels previously was 9
    
    
    wts = zeros(length(EigVal),1);
    m = round(sqrt(NUMDISPLAY)); n = ceil(NUMDISPLAY/m);
    Fin_Seg=zeros(size(Mask));
    %Carried_Object=zeros(size(Mask));
    Num_seg=0; % Number of detected Hypo_Prob
    for j=1: length(Seed_Point) % Nember of contour
       figure(1); clf;
        imshow(Img);
        hold on;
        cx=Seed_Point{1,j}(:,1);
        cy=Seed_Point{1,j}(:,2);
        %[cx,cy] = ginput(numinput);
        
        bw = zeros(nr,nc);
        cx = round(cx);
        cy = round(cy);
        idx_cx=find(cx>nc | cx<1);
        idx_cy=find(cy>nr | cy<1);
        cx(idx_cx)=[];
        cy(idx_cx)=[];
        cy(idx_cy)=[];
        cx(idx_cy)=[];
                

        for k = 1:length(cx),
            bw(cy(k),cx(k)) = 1;
        end
        plot(cx,cy,'b+');
         plot(cx,cy,'bo');
        
        %make the normalized seed vector
        s = D*bw(:);
      % figure(2); clf;
        
        gamma = 0;
        for ee = 2:length(EigVal) % ignore the all 1's vector
            wts(ee) = (EigVect(:,ee)'*s)/(EigVal(ee) - gamma);
        end
        WEigVect = EigVect*wts;
%         
       subplot(m,n,1);
        imagesc(reshape(WEigVect,nr,nc)); axis image; axis off;
         title(sprintf('\\gaFin_Segmma=%.2d', gamma));
% %         
        fprintf('press: n next, t toggle-mode, <-,-> decrease,increase points\n');
        Seg_res=reshape(WEigVect,nr,nc);
        Normalized_seg=mat2gray(Seg_res);
        Seg=Normalized_seg>0.9;
        Seg_area=sum(Seg(:));
        stats=regionprops(Seg,'BoundingBox');
        Seg_height=1;
        for L=1:length(stats)
            Seg_height=max(Seg_height,stats(L,1).BoundingBox(1,4));
        end
        if Seg_height<1.1*(Person_Height/2)
        
        hypo_thresh=45; %po_thresh=29
        
        Over_Lap=Seg.*(Hypo>0.7);
        Overlap_area=sum(Over_Lap(:));
        % % % Overlap Percentage % % % %
        Hypo_Overlap=(Overlap_area/Seg_area)*100;
        Overlap_back=Seg.*(1-Mask);
        Overlap_msk_count=sum(Overlap_back(:));
        Mask_Overlap=(Overlap_msk_count/Seg_area)*100;
       
               
      
        
        if  (Hypo_Overlap<hypo_thresh) & (Mask_Overlap<29) % Mask_Overlap=29
           
            Num_seg=Num_seg+1;
            Carried_Object{Num_seg}.Seg=Seg;
            Carried_Object{Num_seg}.Hypo_Overlap=Hypo_Overlap;
            Carried_Object{Num_seg}.Mask_Overlap=Mask_Overlap;
         
            
        end 
        end
    end
    toc
    kk=1;
    %jj=2;
    SE=strel('rectangle',[10,10]); 
    NumSeg=length(Carried_Object);
    while kk<=NumSeg-1
        jj=kk+1;
        while jj<=NumSeg
            Ca_kk=imdilate(Carried_Object{kk}.Seg,SE);
            Ca_jj=imdilate(Carried_Object{jj}.Seg,SE);
       Overlap_Seg =sum(sum(Ca_kk.*Ca_jj));
       Per_Overlap=(Overlap_Seg/sum(Ca_jj(:)))*100;
        Hyp_Con=(Carried_Object{kk}.Hypo_Overlap-Carried_Object{jj}.Hypo_Overlap);
            Msk_Con=Carried_Object{kk}.Mask_Overlap-Carried_Object{jj}.Mask_Overlap;
              
                        
            if  Per_Overlap>0 
             
            M_Hypo_kk=mean(Hypo(Carried_Object{kk}.Seg));
            M_Hypo_jj=mean(Hypo(Carried_Object{jj}.Seg));
            
            if abs(Hyp_Con)>15 & Carried_Object{jj}.Hypo_Overlap>=30 
                
                Carried_Object(jj)=[];
                jj=jj-1;
                NumSeg=length(Carried_Object);
                               
            else
                    if abs(Hyp_Con)>15 & Carried_Object{kk}.Hypo_Overlap>=30 
                Carried_Object(kk)=[];
                kk=kk-1;
                NumSeg=length(Carried_Object);
                break;
                
                    end
                 
           
               
                
            end
                
                       
            end
            jj=jj+1;
        end
            kk=kk+1;
            
    end  
    if ~isempty(Carried_Object)
     for h=1:length(Carried_Object)
         Fin_Seg=Fin_Seg|Carried_Object{h}.Seg;
     end
     Fin_Seg = bwareaopen(Fin_Seg, 20);%250
    Fin_Seg =imclose(Fin_Seg,SE);
    st=regionprops(Fin_Seg,'BoundingBox');
    Object_Coordinate(i).Cor=st;
    imshow(frame1)
    Mask_Boundary=zeros(size(frame1,1),size(frame1,2));
    for k=1:length(st)
        
                
        w=Cor_img(1,3);
        h=Cor_img(1,4);
        rec=[st(k,1).BoundingBox(1,1)*(1/w) st(k,1).BoundingBox(1,2)*(1/h) st(k,1).BoundingBox(1,3)*(1/w) st(k,1).BoundingBox(1,4)*(1/h)];
        rec(1,1)=rec(1,1)+Cor_img(1,1);
        rec(1,2)=rec(1,2)+Cor_img(1,2);
        st(k).BoundingBox=rec;
        [Rf,Cf,Df]=size(frame1);
        
        rectangle('position', rec, 'edgecolor', 'b');
        % draw Line around object
        rec=floor(rec);
        Mask_Boundary(rec(2)-2:rec(2)+rec(4)-2,rec(1):rec(1)+1)=1;
        Mask_Boundary(rec(2)-2:rec(2)+rec(4)-2,rec(1)+rec(3):rec(1)+rec(3)+1)=1;
        Mask_Boundary(rec(2)-1:rec(2),rec(1):rec(1)+rec(3))=1;
        Mask_Boundary(rec(2)+rec(4)-3:rec(2)+rec(4)-2,rec(1):rec(1)+rec(3)+1)=1;
        RGB=imoverlay(frame1,Mask_Boundary,[0,1,0]);
    end
    Cor_img=floor(Cor_img);
    im=imcrop(RGB,[Cor_img(1,1)  Cor_img(1,2)  size(Fin_Seg,2)*(1/w) size(Fin_Seg,1)*(1/h)]);
    imwrite(im,['/home/farnoosh/Projects/CO_Detector/Result/rec_',num2str(hh),'.png'])       
    Object_Coordinate(i).Cor=st;
    if ~isempty(st)
    Msk_Seg=imresize(Fin_Seg,[R*1/h C*1/w]);
        [R,C]=size(Msk_Seg);
        Mask_F=zeros(Rf,Cf);
        Mask_F(Cor_img(1,2):R+Cor_img(1,2)-1,Cor_img(1,1):C+Cor_img(1,1)-1)=Msk_Seg;
       Mask_F= imresize(Mask_F,[size(frame1,1),size(frame1,2)]);
        RGB=imoverlay(frame1,Mask_F,[0,1,0]);
        im=imcrop(RGB,[Cor_img(1,1)  Cor_img(1,2)  size(Fin_Seg,2)*(1/w) size(Fin_Seg,1)*(1/h)]);
        imwrite(im,['/home/farnoosh/Projects/CO_Detector/Result/Seg_',num2str(hh),'.png'])       
    else
         Object_Coordinate(i).Cor=[];
    end
     
      
    else
        Object_Coordinate(i).Cor=[];
    end
     
end
save('/home/farnoosh/Projects/BSR/grouping/Carried_Object/Co_Hypo9_th.mat','Object_Coordinate');
