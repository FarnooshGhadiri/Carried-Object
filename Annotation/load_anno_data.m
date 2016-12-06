function anno_list  = load_ann_data(imgdir,para_sc)
%anno_list  = load_ann_data(imgdir,mask_dir)
%INPUT:
%   IMGDIR:     the directory contains the training images.
%   MASK_DIR:   the directory contains mask data, the mask data is also
%   stored in image format. 
%       Notes:  1) Mask image should have the same name as the
%   image name and is a gray image.
%               2) A mask image could have multiple object
%   in it, each object has a different value as its mask.
%               3) zeros value is background
%OUTPUT:
%   ANNO_LIST:  A 2D structure list which has two fields: {'img','mask'}
 
anno_list   = [];

Num_Item=0;
for ff=1:8
    Data_Path=[imgdir,num2str(ff),'/Train/'];
    Mask_Path=[imgdir,num2str(ff),'/Mask/'];
    Annot_Path=[imgdir,num2str(ff),'/Annot/'];
    srcFiles=dir(Data_Path);
    for i=3:length(srcFiles)
        
    filename=strcat(Data_Path,srcFiles(i).name);
    img=imread(filename);
    image_name=srcFiles(i).name(1:end-4);
    image_name
    [fd,syserrmsg]=fopen([Annot_Path,image_name,'.txt'],'r');
    bbox=[];
     while ~feof(fd)
            tline = fgetl(fd);
            tspace = strfind(tline,' ');
             x = str2num(tline(1:tspace(1)-1));
             y = str2num(tline(tspace(1)+1:tspace(2)-1));
             w = str2num(tline(tspace(2)+1:tspace(3)-1));
             h = str2num(tline(tspace(3)+1:tspace(4)-1));
             bbox=[bbox;x y w h];
     end
     bbox=round(bbox);
     Hight_Rate=para_sc.model_height/bbox(1,4);
     Width_Rate=para_sc.model_width/bbox(1,3);
     ROI_img=imcrop(img, bbox(2,:));
     [R,C,D]=size(ROI_img);
 F_img=imresize(ROI_img,[Hight_Rate*R Width_Rate*C]);

% -------------------------------------

% % % % % % % Praper Person Mask without carried object % % % % % %
Person_Mask=imread([Mask_Path,srcFiles(i).name]);
Person_Mask=mat2gray(Person_Mask)>0.7;
se=strel('disk',2);
Person_Mask=imdilate(Person_Mask,se);
se=strel('disk',4);

    Sub_Mask_Path=dir([Mask_Path,'bag*',srcFiles(i).name]);
    for k=1:length(Sub_Mask_Path)
        filename2=strcat(Mask_Path,Sub_Mask_Path(k).name);
        temp=imread(filename2);
        temp=mat2gray(temp)>0.7;
        temp=imdilate(temp,se);
        Person_Mask=(1-temp).*Person_Mask;
    end
   % figure(2);imshow(F_img)
     ROI_msk=imcrop(Person_Mask, bbox(2,:));
        

     ROI_msk=mat2gray(ROI_msk)>0.7;
     
     F_msk=imresize(ROI_msk,[Hight_Rate*R Width_Rate*C]);
      %figure(1);imshow(F_msk)
  
[yi,xi]=find(F_msk==1);
yc=mean2(yi);
xc=mean2(xi);

%% Calculating centroid of the person by averaging foreground coordinate 
Mask_Center=[xc yc];

    Num_Item=Num_Item+1;
        anno_list(Num_Item).img    = F_img;
        anno_list(Num_Item).mask   = F_msk;
        anno_list(Num_Item).view=ff;
        anno_list(Num_Item).Center=Mask_Center;
        

    end
end
 

