function Image_IF = fft2_filt(Image,FilterBank,edge_method)
% F= ff2_filt(Image, FilterBank)
% Image: one channel image or gray image
% FilterBank: n by n by (#nb)  could have (#nb) filters
% Image_IF: (#nb) filtered results
%
% reference: fft_filt_2.m by Jianbo Shi
% 
% edge_method:
% 'ZERO'  which pack the outboard with 0
% 'REFLECT' using reflection
% 'EXPAND' using reflection and revert
%


if(~exist('edge_method'))
    edge_method = 'zero';
end

if(strcmp(edge_method,'zero'))
    edge_method =1;
end
if(strcmp(edge_method,'reflect'))
    edge_method =2;
end
if(strcmp(edge_method,'expand'))
    edge_method =3;
end

[fr,fc, fn] = size(FilterBank);
[imr,imc]=size(Image);
if(~isa(Image,'double'))
    Image=double(Image);
end
fwin_r = floor(fr/2);
fwin_c = floor(fc/2);
Image_expand = zeros(size(Image)+[fwin_r*2 fwin_c*2]);
Image_expand(fwin_r+1:fwin_r+imr,fwin_c+1:fwin_c+imc)=Image;
% TODO: do more edge method
if(edge_method ==2)
%     this method use reflection
%     altogether, eight regions should be taken care of
%     1) left region
    Image_expand(fwin_r+1:fwin_r+imr,1:fwin_c) =...
        fliplr(Image(:,1:fwin_c));
%     2) right region
    Image_expand(fwin_r+1:fwin_r+imr,fwin_c+1+imc:end) =...
        fliplr(Image(:,end-fwin_c+1:end));
%     3) up region
    Image_expand(1:fwin_r,fwin_c+1:fwin_c+imc) =...
        flipud(Image(1:fwin_r,:));
%     4) down region
    Image_expand(fwin_r+1+imr:end,fwin_c+1:fwin_c+imc) =...
        flipud(Image(end-fwin_r+1:end,:));
%     5) left-up region
    Image_expand(1:fwin_r,1:fwin_c) = flipud(fliplr(Image(1:fwin_r,1:fwin_c)));
%     6) right_up region
    Image_expand(1:fwin_r,fwin_c+1+imc:end) = flipud(fliplr(Image(1:fwin_r,end-fwin_c+1:end)));
%     7) left-down region
    Image_expand(fwin_r+1+imr:end,1:fwin_c)=flipud(fliplr(Image(end-fwin_r+1:end,1:fwin_c)));
%     8) right-down region
    Image_expand(fwin_r+1+imr:end,fwin_c+1+imc:end)=flipud(fliplr(Image(end-fwin_r+1:end,end-fwin_c+1:end)));    
end
if(edge_method==3)
%     this method use reflection, then revert
%     altogether, eight regions should be taken care of
%     1) left region
    Image_expand(fwin_r+1:fwin_r+imr,1:fwin_c) =...
        2*( Image(:,1)*ones(1,fwin_c))-...
        fliplr(Image(:,2:fwin_c+1));
%     2) right region
    Image_expand(fwin_r+1:fwin_r+imr,fwin_c+1+imc:end) =...
        2*( Image(:,end)*ones(1,fwin_c) )-... 
        fliplr(Image(:,end-fwin_c:end-1));
%     3) up region
    Image_expand(1:fwin_r,fwin_c+1:fwin_c+imc) =...
        2*( ones(fwin_r,1)*Image(1,:) )-...
        flipud(Image(2:fwin_r+1,:));
%     4) down region
    Image_expand(fwin_r+1+imr:end,fwin_c+1:fwin_c+imc) =...
        2*( ones(fwin_r,1)*Image(end,:) )-...
        flipud(Image(end-fwin_r:end-1,:));
%     5) left-up region
    Image_expand(1:fwin_r,1:fwin_c) = ...
        2*Image(1,1) - flipud(fliplr(Image(2:fwin_r+1,2:fwin_c+1)));
%     6) right_up region
    Image_expand(1:fwin_r,fwin_c+1+imc:end) =...
        2*Image(1,end) - flipud(fliplr(Image(2:fwin_r+1,end-fwin_c:end-1)));
%     7) left-down region
    Image_expand(fwin_r+1+imr:end,1:fwin_c)=...
        2*Image(end,1) - flipud(fliplr(Image(end-fwin_r:end-1,2:fwin_c+1)));
%     8) right-down region
    Image_expand(fwin_r+1+imr:end,fwin_c+1+imc:end)=...
        2*Image(end,end) - flipud(fliplr(Image(end-fwin_r:end-1,end-fwin_c:end-1)));   
end


Image_IF=zeros(imr,imc,fn);


%% method 1: convolve

for b=1:fn
    tmp=conv2(Image_expand,FilterBank(:,:,b),'same');
    Image_IF(:,:,b)=tmp(fwin_r+1:fwin_r+imr,fwin_c+1:fwin_c+imc);
end


%% method 2: fft2 and ifft2

% % % Image_F = fft2(Image_expand);
% % % 
% % % Image_IF=zeros(imr,imc,fn);
% % % 
% % % for b=1:fn
% % %   FB_F=FilterBank(:,:,b);
% % %   FB_F=fliplr(FB_F);
% % %   FB_F=flipud(FB_F);
% % %   tmp=zeros(imr+fr-1,imc+fc-1);
% % %   tmp(floor((imr-size(FB_F,1))/2):floor((imr-size(FB_F,1))/2)+size(FB_F,1)-1,...
% % %       floor((imc-size(FB_F,2))/2):floor((imc-size(FB_F,2))/2)+size(FB_F,2)-1);
% % %    FB_F=fft2(FB_F,imr+fr-1,imc+fc-1);
% % %    tmp = Image_F.*FB_F;
% % %    tmp = ifft2(tmp);
% % %    tmp = real(tmp);
% % %    Image_IF(:,:,b) = tmp(fwin+1:fwin+imr,fwin+1:fwin+imc);
% % % end