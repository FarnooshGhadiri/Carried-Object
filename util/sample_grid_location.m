function locations = sample_grid_location(imgsz, sample_step)
% locations = sample_grid_location(imgsz, sample_step)
% sampling locations in the image 
% INPUT:
%   IMGSZ:  size of image
% 	SAMPLE_STEP: sampling step
% OUTPUT:
%   LOCATIONS : [#pos x 2], sampled positions nb_pos *2, each row [x, y]
% 
%   Liming Wang, Jan 2008
%

if(length(sample_step)==1)
    sample_step = [sample_step,sample_step];
end

imgh = imgsz(1);
imgw = imgsz(2);

[x,y]= meshgrid(3:sample_step(1):imgw,3:sample_step(2):imgh);


locations = [x(:),y(:)];
