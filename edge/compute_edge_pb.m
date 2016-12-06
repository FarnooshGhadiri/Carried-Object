function [ pb, theta ] = compute_edge_pb( img )
%[ pb, theta ] = compute_edge_pb( img )
%
%
%   Liming Wang, Jan 2008
%

if(exist('pbCGTG.m','file')~=2)
    disp('You got this error because you have no pb edge detector installed.');
    disp('Please go to ');
    disp('http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/segbench/');
    disp(' and download the latest the tar ball and install it before running');
    disp(' this script again.');
    error('--------------------------------------------------------------');
end

img = im2double(img);

img = rescaleImage(img);


if (size(img, 3) > 1)
    [pb, theta] = pbCGTG(img);
else
    [pb, theta] = pbBGTG(img);
end

