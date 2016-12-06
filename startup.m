function startup(rootpath)
% initialize matlab path
%
%
% Liming Wang, Jan 2008
%

if (~exist('rootpath','var'))
    rootpath = cd;
end;

% cd('..');
% cd_up   = pwd;
% cd(rootpath);
addpath(genpath(rootpath));

