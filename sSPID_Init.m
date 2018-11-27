% This script is to set search paths.

% Copyright [2017] <oracleyue>
% Last modified on 31 Jan 2018



% project root path
libRootPath = [fileparts(mfilename('fullpath')) '/'];
addpath(libRootPath);

% functions
addpath([libRootPath 'functions/']);

% member functions
addpath([libRootPath 'functions/members/']);


% -------------------------------------------------
% Include external toolboxes and scripts
% -------------------------------------------------

% toolbox for matrix functions
addpath([libRootPath 'supports/mftoolbox/']);

% auxiliary scripts
addpath([libRootPath 'supports/others/']);