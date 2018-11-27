% Automatically generate a number datasets of the continuous-time system
% with low sampling frequencies.



close all
clear all

% basic settings
n = 24;             % dimension of system (A)
m = 40;             % length of time series
dataIDList = [48 49];
fpath = './Data';

% data generations
fprintf('Generation of data sets:\n')
for dataID = dataIDList
    spid_sys_sim(n, m, dataID, fpath);
    fprintf(['    ... No.' num2str(dataID) ', done\n'])
end
fprintf('Finished.\n')
