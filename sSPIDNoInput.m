% Sparse network reconstruction using full-state measurement time-series
% with low sampling frequencies.
% Algorithm: constrained Gauss-Netwon Method
%
% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 16 May 2017



clear all
close all

%% Solver Settings
spidOptions.optvar_precision = 1e-8;   % |A_k+1 - A_k|_2 <= ...
spidOptions.optval_precision = 1e-12;  % |fA_k+1 - fA_k|_2 <= ...
spidOptions.ls_precision = 1e-12;  % ls - line search
spidOptions.ls_alpha = .2;   % default: <1/2
spidOptions.iteration_loops = 50;
spidOptions.init_gamma = 1;  % choose initial point in nonlinear lsq
spidOptions.sparsity_method = 'lasso';

spidOptions.lambda = .05;
spidOptions.zero = 1e-6;


%% Load Datasets
% dataset_sim gives X1, X2, A_gt, Ad_gt, m, n, Ts
dataID = 4;
dataName = ['spfreq_p24_N40_ID' num2str(dataID)];
dataPath = './Data/';
resPath = './Results/';


%% Perform System Id

% Start timer
netTimer = tic;

% Call sysid function
resName = spid_sysid_noinput(dataName, dataPath, resPath, spidOptions);

% End timer
netEtime = toc(netTimer);
fprintf('Elapsed time (network inference) is %f seconds \n', ...
        netEtime);


%% Plot results
spid_plot(1, resName, 'fname', [resPath dataName '_A'   ]);
spid_plot(2, resName, 'fname', [resPath dataName '_Ad'  ]);
spid_plot(4, resName, 'fname', [resPath dataName '_Conv']);
spid_plot(5, resName, 'fname', [resPath dataName '_ROC' ]);
spid_plot(6, resName, 'fname', [resPath dataName '_PR'  ]);
