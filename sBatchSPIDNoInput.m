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



clear all; close all

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
dataIDList = 1:50;
dataPath = './Data/';
resPath = './Results/';

fid = fopen([resPath 'summary.log'], 'w');
fprintf(fid, 'Data ID\tTP\tFP\tFN\tTN\tlambda\tType I\tType II\n');
fspec = '%d\t%d\t%d\t%d\t%d\t%.2e\t%.2f\t%.2f\n';

%% Perform System Id

% Start timer
netTimer = tic;

fprintf('Network reconstruction start:\n')
for dataID = dataIDList
    dataTimer = tic;
    dataName = ['spfreq_p24_N40_ID' num2str(dataID)];

    % call sysid function
    warning('off', 'MATLAB:logm:nonPosRealEig');
    [resName, perfMetric] = ...
        spid_sysid_noinput(dataName, dataPath, resPath, spidOptions);

    % print performance results
    fprintf(fid, fspec, [dataID perfMetric]);

    % export figures
    spid_plot(1, resName, 'fname', [resPath dataName '_A'   ], 'visible', 'off');
    spid_plot(2, resName, 'fname', [resPath dataName '_Ad'  ], 'visible', 'off');
    spid_plot(4, resName, 'fname', [resPath dataName '_Conv'], 'visible', 'off');
    spid_plot(5, resName, 'fname', [resPath dataName '_ROC' ], 'visible', 'off');
    spid_plot(6, resName, 'fname', [resPath dataName '_PR'  ], 'visible', 'off');

    dataEtime = toc(dataTimer);
    fprintf(['    ... dataset ' num2str(dataID) ': done; ' ...
             'process time: %f sec\n'], dataEtime);
end
fclose(fid);
fprintf('All datasets processed.\n\n')

% End timer
netEtime = toc(netTimer);
fprintf('Elapsed time (network inference) is %f seconds \n', ...
        netEtime);
