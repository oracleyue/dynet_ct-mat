function spid_sys_sim(n, m, dataID, fpath)
% SPID_SYS_SIM Simulate the continuous-time system to generate data for system
% identification.
%
% Inputs:
%   - n      : dimension of state variable (i.e. dim(A))
%   - m      : length of time series
%   - dataID : label used in prefixing file names
%   - fpath  : path to save datasets
%
% Example:
%   spid_sys_sim(24, 40, 0, './datasets/')

% Copyright [2017] <oracleyue>
% Last modified on 30 Jan 2018



%% --- arguments parsing ---
if nargin == 3
    fpath = './datasets/';
end


%% --- systems generation ---

% name prefix of data file
fname_prefix = ['p' num2str(n) '_' 'N' num2str(m) '_ID' num2str(dataID)];
if fpath(end) ~= '/', fpath = [fpath '/']; end
fname = [fpath 'spfreq_' fname_prefix '.mat'];

% generate sparse stable A
[A, fig_hlA, fig_hlNet] = spid_sprandstab(n, 'nicolo-overlap', .05, 8);

% determine Ts (low sampling frequency)
imagEigA = abs(imag(eig(A)));
[~,~,imagEigA] = find(imagEigA);
Ts = .5*min(pi./imagEigA);


%% --- system simulations ---
C = eye(n);
K = eye(n);
sys_ss = idss(A, [], C, [], K, 'ts', 0);
init_val = randn(n,1)*1;   % *1
noise = randn(m, n)*.1;   % *.01
opt = simOptions('AddNoise', true, 'NoiseData', noise, ...
                 'InitialCondition',init_val);

output_sim = sim(sys_ss, iddata([], zeros(m,0), Ts), opt);


%% --- format and export datasets ---
dataset = output_sim.y';
X1 = dataset(:, 2:end);
X2 = dataset(:, 1:end-1);
sys_dim = n;        % dimension of ss sys, i.e. A
ts_dim = m;         % length of time series
Ad_gt = expm(Ts*A); % groud truth of Ad
A_gt = A;           % groud truth of A

save(fname, 'X1', 'X2', 'sys_dim', 'ts_dim', ...
    'A_gt', 'Ad_gt', 'Ts', 'output_sim');


%% --- plot time series ---
fig_hlSig = figure('visible','off');
set(fig_hlSig,'Units','Inches', 'position',[2.1528 3.2639 16.9861 9.1111]);
time = (1:m)*Ts;
for i = 1:n
    subplot(4,ceil(n/4),i)
    plot(time, output_sim.y(:,i), '*-')
end


%% --- export plot of signals in PDF ---
figname = [fpath 'figures/spfreq_' fname_prefix '_signals' '.pdf'];
pos = get(fig_hlSig,'Position');
set(fig_hlSig,'PaperPositionMode','Auto','PaperUnits','Inches',...
              'PaperSize',[pos(3), pos(4)])
print(fig_hlSig, figname,'-dpdf','-r0')

figname = [fpath 'figures/spfreq_' fname_prefix '_A' '.pdf'];
pos = get(fig_hlA,'Position');
set(fig_hlA,'PaperPositionMode','Auto','PaperUnits','Inches',...
            'PaperSize',[pos(3), pos(4)])
print(fig_hlA, figname,'-dpdf','-r0')

figname = [fpath 'figures/spfreq_' fname_prefix '_network' '.pdf'];
pos = get(fig_hlNet,'Position');
set(fig_hlNet,'PaperPositionMode','Auto','PaperUnits','Inches',...
              'PaperSize',[pos(3), pos(4)])
print(fig_hlNet, figname,'-dpdf','-r0')