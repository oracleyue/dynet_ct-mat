% This script is to produce averaged ROC and P-R curves from network
% reconstruction results.

% Copyright [2018] <oracleyue>
% Last modified on 31 Jan 2018


clear; close all

%% Load Data

dataIDList = 1:50;
resPath = './Results/';
prefix = 'spfreq_p24_N40_ID';
postfix = '_results';

for dataID = dataIDList
    resName = [resPath prefix num2str(dataID) postfix];
    load(resName);

    % pack up results
    vecAgt   = reshape(A_gt, [],1);
    vecAk    = reshape(A_k, [],1);
    vecAlogm = reshape(A_logm, [],1);

    cAgt(:,dataID) = vecAgt;
    cAk(:,dataID) = vecAk;
    cAlogm(:,dataID) = vecAlogm;

    % clear the rest
    clear -regexp '^(?!resPath|prefix|postfix|dataIDList|c).*';
end


%% Performance Curves

% ----------------------------------------------------------------
% ROC curve
% ----------------------------------------------------------------
for dataID = dataIDList
    vecAgt   = cAgt(:,dataID);
    vecAk    = cAk(:,dataID);
    vecAlogm = cAlogm(:,dataID);

    labels = {};
    for i = 1:length(vecAgt)
        if vecAgt(i) == 0, labels = [labels 'zero'];
        else labels = [labels 'nonzero']; end
    end

    XVals = 0:.01:1;
    [xAk, yAk, tAk, aucAk] = ...
        perfcurve(labels, abs(vecAk), 'nonzero', 'XVals', XVals);
    [xAlogm, yAlogm, tAlogm, aucAlogm] = ...
        perfcurve(labels, abs(vecAlogm), 'nonzero', 'XVals', XVals);

    xValsLength = length(XVals) + 1;
    xAk = mAppendOnes(xAk, xValsLength);
    yAk = mAppendOnes(yAk, xValsLength);
    tAk = mAppendOnes(tAk, xValsLength, 0);
    xAlogm = mAppendOnes(xAlogm, xValsLength);
    yAlogm = mAppendOnes(yAlogm, xValsLength);
    tAlogm = mAppendOnes(tAlogm, xValsLength, 0);

    cxAk(:, dataID)      = xAk;
    cyAk(:, dataID)      = yAk;
    ctAk(:, dataID)      = tAk;
    caucAk(:, dataID)    = aucAk;
    cxAlogm(:, dataID)   = xAlogm;
    cyAlogm(:, dataID)   = yAlogm;
    ctAlogm(:, dataID)   = tAlogm;
    caucAlogm(:, dataID) = aucAlogm;
end

% compute statistics
xAk_mean = mean(cxAk, 2);
yAk_mean = mean(cyAk, 2);
tAk_mean = mean(ctAk, 2);
aucAk_mean = mean(caucAk, 2);

xAk_std = std(cxAk, 0, 2);
yAk_std = std(cyAk, 0, 2);
tAk_std = std(ctAk, 0, 2);
aucAk_std = std(caucAk, 0, 2);

xAlogm_mean = mean(cxAlogm, 2);
yAlogm_mean = mean(cyAlogm, 2);
tAlogm_mean = mean(ctAlogm, 2);
aucAlogm_mean = mean(caucAlogm, 2);

xAlogm_std = std(cxAlogm, 0, 2);
yAlogm_std = std(cyAlogm, 0, 2);
tAlogm_std = std(ctAlogm, 0, 2);
aucAlogm_std = std(caucAlogm, 0, 2);

% plot roc curves
fig_hl = figure;   % ('visible','off');
figax1 = plot(xAk_mean, yAk_mean);
figax1.LineWidth = 1.7;
hold on
hFill = fill([xAk_mean; flipud(xAk_mean)], ...
             [yAk_mean + yAk_std; flipud(yAk_mean - yAk_std)], ...
             'k', 'LineStyle', 'none');
set(hFill,'FaceAlpha',0.1)

figax2 = plot(xAlogm_mean, yAlogm_mean);
figax2.LineWidth = 1.7;
hFill= fill([xAlogm_mean; flipud(xAlogm_mean)], ...
            [yAlogm_mean + yAlogm_std; flipud(yAlogm_mean - yAlogm_std)], ...
            'k', 'LineStyle', 'none');
set(hFill,'FaceAlpha',0.15)

ax = gca;
ax.FontSize = 14;
xlim([0 1]); ylim([0 1]);
legend([figax1 figax2], {'regularized', 'principal logarithm'},...
       'Location','SouthEast');
xlabel('False positive rate'); ylabel('True positive rate');
title(['ROC Curve for DNR (AUC: ' num2str(aucAk_mean) ')'])
hold off
set(fig_hl,'Units','Inches', 'position',[6 4.8611 5.5 4.3]);

save_fig('./results/ROCcurve-avg.pdf')

% ----------------------------------------------------------------
%% Prec-Recall curve
% ----------------------------------------------------------------
cxAk = []; cyAk = []; ctAk = []; caucAk = [];
cxAlogm = []; cyAlogm = []; ctAlogm = []; caucAlogm = [];

for dataID = dataIDList
    vecAgt   = cAgt(:,dataID);
    vecAk    = cAk(:,dataID);
    vecAlogm = cAlogm(:,dataID);

    labels = {};
    for i = 1:length(vecAgt)
        if vecAgt(i) == 0, labels = [labels 'zero'];
        else labels = [labels 'nonzero']; end
    end

    XVals = 0:.05:1;
    [xAk, yAk, tAk, aucAk] = ...
        perfcurve(labels, abs(vecAk), 'nonzero', 'XVals', XVals,...
                  'xCrit', 'reca', 'yCrit', 'prec');
    [xAlogm, yAlogm, tAlogm, aucAlogm] = ...
        perfcurve(labels, abs(vecAlogm), 'nonzero', 'XVals', XVals,...
                  'xCrit', 'reca', 'yCrit', 'prec');

    xValsLength = length(XVals);
    xAk = mRemoveZeros(xAk, xValsLength);
    yAk = mRemoveZeros(yAk, xValsLength);
    tAk = mRemoveZeros(tAk, xValsLength);
    xAlogm = mRemoveZeros(xAlogm, xValsLength);
    yAlogm = mRemoveZeros(yAlogm, xValsLength);
    tAlogm = mRemoveZeros(tAlogm, xValsLength);

    cxAk(:, dataID)      = xAk;
    cyAk(:, dataID)      = yAk;
    ctAk(:, dataID)      = tAk;
    caucAk(:, dataID)    = aucAk;
    cxAlogm(:, dataID)   = xAlogm;
    cyAlogm(:, dataID)   = yAlogm;
    ctAlogm(:, dataID)   = tAlogm;
    caucAlogm(:, dataID) = aucAlogm;
end

% compute statistics
xAk_mean = mean(cxAk, 2);
yAk_mean = mean(cyAk, 2);
tAk_mean = mean(ctAk, 2);
aucAk_mean = mean(caucAk, 2);

xAk_std = std(cxAk, 0, 2);
yAk_std = std(cyAk, 0, 2);
tAk_std = std(ctAk, 0, 2);
aucAk_std = std(caucAk, 0, 2);

index = find(isnan(yAk_mean));
xAk_mean(index) = [];   xAk_std(index) = [];
yAk_mean(index) = [];   yAk_std(index) = [];
tAk_mean(index) = [];   tAk_std(index) = [];

xAlogm_mean = mean(cxAlogm, 2);
yAlogm_mean = mean(cyAlogm, 2);
tAlogm_mean = mean(ctAlogm, 2);
aucAlogm_mean = mean(caucAlogm, 2);

xAlogm_std = std(cxAlogm, 0, 2);
yAlogm_std = std(cyAlogm, 0, 2);
tAlogm_std = std(ctAlogm, 0, 2);
aucAlogm_std = std(caucAlogm, 0, 2);

index = find(isnan(yAlogm_mean));
xAlogm_mean(index) = [];   xAlogm_std(index) = [];
yAlogm_mean(index) = [];   yAlogm_std(index) = [];
tAlogm_mean(index) = [];   tAlogm_std(index) = [];

uncertainX = min([xAk_mean(1) xAlogm_mean(1)]);
% used to mark the region in P-R curve that is uncertain

% plot P-R curves
fig_hl = figure;   % ('visible','off');
figax1 = plot(xAk_mean, yAk_mean, '-');
figax1.LineWidth = 1.7;
hold on
hFill = fill([xAk_mean; flipud(xAk_mean)], ...
             [yAk_mean + yAk_std; flipud(yAk_mean - yAk_std)], ...
             'k', 'LineStyle', 'none');
set(hFill,'FaceAlpha', 0.1);

figax2 = plot(xAlogm_mean, yAlogm_mean, '-');
figax2.LineWidth = 1.7;
hFill = fill([xAlogm_mean; flipud(xAlogm_mean)], ...
             [yAlogm_mean + yAlogm_std; flipud(yAlogm_mean - yAlogm_std)], ...
             'k', 'LineStyle', 'none');
set(hFill,'FaceAlpha', 0.15);

colortbl = colormap('colorcube');
hFill = fill([0 uncertainX uncertainX 0], [0 0 1 1], ...
             colortbl(end-3,:), 'LineStyle', '--');
set(hFill, 'FaceAlpha', .8);

ax = gca;
ax.FontSize = 14;
xlim([0 1]); ylim([0 1]);
legend([figax1 figax2], {'regularized', 'principal logarithm'},...
       'Location','NorthEast');
xlabel('Recall'); ylabel('Precision');
title(['P-R Curve for DNR (AUC: ' num2str(aucAk_mean) ')'])
hold off
set(fig_hl,'Units','Inches', 'position',[6 4.8611 5.5 4.3]);

save_fig('./results/PRcurve-avg.pdf')
