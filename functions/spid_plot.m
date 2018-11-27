function spid_plot(figNum, figResMat, varargin)
% SPID_PLOT Plot figures from sparse continuous-time identifications
% spid_plot(figNum) plot the "figNum"-th figure.
%
% INPUT:
%   figNum    : an integer; refer to source codes for which type of figure
%                  it indicates
%   figResMat : full path of .mat file that saves results
%
%   Name-Value Pairs:
%       - 'figName':string    : set the name for PDF export
%       - 'pdf'    :{0,1}     : save the figure in PDF; default 0
%       - 'tikz'   :{0,1}     : set 1 to reprocess the figure via TikZ in LaTeX
%       - 'visible':'off'     : disable showing figures
%
% Examples:
%   spid_plot(4, './workspace_figplot.mat')
%   spid_plot(4, './workspace_figplot.mat', 'fname', './results/output')
%   spid_plot(4, './workspace_figplot.mat', 'tikz', 1, 'fname', './results/output')
%

% Copyright (c) 2015-2016, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Modified on 16 May 2017



%% Argument Parsing
parser = inputParser;
isFigNum = @(x) isnumeric(x) && isscalar(x) && rem(x,1)==0;
addRequired(parser, 'figNum', isFigNum);
addRequired(parser, 'figResMat', @ischar);
isBool = @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1);
addParameter(parser, 'pdf', 1, isBool);
addParameter(parser, 'tikz', 0, isBool);
addParameter(parser, 'fname', './results/output', @ischar);
addParameter(parser, 'visible', 'on', @ischar);
parse(parser, figNum, figResMat, varargin{:});
figNum = parser.Results.figNum;
figResMat = parser.Results.figResMat;
pdfFlag = parser.Results.pdf;
tikzFlag = parser.Results.tikz;
fname = parser.Results.fname;
visible = parser.Results.visible;
if strcmp(visible, 'off'), showFig = 0;
else  showFig = 1; end


%% Load Result .mat File
load(figResMat);


%% Figure Plotting

% Default figure settings
% set(0, 'DefaultAxesFontSize', 16);

% Load Datasets
load(figResMat);

% Plot the 'figNum'-th figure
switch figNum

  case 1 % matrix A's
    if showFig, fig_hl = figure(figNum);
    else        fig_hl = figure('visible','off'); end
    set(fig_hl, 'Units', 'inches', 'Position', [0.7222 6.8472 11 3.0556]);
    upval = max([max(max(abs(A_gt))), max(max(abs(A_k)))]);%,...
                 %max(max(abs(A_logm)))]);
    clims = [-upval, upval];
    subplot(1,3,1)
    imagesc(real(A_gt), clims)
    ax = gca;
    ax.FontSize = 14;
    title('$A$ (ground truth)','Interpreter','latex')

    subplot(1,3,2)
    imagesc(real(A_k), clims)
    ax = gca;
    ax.FontSize = 14;
    title('$\hat{A}_\mathrm{alg}$','Interpreter','latex')

    subplot(1,3,3)
    imagesc(real(A_logm), clims)
    ax = gca;
    ax.FontSize = 14;
    title('$\hat{A}_\mathrm{logm}$','Interpreter','latex')

    c_hl = colorbar;
    ymap = redblue(160);
    colormap(ymap);
    set(c_hl, 'Units', 'normalized', 'Position', [0.9323 0.1111 0.0192 0.8138]);


  case 2 % matrix Ad's
    if showFig, fig_hl = figure(figNum);
    else        fig_hl = figure('visible','off'); end
    set(fig_hl, 'Units', 'inches', 'Position', [0.1111 8.1250 10.9861 3.0417]);
    upval = max([max(max(abs(Ad_gt))), max(max(abs(Ak_d))),...
                 max(max(abs(Ad)))]);
    clims = [-upval, upval];
    subplot(1,3,1)
    imagesc(real(Ad_gt), clims)
    ax = gca;
    ax.FontSize = 14;
    title('$A_d$ (ground truth)','Interpreter','latex')
    subplot(1,3,2)
    imagesc(real(Ak_d), clims)
    ax = gca;
    ax.FontSize = 14;
    title('$\hat{A}_{d,\mathrm{alg}}$','Interpreter','latex')
    subplot(1,3,3)
    imagesc(real(Ad), clims)
    ax = gca;
    ax.FontSize = 14;
    title('$\hat{A}_d$','Interpreter','latex')

    c_hl = colorbar;
    ymap = redblue(160);
    colormap(ymap);
    set(c_hl, 'Units', 'normalized', 'Position', [0.9323 0.1111 0.0192 0.8138]);


  case 3  % plot time series
    if showFig, fig_hl = figure(figNum);
    else        fig_hl = figure('visible','off'); end
    [n,m] = size(output_sim);
    Ts = output_sim.Ts;
    set(fig_hl,'Units','Inches', 'position',[2.1528 3.2639 16.9861 9.1111]);
    time = (1:m)*Ts;
    for i = 1:n
        subplot(4,ceil(n/4),i)
        plot(time, output_sim.y(:,i), '*-')
        xlabel('Time (s)');
        ylabel(['x_{' num2str(i) '}(t)']);
    end

  case 4 % convergency curves
    if showFig, fig_hl = figure(figNum);
    else        fig_hl = figure('visible','off'); end
    set(fig_hl, 'Units', 'inches', 'Position', [3.4028 1.6528 5.4167 3.9444]);
    loop_ind = 0:1:length(residual_list)-1;
    figax = semilogy(loop_ind, residual_list, '*-');
    figax.LineWidth = 1.5;
    residule_gt = norm(X1 - expm(h*A_gt)*X2);
    residule_logm = norm(X1 - expm(h*A_logm)*X2);
    hold on;
    semilogy(loop_ind, ones(size(loop_ind))*residule_gt);
    semilogy(loop_ind, ones(size(loop_ind))*residule_logm);
    hold off;
    ax = gca;
    ax.FontSize = 14;
    xlabel('iteration');
    ylabel('$\|X_1 - \mathrm{exp}(h\hat{A}) X_2\|_F$','Interpreter','latex')
    legend('pred. error', 'pred. error of groundtruth', ...
           'pred. error of logm', ...
           'Location', 'best')
    grid on



  case 5  % ROC curve
    vecAgt = reshape(A_gt, [],1);
    vecAk = reshape(A_k, [],1);
    vecAlogm = reshape(A_logm, [],1);

    labels = {};
    for i = 1:length(vecAgt)
        if vecAgt(i) == 0
            labels = [labels 'zero'];
        else
            labels = [labels 'nonzero'];
        end
    end
    [xAk, yAk, tAk, aucAk] = perfcurve(labels, abs(vecAk), 'nonzero');
    [xAlogm, yAlogm, tAlogm, aucAlogm] = perfcurve(labels, abs(vecAlogm), 'nonzero');
    % plot roc curves
    if showFig, fig_hl = figure(figNum);
    else        fig_hl = figure('visible','off'); end
    figax = plot(xAk, yAk);
    figax.LineWidth = 1.7;
    hold on
    figax = plot(xAlogm, yAlogm);
    figax.LineWidth = 1.7;
    ax = gca;
    ax.FontSize = 16;
    legend('regularized','principal logarithm','Location','Best')
    xlabel('False positive rate'); ylabel('True positive rate');
    title(['ROC Curve for DNR (AUC: ' num2str(aucAk) ')'])
    hold off
    set(fig_hl,'Units','Inches', 'position',[6 4.8611 5.5 4.3]);


  case 6  % Prec-Recall curve
    vecAgt = reshape(A_gt, [],1);
    vecAk = reshape(A_k, [],1);
    vecAlogm = reshape(A_logm, [],1);

    labels = {};
    for i = 1:length(vecAgt)
        if vecAgt(i) == 0
            labels = [labels 'zero'];
        else
            labels = [labels 'nonzero'];
        end
    end
    [xAk, yAk, tAk, aucAk] = perfcurve(labels, abs(vecAk), 'nonzero',...
        'xCrit', 'reca', 'yCrit', 'prec');
    [xAlogm, yAlogm, tAlogm, aucAlogm] = perfcurve(labels, abs(vecAlogm), 'nonzero',...
        'xCrit', 'reca', 'yCrit', 'prec');
    % plot P-R curves
    if showFig, fig_hl = figure(figNum);
    else        fig_hl = figure('visible','off'); end
    figax = plot(xAk, yAk);
    figax.LineWidth = 1.7;
    hold on
    figax = plot(xAlogm, yAlogm);
    figax.LineWidth = 1.7;
    ax = gca;
    ax.FontSize = 16;
    legend('regularized','principal logarithm','Location','Best')
    xlabel('Recall'); ylabel('Precision');
    title(['P-R Curve for DNR (AUC: ' num2str(aucAk) ')'])
    hold off
    set(fig_hl,'Units','Inches', 'position',[12.6806 4.8611 5.5 4.3]);


  otherwise
    error('Your choice of "figNum" is not supported!')
end


%% Export figure in PDF with auto-resize page size
fname = regexprep(fname, '\.pdf$', '');
figName = [fname '.pdf'];

if pdfFlag
    set(fig_hl,'Units','Inches');
    pos = get(fig_hl,'Position');
    set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches',...
               'PaperSize',[pos(3), pos(4)])
    print(fig_hl, figName,'-dpdf','-r0')
end


%% Use matlab2tikz plot to export figure in TikZ and PDF
if tikzFlag
    addpath('./mat2tikz');
    matlab2tikz_run('large')
    close(fig_hl)
    system('./mat2tikz/preview.sh')
end
