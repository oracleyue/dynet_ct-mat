function [A, fig_hl1, fig_hl2] = spid_sprandstab(n, method, density, varargin)
% SPID_SPRANDSTAB generates large-scale, sparse, stable random A-matrices.
% A = SPID_SPRANDSTAB(n, 'nicolo', density) use "nicolo" method to generate
%     large-scale random sparse stable matrices; by default: "subn"= 4, "nperm" = 5.
%     Note that "density" is the sparsity density of the sub-matrices,
%     which cannot be too small (> 0.5), otherwise the sub-matrices are
%     always unstable.
% A = SPID_SPRANDSTAB(n, 'alex', density) use 'alex' method to generate
%     large-scale random sparse stable matrices; by default: "stabl_margin = 1e-2".
%     Note that "density" is the sparsity density of the whole matrix.
% A = SPID_SPRANDSTAB(n, 'sprandn', density) uses "sprandn(n)" (or "sprand(n)"
%     if use 'sprand') to generate random block-diagonal sparse stable matrices.
%     By default, "subn" = 4;
% A = SPID_SPRANDSTAB(n, 'nicolo', density, [subn nperm]) use "nicolo" method.
% A = SPID_SPRANDSTAB(n, 'alex', density, stabl_margin] use "alex" method.
% A = SPID_SPRANDSTAB(n, 'sprandn', density, subn) uses "sprandn" or "sprand".
%
% Inputs:
%   n      : scalar;
%            dimension of A (n x n)
%   method : string;
%            'nicolo', 'alex', 'sprandn', 'sprand'
%   density: scalar in (0,1];
%            - "nicolo": the sparsity density of sub-matrix of dimension "subn"
%            - "alex"  : the sparsity density of the whole matrix
%   [subn, nperm]  : (1x2) vector of integers;
%            - "subn" : the dimension of sub-matrices
%            - "nperm": the number of permutation on A to generate general
%                       sparse matrices, instead of the block-diagonal
%   subn   : integer
%   stabl_margin : doulbe
%                  a positive real number, close to zero; the stability margin
%
% Outputs:
%   A: (n x n) random matrix
%

% Interfaces
if nargin < 3
    error('Not enough arguments');
elseif nargin == 3
    if strcmp(method, 'nicolo')
        subn = 4;
        nperm = 5;
    end
    if strcmp(method, 'alex')
        stabl_margin = 1e-2;
    end
    if sum(strcmp('sprand', {'sprandn', 'sprand'}))
        subn = 4;
    end
elseif nargin == 4
    if strcmp(method, 'nicolo')
        if isscalar(varargin{1})
            error(['The argument must be 2-d vector: [subn, nperm]. ' ...
                   '"subn" is the dimension of sub-matrix; '...
                   '"nperm" is the number of permutation.']);
        end
        subn = varargin{1}(1);
        nperm = varargin{1}(2);
    end
    if strcmp(method, 'nicolo-overlap')
        nperm = varargin{1};
    end
    if strcmp(method, 'alex')
        if ~isscalar(varargin{1})
            error(['The last argument is the stability margin, ' ...
                   'which is a positive real number close to zero.']);
        end
        stabl_margin = varargin{1};
    end
    if sum(strcmp(method, {'sprandn', 'sprand'}))
        if ~isscalar(varargin{1})
            error(['The last argument is the sub-matrix dimension, ' ...
                   'which must be an integer.']);
        end
        subn = varargin{1};
    end
elseif nargin >= 5
    error('Too many arguments');
end


% main body
switch method
  case {'rss', 'nicolo', 'sprand', 'sprandn'}
    num_blk = floor(n/subn);
    last_blk_dim = n - subn*num_blk;
    A = [];
  case {'alex', 'nicolo-overlap'}
    ;
  otherwise
    error('Invalid argument for "method"!');
end

switch method
%   case 'rss'
%     for iblk = 1:(num_blk+1)
%
%         if iblk <= num_blk
%             sys = rss(subn);
%         else
%             sys = rss(last_blk_dim);
%         end
%         A = blkdiag(A, sys.A);
%     end

  case 'nicolo'
    for iblk = 1:(num_blk+1)
        while 1
            seed_randn = rng;
            if iblk <= num_blk
                Ablk = full(sprandn(subn,subn, density));
            else
                Ablk = full(sprandn(last_blk_dim, last_blk_dim, density));
            end
            % testing if A is stable
            eig_Ablk = eig(Ablk);
            cond_eigA = sum(real(eig_Ablk) >= 0);
            if ~cond_eigA
                break
            end
        end
        A = blkdiag(A, Ablk);
    end
    % generate random permutation matrix
    for iperm = 1:nperm
        P = eye(n);
        idx = randperm(n);
        P = P(idx,:);
        A = P*A*inv(P);
    end

  case 'nicolo-overlap'
    while 1
%         cornerVal = abs(randn(1));
        cornerVal = .53;
        Abase = -1*eye(n) + [[zeros(1,n-1) cornerVal]; [eye(n-1) zeros(n-1,1)]];
        Alayer = sprandn(n,n, density);
        A = Abase + Alayer;

        eig_A = eig(A);
        cond_eigAa = sum(real(eig_A) >= 0);
        if ~cond_eigAa
            break
        end
    end
    % generate random permutation matrix
    for iperm = 1:nperm
        P = eye(n);
        idx = randperm(n);
        P = P(idx,:);
        A = P*A*inv(P);
    end

  case 'alex'
    A = sprandn(n,n, density);
    eigv_A = eig(full(A));
    % find the largest unstable eigenvalue
    idx_unstbl = find(eigv_A >= 0);
    if ~isempty(idx_unstbl)
        max_eigv = max(eigv_A(idx_unstbl));
        min_eigv = min(eigv_A);
        % shift A such that all eigenvalues are negative
        %stabl_margin = 1e-2;  % otherwise, marginal stalbe
        A = A - (max_eigv + stabl_margin)*eye(n);
        % scale A to avoid too large negative eigv
        ratio = abs(real(min_eigv))/(real(max_eigv) - real(min_eigv));
        if ~isnan(ratio) && ratio
            A = 1/ratio * A;
        end
    end

  case 'sprandn'  % only block-diagonal matrices
    for iblk = 1:(num_blk+1)
        while 1
            seed_randn = rng;
            if iblk <= num_blk
                Ablk = full(sprandn(subn,subn, density));
            else
                Ablk = full(sprandn(last_blk_dim, last_blk_dim, density));
            end
            % testing if A is stable
            eig_Ablk = eig(Ablk);
            cond_eigA = sum(real(eig_Ablk) >= 0);
            if ~cond_eigA
                break
            end
        end
        A = blkdiag(A, Ablk);
    end
  case 'sprand' % only block-diagonal matrices
    for iblk = 1:(num_blk+1)
        while 1
            seed_randn = rng;
            if iblk <= num_blk
                Ablk = full(sprand(subn,subn, density));
            else
                Ablk = full(sprand(last_blk_dim, last_blk_dim, density));
            end
            % testing if A is stable
            eig_Ablk = eig(Ablk);
            cond_eigA = sum(real(eig_Ablk) >= 0);
            if ~cond_eigA
                break
            end
        end
        A = blkdiag(A, Ablk);
    end
end




% --------------------------------------------------------
% figure plot
plot_flag = 1;
if plot_flag
    fig_hl1 = figure('visible', 'off');
    set(fig_hl1, 'Units', 'inches', 'Position', [6.3472 2.8194 10.7917 9.4722]);
    upval = max(max(abs(A)));
    clims = [-upval, upval];
    imagesc(real(A), clims)
    ax = gca;
    ax.FontSize = 14;
    title('$A$ (ground truth)','Interpreter','latex')

    c_hl = colorbar;
    ymap = redblue(160);
    colormap(ymap);
    set(c_hl, 'Units', 'normalized', 'Position', [0.9323 0.1111 0.0192 0.8138]);

    fig_hl2 = figure('visible', 'off');
    set(fig_hl2, 'Units', 'inches', 'Position', [6.3472 6.5139 7.0556 5.7778]);
    pgt = plot(digraph(A'>= 1e-8));
    pgt.NodeColor = 'red';
    title('Ground Truth of the Dynamic Network');
end
