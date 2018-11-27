function [resName, perfMetric] = spid_sysid_noinput(dataName, dataPath, resPath, spidOptions)
% SPID_SYSID_NOINPUT perform network reconstruction with
% low-sampling-frequency time series.
%
% Inputs:
%   - dataName: name/identifier of data file
%   - dataPath: path to keep data sets
%   - resPath:  path to save results
%   - spidOptions: options of solvers, see comments below.
%
% Outputs:
%   - resName:    full path of .mat file of results
%   - perfMetric: 7x1 vector, [TP, FP, FN, TN, lambda, TypeIError, TypeIIError]

% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 30 Jan 2018



debug_flag = 0;     % flag to print debug info

%% Parse input arguments
spid_optvar_precision = spidOptions.optvar_precision; % |A_k+1 - A_k|_2 <= ...
spid_optval_precision = spidOptions.optval_precision; % |fA_k+1 - fA_k|_2 <= ...
spid_ls_precision     = spidOptions.ls_precision;     % ls - line search
spid_ls_alpha         = spidOptions.ls_alpha;         % default: <1/2
spid_iteration_loops  = spidOptions.iteration_loops;
spid_init_gamma       = spidOptions.init_gamma;       % set initial point
spid_sparsity_method  = spidOptions.sparsity_method;

lambda = spidOptions.lambda;
zero   = spidOptions.zero;


%% Load Dataset
load([dataPath dataName '.mat']);
n = sys_dim;   % dim of A
m = ts_dim;    % length of time series
h = Ts;        % sampling period h


%% Problem Setup
y = reshape(X1, [], 1);     % vec(X1)
Phi = kron(X2', eye(n));    % (X_2^T \otims I_n)
Ad = X1/X2;


%% Variables to Check Convergence
residual_list = [];
gap_gt_list = [];
% vecAk_intermdiate_list = [];
vecA_gt = reshape(A_gt, [], 1);


%% Optimization Solver
A_logm = logm(Ad)/h;     % estimation via logm of Ad
% A_k = spid_init_gamma*(A_logm);     % initial point A_0
A_k = zeros(size(Ad));

f = @(A) norm(y - Phi*reshape(expm(h*A),[],1))^2 ...
         + lambda*norm(reshape(A,[],1),1);
f_A = f(A_k);  % value of objective function at A0

for k = 0:1:spid_iteration_loops
    [cond_fre, K_A] = expm_cond(h*A_k);
    r_A = y - Phi * reshape(expm(h*A_k), [], 1);
    J_A = - h * Phi * K_A;
    naphi_A = 2*J_A'*r_A;    % \nabla \phi(A)
    vec_Ak = reshape(A_k, [], 1);
    gBar_A = naphi_A + lambda*sign(vec_Ak);
    W_A = eye(n^2) - diag(abs(sign(vec_Ak)));

    % Save for checking convergence
    residual_list = [residual_list norm(r_A)^2];
    gap_gt_list = [gap_gt_list norm(vec_Ak - vecA_gt, 1)];

    % Optimize the local convex problem (Gauss-Newton-like Method)
    switch spid_sparsity_method
      case 'lasso'
        % using LASSO
        cvx_begin quiet
            variable p_k(n^2)
            minimize( norm(r_A + J_A*p_k) + lambda*norm(p_k + vec_Ak,1) )
            subject to
                gBar_A'*p_k + lambda*norm(W_A*p_k, 1) <= 0
        cvx_end
      case 'sbl'
        % using SBL
        learn_Lambda = 0;
        y_sbl = r_A - J_A * vec_Ak;
        Phi_sbl = -J_A;
        [p_k, gamma_est, gamma_used, count] ...
            = MSBL(Phi_sbl, y_sbl, lambda, learn_Lambda);
    end
    E_k = reshape(p_k, [], n);

    % Line Search: backtracking line search
    s_k = 1;  % initialize s_k
    alpha = spid_ls_alpha;
    beta = .5;
    fPrime_A = gBar_A'*p_k + lambda*norm(W_A*p_k, 1);

    flagLineSearch = 0;   % line search fails when flag = 1
    if debug_flag
        fprintf('\n------ Loop %i ------\n', k);
        fprintf('f(A_%i): %d\n', k,f_A);
        fprintf('-- Line Search: Start\n')
    end
    while 1
        f_A_next = f(A_k + s_k*E_k);
        % f_A = f(A_k);    % the saved f_A_next in the previous round

        if debug_flag
            fprintf('f(A_%i + s_k*p_%i): %d\n', k,k,f_A_next);
        end

        if f_A_next - f_A <= alpha*s_k*fPrime_A
            break
        elseif s_k*beta > spid_ls_precision
            s_k = s_k*beta;
            continue
        end

        flagLineSearch = 1;    % line search method failed
        warning(sprintf(['Linear search failed! (the lower bound for ' ...
                         's_k is %.2e)'], spid_ls_precision));
        break
    end
    if debug_flag
        fprintf('-- Line Search: End\n')
        fprintf('f(A_%i): %d\n', k+1, f_A_next);
    end

    % Stopping Criteria
    diffObjectFunction = abs(f_A - f_A_next);
    diffOptimalVariable = norm(s_k * p_k, 2);

    if debug_flag
        fprintf('|f(A_%i) - f(A_%i)|: %d\n', k+1, k, diffObjectFunction);
        fprintf('|A_%i - A_%i|: %d\n', k+1, k, diffOptimalVariable);
    end

%     if diffOptimalVariable < spid_optvar_precision
    if diffOptimalVariable < spid_optvar_precision ...
        ||  diffObjectFunction < spid_optval_precision
        solver_status = 'success: standard stopping criteria';
        break
    end

    % Update
    if flagLineSearch
        solver_status = 'unsuccess: stop at line search conditions';
        warning(sprintf(['solver stops at %d-th iteration due to the failure ' ...
                         'of linear search.'], k))
        break
    else
        A_k = A_k + s_k*E_k;     % update
        f_A = f_A_next;
    end

    % Run out of loops
    if k == spid_iteration_loops
        solver_status = 'underdetermined success: run out loops';
        warning(['Solver Status:    ' solver_status]);
    end
end
Ak_d = expm(h*A_k);

% Print breif summary of results
if debug_flag
    fprintf('\nStatus of Solver\n----------------\n%s\n\n', solver_status);
end


%% Construct Dynamic Boolean Network
netGT = double(abs(A_gt) > zero);
netRes = double(abs(A_k) > zero);
TP = sum(sum(netGT & netRes));
FP = sum(sum(~netGT & netRes));
FN = sum(sum(netGT & ~netRes));
TN = sum(sum(~netGT & ~netRes));
Prec = TP / (TP + FP);       % 1 - Prec = numWrongLinks / numTotalLinksRes
TPR  = TP / (TP + FN);       % 1 - TPR =  numMissedLinks / numTotalLinksGT
TypeIError = (1-Prec)*100;
TypeIIError = (1-TPR)*100;

% fprintf('Type I  Error: %f%% \n', TypeIError);
% fprintf('Type II Error: %f%% \n\n', TypeIIError);


%% Export results
perfMetric = [TP FP FN TN lambda TypeIError TypeIIError];
resName = [resPath dataName '_results.mat'];
save(resName, 'A_gt', 'A_k', 'A_logm', 'Ad_gt', 'Ak_d', 'Ad', ...
     'residual_list', 'X1','X2','h','gap_gt_list');
