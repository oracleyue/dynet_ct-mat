function intVal = mIntVecExpAs(A, h)
% MINTVECEXPAS - to calculate the integration of expm(sA) over [0,h],
% using Trapezoidal numerical integration.

% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 19 Jun 2017


% Set precision
AbsTol = 1e-12;
RelTol = 1e-8;
numGrid = 1e3;

% Dimensions
n = size(A, 1);  % dim of A
vn = n^2;  % dim of vecA
AbsErr = 1;  % any value larger than AbsTol

while AbsErr > AbsTol || RelErr > RelTol
    % split [0,h] with numSteps and integrate
    funValList = zeros(numGrid, vn);
    k = 1;
    unitGridList = linspace(0,h, numGrid);
    for s = unitGridList
        funValList(k, :) = reshape(exp(s*A), [], 1)';
        k = k + 1;
    end
    intVal = trapz(unitGridList, funValList);

    % split [0,h] with numSteps+1 and integrate
    funValList = zeros(numGrid+1, vn);
    k = 1;
    unitGridList = linspace(0,h, numGrid+1);
    for s = unitGridList
        funValList(k, :) = reshape(exp(s*A), [], 1)';
        k = k + 1;
    end
    intValNew = trapz(unitGridList, funValList);

    % calculate error
    AbsErr = max(abs(intVal - intValNew));
    RelErr = max(abs((intVal - intValNew)./intVal));
    
    % update grids
    numGrid = numGrid*2;
end

intVal = intValNew';