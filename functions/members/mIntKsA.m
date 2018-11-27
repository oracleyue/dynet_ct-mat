function intVal = mIntKsA(A, h)
% MINTSKSA - to calculate the integration of sK(sA) over [0,h],
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
vn = n^4;  % dim of vecA
AbsErr = 1;  % any value larger than AbsTol

while AbsErr > AbsTol || RelErr > RelTol
    % -- debug
    numGrid

    tic
    % split [0,h] with numSteps and integrate
    funValList = zeros(numGrid, vn);    
    unitGridList = linspace(0,h, numGrid);
    parfor k = 1:1:numGrid
        s = unitGridList(k);
        funValList(k,:) = calcKsA(A, s)';
    end
    intVal = trapz(unitGridList, funValList);

    % split [0,h] with numSteps+1 and integrate
    funValList = zeros(numGrid+1, vn);
    unitGridList = linspace(0,h, numGrid+1);
    parfor k = 1:1:numGrid+1
        s = unitGridList(k);
        funValList(k,:) = calcKsA(A, s)';
    end
    intValNew = trapz(unitGridList, funValList);

    % calculate error
    AbsErr = max(abs(intVal - intValNew));
    RelErr = max(abs((intVal - intValNew)./intVal));

    % update grids
    numGrid = numGrid*2;

    AbsErr, RelErr
    toc
end

intVal = reshape(intValNew', n^2, n^2);


end



% Local Functions

function y = calcKsA(A, s)
% To calculate the value of $sK(sA)$

[~, KsA] = expm_cond(s*A);
y = s*reshape(KsA, [], 1);

end