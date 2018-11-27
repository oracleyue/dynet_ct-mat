function y = mObjFuncRPEM(A, B, h, PhiA, PhiB)
% MOJBFUNCRPEM - This function is to calculate the objective value of the
% optimizaiton problem in regularized PEM.

% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 19 Jun 2017


funExpAs = @(s) expm(s*A);
intExpAh = integral(funExpAs, 0, h);

y = norm(y - PhiA*reshape(expm(h*A),[],1) ...
           - PhiB * reshape(intExpAh, [], 1) )^2 ...
    + lambda*norm(reshape(A,[],1),1);
