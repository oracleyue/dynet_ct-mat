function vec = mAppendOnes(rawVec, xLength, val)

if nargin == 2, val = 1; end

diff = xLength - length(rawVec);
if diff
    vec = [rawVec; ones(diff,1)*val];
else
    vec = rawVec;
end
