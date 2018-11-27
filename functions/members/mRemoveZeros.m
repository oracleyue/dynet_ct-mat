function vec = mRemoveZeros(rawVec, xLength)

diff = length(rawVec) - xLength;
vec = rawVec;
if diff
    vec(1:diff) = [];
end
