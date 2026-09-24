function M = dc_time_match(tA, tB, opts)
%DC_TIME_MATCH  For each time in A, the nearest time in B — or nothing, when B was not looking.
%
%   M = dc_time_match(tA, tB)
%   M = dc_time_match(tA, tB, opts)
%
% The arithmetic under dual-colour tracking when the two colours were NOT acquired at the same rate.
% Colour B carrying half the frames of colour A is the common case: interleaved acquisition, or a
% strobed second channel. Half of A's localizations then have no simultaneous B frame at all, and
% what is done about that decides the answer.
%
% THE ONE RULE THIS ENFORCES: a localization with no partner in time is UNMEASURED, not far. It
% comes back with index 0, and .nUnmatched says how many there were. Treating it as a large distance
% would turn "B was not imaged then" into "B was not there then", which is the opposite conclusion
% and is invisible once it reaches a histogram.
%
% opts
%   .tol_s   how far in time a partner may be. Default: HALF of B's median spacing, so a frame
%            exactly between two B frames still matches and nothing matches across a gap. Pass 0 to
%            require exact simultaneity, or Inf to always take the nearest.
%   .rule    'nearest' (default) | 'previous' (the last B frame at or before tA) | 'next'.
%            'previous' is the honest rule when B is a structure that persists between its frames;
%            'nearest' is right when B is a molecule that moved.
%
% OUTPUT M
%   .idx        nA x 1, the index into tB, or 0 where nothing was within tol_s
%   .dt_s       nA x 1, the SIGNED gap (tB(idx) - tA), NaN where unmatched. Its spread is what says
%               how much of a measured distance is really the time gap.
%   .matched    nA x 1 logical
%   .nUnmatched how many of A had no partner
%   .tol_s .rule .medianSpacingB   what the numbers were produced under

if nargin < 3 || ~isstruct(opts), opts = struct(); end
tA = tA(:); tB = tB(:);
nA = numel(tA); nB = numel(tB);
rule = lower(getf(opts,'rule','nearest'));
assert(ismember(rule,{'nearest','previous','next'}), 'dc_time_match:rule', ...
    'rule must be nearest | previous | next, got ''%s''', rule);

spB = NaN;
if nB > 1
    d = diff(sort(tB(isfinite(tB))));
    d = d(d > 0);
    if ~isempty(d), spB = median(d); end
end
tol = getf(opts,'tol_s', []);
if isempty(tol)
    % Half of B's spacing: a frame exactly between two B frames still matches, and nothing matches
    % across a gap in B. With B at half of A's rate this matches EVERY A localization, which is the
    % intent — the time gap is then reported in .dt_s rather than hidden by dropping the frame.
    if isfinite(spB), tol = spB/2; else, tol = Inf; end
end

% THE BOUNDARY IS THE COMMON CASE, so it cannot be left to rounding. The default tolerance is
% exactly half of B's spacing, and interleaved acquisition puts every partner at exactly that
% distance — 0.14 - 0.13 evaluates to 0.010000000000000009, which a bare > rejects. Half the frames
% then come back unmatched for no reason a user could ever see. A relative slack of 1e-9 is far
% below any real timing and far above the arithmetic.
tolEff = tol * (1 + 1e-9) + 1e-12;

M = struct('idx',zeros(nA,1), 'dt_s',nan(nA,1), 'matched',false(nA,1), ...
           'nUnmatched',nA, 'tol_s',tol, 'rule',rule, 'medianSpacingB',spB);
if nA == 0 || nB == 0, return; end

[ts, ord] = sort(tB);                      % B need not arrive sorted, and a track rarely is
okB = isfinite(ts);
ts = ts(okB); ord = ord(okB);
if isempty(ts), return; end

for i = 1:nA
    t = tA(i);
    if ~isfinite(t), continue; end
    j = binSearch(ts, t);                  % last index with ts(j) <= t, 0 when t is before all
    switch rule
        case 'previous'
            cand = j;
        case 'next'
            cand = j + 1; if cand > numel(ts), cand = 0; end
        otherwise
            cand = 0; best = Inf;
            for c = [j, j+1]
                if c >= 1 && c <= numel(ts)
                    g = abs(ts(c) - t);
                    if g < best, best = g; cand = c; end
                end
            end
    end
    if cand < 1 || cand > numel(ts), continue; end
    gap = ts(cand) - t;
    if abs(gap) > tolEff, continue; end
    M.idx(i) = ord(cand);                  % back to the caller's own ordering
    M.dt_s(i) = gap;
    M.matched(i) = true;
end
M.nUnmatched = nnz(~M.matched);
end

% =================================================================================================
function j = binSearch(ts, t)
% Last index with ts(j) <= t; 0 when t is before every element. ts is sorted ascending.
lo = 1; hi = numel(ts); j = 0;
while lo <= hi
    mid = floor((lo+hi)/2);
    if ts(mid) <= t, j = mid; lo = mid + 1; else, hi = mid - 1; end
end
end

function v = getf(s,f,d), if isstruct(s)&&isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
