function R = dc_pbsa_steps(y, varargin)
%DC_PBSA_STEPS  Find discrete photobleaching steps in an intensity trace (Kalafut-Visscher + SIC).
%
%   R = dc_pbsa_steps(y)
%   R = dc_pbsa_steps(y, 'MinStep', h, 'MaxSteps', n, 'Verbose', tf)
%
% Counting the discrete drops in a photobleaching trace gives the number of fluorophores that were
% emitting, and hence the stoichiometry of whatever they label. This is the first stage of the
% method in Hummert et al., Mol Biol Cell 32(21):ar35 (2021) — a reimplementation from the published
% algorithm, not a port of their quickPBSA package (which is GPLv3; nothing of it is used here).
%
% THE ALGORITHM (Kalafut & Visscher, Comput Phys Commun 179:716, 2008). Fit the trace as a
% piecewise-constant function. Start with one segment. Repeatedly try inserting one more step at
% every remaining position, keep whichever insertion most reduces the Schwarz Information Criterion,
% and stop when no insertion reduces it further. SIC punishes each extra step by log(N), so the fit
% stops adding steps exactly when the variance they explain stops paying for their cost — which is
% what makes the method model-independent: nothing here assumes how many fluorophores there are.
%
%   SIC = (k + 2)*log(N) + N*log(sigma^2)
%
% with k steps, N points, and sigma^2 the residual variance of the piecewise-constant fit. The two
% extra parameters are the first level and the variance itself.
%
% WHY THIS IS FAST. The naive search is O(N^2) per step because each candidate needs the residual of
% the whole trace. Cumulative sums of y and y^2 make any segment's sum-of-squared-error O(1), and
% inserting a step only ever splits ONE segment — so each round rescans only the segment that
% changed, not all of them. On a 400-point trace this is the difference between milliseconds and
% seconds, which matters because the viewer re-runs it on every track you click.
%
% INPUT
%   y : intensity trace, a vector. NaNs are dropped (and reported in R.nDropped) — a gap-filled
%       track legitimately has them, and a NaN would poison every segment mean containing it.
%
% OPTIONS
%   'MinStep'  ('auto') reject fitted steps smaller than this in |height| after the search.
%                      'auto' (the default) uses 3x a robust noise estimate; pass a number to set it
%                      by eye at about half a single-fluorophore step; pass 0 to disable it and get
%                      the unmodified SIC result.
%
%                      WHY A THRESHOLD IS THE DEFAULT AND NOT ZERO. SIC charges log(N) per step, so
%                      on an N=120 trace a spurious step only has to cut the residual by 3.9% to pay
%                      for itself — and the search picks the best of ~N candidate positions, so pure
%                      noise clears that bar routinely. Measured here: a clean 1-step trace at noise
%                      sd 8 reliably gains a second "step" of height ~8, i.e. exactly the noise. The
%                      published method has the same property, which is why quickPBSA makes the
%                      threshold a user parameter set by visual inspection. Defaulting to 0 would
%                      hand back a dimer where the data shows a monomer, so the default estimates
%                      the noise instead of ignoring it.
%   'MaxSteps' (200)   hard cap on the search, so a pathological trace cannot spin.
%   'Verbose'  (false) print the SIC at each accepted step.
%
% OUTPUT  R, a struct:
%   k         number of steps kept
%   idx       1xk step positions, as the index of the LAST point before each step
%   levels    1x(k+1) fitted level of each segment
%   heights   1xk signed step heights (negative = a drop, which is what bleaching gives)
%   fit       piecewise-constant fit, same length as the input (NaN where the input was NaN)
%   sic       SIC after each accepted step, starting with the no-step value
%   sigma     residual standard deviation of the final fit
%   nDropped  how many NaNs were removed before fitting
%
% A photobleaching trace should give all-negative heights. A mixture of signs means the trace is not
% a clean bleach — blinking, a tracking error, or two particles crossing — and R.heights is the
% place to see that rather than a count that silently averages it away.
p = inputParser;
p.addParameter('MinStep',  'auto', @(v) (isnumeric(v) && isscalar(v) && v >= 0) || ...
                                        (ischar(v) || isstring(v)) && strcmpi(char(v),'auto'));
p.addParameter('MaxSteps', 200,   @(v) isnumeric(v) && isscalar(v) && v >= 1);
p.addParameter('Verbose',  false, @(v) islogical(v) || isnumeric(v));
p.parse(varargin{:});
minStep  = p.Results.MinStep;
maxSteps = p.Results.MaxSteps;
autoMin  = (ischar(minStep) || isstring(minStep)) && strcmpi(char(minStep),'auto');
verbose  = logical(p.Results.Verbose);

R = struct('k',0,'idx',[],'levels',[],'heights',[],'fit',[],'sic',[],'sigma',NaN,'nDropped',0, ...
           'minStepUsed',NaN,'noiseSd',NaN);
if nargin < 1 || isempty(y), return; end

yIn  = y(:);
good = isfinite(yIn);
R.nDropped = sum(~good);
x = yIn(good);
N = numel(x);
if N < 4                                    % below this a "step" is indistinguishable from noise
    R.fit = nan(size(yIn));
    if N >= 1, R.levels = mean(x); R.fit(good) = mean(x); R.sigma = std(x); end
    return
end

% ---- robust noise estimate, for the automatic step threshold ------------------------------------
% From the first differences rather than the trace itself, so the real steps do not inflate it: a
% difference of independent noise has variance 2*sigma^2, hence the 1/sqrt(2). MAD (not sd) because
% the handful of genuine bleaching steps ARE outliers here and must not count toward the noise.
d = diff(x);
noiseSd = 1.4826 * median(abs(d - median(d))) / sqrt(2);
% A MAD of zero means a genuinely flat-between-steps trace (synthetic, or a very bright emitter):
% the noise really is ~0, so the threshold should be 0 and every step SIC finds is real. Falling
% back to std(x) here would be backwards — std includes the steps themselves, so a clean 2-step
% trace would set a threshold larger than its own steps and return nothing.
if ~isfinite(noiseSd) || noiseSd < 0, noiseSd = 0; end
R.noiseSd = noiseSd;
if autoMin, minStep = 3 * noiseSd; end
R.minStepUsed = minStep;

% ---- cumulative sums: any segment's SSE in O(1) -------------------------------------------------
cs  = [0; cumsum(x)];
cs2 = [0; cumsum(x.^2)];
    function s = sse(a, b)
        % SSE of x(a:b) about its own mean. The max(.,0) guards the catastrophic cancellation you
        % get from the sum-of-squares form when a segment is long and nearly constant — without it
        % this returns a small negative number and log(sigma^2) becomes complex.
        n = b - a + 1;
        s = max((cs2(b+1) - cs2(a)) - (cs(b+1) - cs(a))^2 / n, 0);
    end
    function [bp, gain] = bestSplit(a, b)
        % Best single split of x(a:b), by SSE reduction. Returns bp = last index of the left part.
        % Vectorised over every interior cut at once — this is the inner loop of the whole method.
        bp = 0; gain = 0;
        if b - a + 1 < 2, return; end
        m  = (a:b-1)';                                   % candidate left-end indices
        nL = m - a + 1;              nR = b - m;
        sL = cs(m+1) - cs(a);        sR = cs(b+1) - cs(m+1);
        qL = cs2(m+1) - cs2(a);      qR = cs2(b+1) - cs2(m+1);
        tot = max(qL + qR - (sL+sR).^2 ./ (nL+nR), 0);   % SSE with no split (constant across m)
        spl = max(qL - sL.^2./nL, 0) + max(qR - sR.^2./nR, 0);
        [g, i] = max(tot - spl);
        if isfinite(g) && g > 0, bp = m(i); gain = g; end
    end

% ---- greedy search, one step at a time ----------------------------------------------------------
bounds = [1 N];                       % segment start/end pairs, one row each
totSSE = sse(1, N);
sicOf  = @(k, s) (k + 2)*log(N) + N*log(max(s, realmin)/N);
curSIC = sicOf(0, totSSE);
sicTrace = curSIC;
steps  = [];

% Cache each segment's best split so a round only recomputes the segment that actually changed.
nSeg = 1;
bp   = zeros(1,1); gn = zeros(1,1);
[bp(1), gn(1)] = bestSplit(1, N);

for it = 1:maxSteps
    [bestGain, seg] = max(gn(1:nSeg));
    if ~(bestGain > 0), break; end
    cut = bp(seg);
    newSSE = totSSE - bestGain;
    newSIC = sicOf(numel(steps) + 1, newSSE);
    if ~(newSIC < curSIC), break; end             % SIC stopped paying for the extra step

    a = bounds(seg,1); b = bounds(seg,2);
    bounds = [bounds(1:seg-1,:); a cut; cut+1 b; bounds(seg+1:end,:)];
    steps  = sort([steps cut]);
    totSSE = newSSE; curSIC = newSIC; sicTrace(end+1) = newSIC; %#ok<AGROW>

    % only the split segment's two halves need re-evaluating
    [l1, g1] = bestSplit(a, cut);
    [l2, g2] = bestSplit(cut+1, b);
    bp = [bp(1:seg-1) l1 l2 bp(seg+1:nSeg)];
    gn = [gn(1:seg-1) g1 g2 gn(seg+1:nSeg)];
    nSeg = nSeg + 1;
    if verbose, fprintf('  step %d at %d -> SIC %.2f\n', numel(steps), cut, newSIC); end
end

% ---- levels, then the height filter -------------------------------------------------------------
[lv, hh] = levelsOf(steps);
if minStep > 0 && ~isempty(steps)
    % Drop steps below the threshold SMALLEST-FIRST, recomputing levels each time: merging two
    % segments changes the neighbouring heights, so a single vectorised mask would drop steps that
    % are only small because of a step next to them that is itself about to go.
    while ~isempty(steps)
        [mn, w] = min(abs(hh));
        if mn >= minStep, break; end
        steps(w) = [];
        [lv, hh] = levelsOf(steps);
    end
end

R.k       = numel(steps);
R.idx     = steps;
R.levels  = lv;
R.heights = hh;
R.sic     = sicTrace;
fitG = zeros(N,1);
e = [steps N];  s = [1 steps+1];
for j = 1:numel(lv), fitG(s(j):e(j)) = lv(j); end
R.fit = nan(size(yIn));
R.fit(good) = fitG;
R.sigma = sqrt(max(sum((x - fitG).^2), 0) / N);

    function [lv, hh] = levelsOf(st)
        e_ = [st N]; s_ = [1 st+1];
        lv = zeros(1, numel(e_));
        for j_ = 1:numel(e_)
            lv(j_) = (cs(e_(j_)+1) - cs(s_(j_))) / (e_(j_) - s_(j_) + 1);
        end
        hh = diff(lv);
    end
end
