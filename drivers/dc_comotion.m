function R = dc_comotion(S, opts)
%DC_COMOTION  Do nearby trajectories move together? Step dot products, binned by separation.
%
%   R = dc_comotion(S)                                   S from dc_steps (drift-corrected)
%   R = dc_comotion(S, struct('rMaxUm', 3, 'nMin', 20))
%
% For every pair of tracks that have a step starting at the SAME timepoint, this takes the two step
% vectors u and v and records
%
%       dot  = u . v                 um^2   — co-movement with its amplitude
%       cos  = u . v / (|u| |v|)            — direction alone, blind to how fast either was going
%       r    = |A(t) - B(t)|         um     — how far apart they were when they took the step
%
% and returns them two ways: every step binned by r (the curve), and every pair pooled over its own
% shared steps (which pairs).
%
% COS OR DOT. Use cos as the primary. A dot product is dominated by whichever pair happened to take
% the biggest steps, so a single fast excursion by two unrelated molecules outweighs a hundred quiet
% steps of a genuinely coupled pair; it is also in um^2, so it cannot be compared between datasets
% with different diffusion coefficients. cos answers "did they go the same way", which is the
% question, and its null does not depend on the step-size distribution. dot is kept because
% co-movement with amplitude is the stronger claim when you can make it.
%
% THE THIRD STATISTIC. corrNorm = <u.v> / sqrt(<|u|^2><|v|^2>) is a Pearson correlation of the step
% vectors. Between cos and dot: it weights by amplitude like dot, but is dimensionless like cos. For
% a pair whose steps vary a lot in size it is less noisy than mean-cos.
%
% THE BINNED CURVE IS DILUTED BY CHANCE ENCOUNTERS, AND THE PAIRS ARE NOT. A separation bin holds
% every pair that happened to be that far apart, so a genuinely coupled pair sits in it alongside all
% the unrelated tracks that merely wandered close. On a synthetic field of 42 tracks in 8 x 8 um with
% 14 coupled pairs, the coupled pairs read cos = 0.50 individually while the r <= 0.5 um bin reads
% 0.13 — roughly three quarters of the near-field step pairs were incidental. So .bins answers "is
% there co-motion at short range in this movie", and its height is NOT the coupling strength; .pairs
% answers "which pairs, and how strongly", and that is where the 0.50 is visible. Reading a coupling
% constant off the curve underestimates it by whatever the local density happens to be.
%
% (For 2D Gaussian steps sharing a fraction rho of their variance, E[cos] ~ (pi/4)*rho for small rho:
% a per-component correlation of 0.6 gives a mean cosine near 0.50, not 0.6.)
%
% WHAT THE NUMBER IS DILUTED BY. Localization error sigma adds an independent random vector to each
% step, of variance 2*sigma^2 per axis. It is independent between tracks, so it cannot CREATE
% correlation — it only shrinks the measured value toward zero, by roughly
% var_motion / (var_motion + 2 sigma^2) for each track. A measured 0.2 with noisy localization can be
% a true 0.4. This matters for interpreting the size of an effect and not at all for whether one
% exists.
%
% PAIRS OF THE SAME COLOUR ARE NOT THE SAME EVIDENCE. Two tracks of one colour that come within about
% a PSF can have their identities swapped by the linker, which manufactures correlation at exactly the
% separations of interest. Two tracks of DIFFERENT colours cannot: neither linker ever saw the other
% channel. R splits the pairs into 'cross' and 'same' for that reason; believe the cross-colour curve
% at small r, and treat the same-colour curve there as an upper bound.
%
% OPTIONS
%   .rMaxUm   3      ignore pairs further apart than this at that step (keeps the pair count sane)
%   .nMin     20     a pair needs this many shared steps to get a row in R.pairs
%   .edgesUm  []     separation bin edges; default 0:0.1:1, then 1.25:0.25:2, then 2.5, 3
%   .maxPairs 5e6    guard: stop and say so rather than filling memory
%
% OUTPUT R
%   .steps   table, one row per simultaneous step pair: trackA trackB chA chB tp r dot cos
%   .bins    table, one row per separation bin: rMid n meanCos seCos meanDot corrNorm, per class
%   .pairs   table, one row per pair with >= nMin shared steps: trackA trackB class n rMedian
%            meanCos seCos z meanDot corrNorm
%   .params  what was used
%
% Judge a number in .bins or .pairs against dc_comotion_null, never against zero: zero is the null
% only if there is no drift, no flow and no shared confinement, and the far field measures whether
% that is true here.

if nargin < 2 || ~isstruct(opts), opts = struct(); end
rMax  = getf(opts,'rMaxUm', 3);
nMin  = getf(opts,'nMin', 20);
maxP  = getf(opts,'maxPairs', 5e6);
edges = getf(opts,'edgesUm', [0:0.1:1, 1.25:0.25:2, 2.5, 3]);
edges = unique([edges(:)' rMax]);
edges = edges(edges <= rMax);

assert(~isempty(S), 'dc_comotion:noSteps', 'no steps to pair');
assert(numel(unique(S.span)) == 1, 'dc_comotion:mixedSpan', ...
    ['the step table mixes spans (%s). A step over two timepoints is not simultaneous with one over ' ...
     'a single timepoint and is on average longer; pair only steps of equal span (dc_steps defaults ' ...
     'to span 1).'], mat2str(unique(S.span)'));

% ---- pair up every pair of steps that start at the same timepoint ------------------------------
[g, tps] = findgroups(S.tp0);
nPer = accumarray(g, 1);
est  = sum(nPer .* (nPer-1) / 2);
assert(est <= maxP, 'dc_comotion:tooManyPairs', ...
    ['%.3g simultaneous pairs, over the %.3g cap. Raise maxPairs, shorten the movie, or raise the ' ...
     'detection threshold — this grows as the square of the number of tracks alive at once.'], est, maxP);

A = zeros(0,1); B = zeros(0,1); TP = zeros(0,1); RR = zeros(0,1); DT = zeros(0,1); CS = zeros(0,1);
uA = zeros(0,2); uB = zeros(0,2);
for k = 1:numel(tps)
    idx = find(g == k);
    m = numel(idx);
    if m < 2, continue; end
    [i1, i2] = find(triu(true(m), 1));
    a = idx(i1); b = idx(i2);
    dx = S.x(a) - S.x(b);  dy = S.y(a) - S.y(b);
    r  = hypot(dx, dy);
    in = r <= rMax;
    if ~any(in), continue; end
    a = a(in); b = b(in); r = r(in);
    ua = [S.ux(a) S.uy(a)];  ub = [S.ux(b) S.uy(b)];
    d  = sum(ua .* ub, 2);
    la = hypot(ua(:,1),ua(:,2)); lb = hypot(ub(:,1),ub(:,2));
    c  = d ./ max(la .* lb, eps);
    c(la <= 0 | lb <= 0) = NaN;            % a zero-length step has no direction
    A = [A; S.trackId(a)]; B = [B; S.trackId(b)]; %#ok<AGROW>
    TP = [TP; repmat(tps(k), numel(a), 1)]; RR = [RR; r]; DT = [DT; d]; CS = [CS; c]; %#ok<AGROW>
    uA = [uA; ua]; uB = [uB; ub]; %#ok<AGROW>
end

[uid, iu] = unique(S.trackId);                 % one row per track, not one per step
chOf = containers.Map(num2cell(uid), cellstr(string(S.ch(iu))));
chA = string(values(chOf, num2cell(A)));      % values() keeps the input's shape: already Nx1
chB = string(values(chOf, num2cell(B)));
cls = repmat("same", numel(A), 1);
cls(chA ~= chB) = "cross";

% The step VECTORS are kept, not just the products: the time-shifted null in dc_comotion_null needs
% to pair A's step with B's from another moment, which cannot be recovered from a cosine.
R.steps = table(A, B, chA, chB, cls, TP, RR, DT, CS, uA(:,1), uA(:,2), uB(:,1), uB(:,2), ...
    'VariableNames', {'trackA','trackB','chA','chB','class','tp','r','dot','cos', ...
                      'uax','uay','ubx','uby'});
R.params = struct('rMaxUm',rMax, 'nMin',nMin, 'edgesUm',edges, 'nSteps',height(S), ...
                  'nTracks',numel(unique(S.trackId)), 'nStepPairs',numel(A));

% ---- the curve: every step, binned by separation ------------------------------------------------
R.bins = binUp(R.steps, uA, uB, edges, "all");
for c = ["cross","same"]
    m = R.steps.class == c;
    if any(m), R.bins = [R.bins; binUp(R.steps(m,:), uA(m,:), uB(m,:), edges, c)]; end
end

% ---- the pairs ----------------------------------------------------------------------------------
assert(max([A;B]) < 1e7, 'dc_comotion:trackIdRange', 'track ids must be below 1e7 to pack a pair key');
key = A * 1e7 + B;
[uk, ~, ic] = unique(key);
n   = accumarray(ic, 1);
sel = n >= nMin;
if ~any(sel)
    R.pairs = emptyPairs(); return
end
mC  = accumarray(ic, R.steps.cos, [], @(v) mean(v,'omitnan'));
sC  = accumarray(ic, R.steps.cos, [], @(v) std(v,'omitnan'));
mD  = accumarray(ic, R.steps.dot, [], @(v) mean(v,'omitnan'));
mR  = accumarray(ic, R.steps.r,   [], @median);
sAA = accumarray(ic, sum(uA.^2,2), [], @mean);
sBB = accumarray(ic, sum(uB.^2,2), [], @mean);
cn  = mD ./ sqrt(max(sAA .* sBB, eps));
first = accumarray(ic, (1:numel(ic))', [], @min);

se = sC ./ sqrt(max(n,1));
R.pairs = table(floor(uk/1e7), mod(uk,1e7), R.steps.class(first), n, mR, mC, se, mC./max(se,eps), mD, cn, ...
    'VariableNames', {'trackA','trackB','class','n','rMedian','meanCos','seCos','z','meanDot','corrNorm'});
R.pairs = sortrows(R.pairs(sel,:), 'z', 'descend');
end

% =================================================================================================
function T = binUp(St, uA, uB, edges, lbl)
b = discretize(St.r, edges);
ok = isfinite(b);
b = b(ok); St = St(ok,:); uA = uA(ok,:); uB = uB(ok,:);
nb = numel(edges) - 1;
n = accumarray(b, 1, [nb 1]);
mC = accumarray(b, St.cos, [nb 1], @(v) mean(v,'omitnan'), NaN);
sC = accumarray(b, St.cos, [nb 1], @(v) std(v,'omitnan'),  NaN);
mD = accumarray(b, St.dot, [nb 1], @(v) mean(v,'omitnan'), NaN);
aa = accumarray(b, sum(uA.^2,2), [nb 1], @mean, NaN);
bb = accumarray(b, sum(uB.^2,2), [nb 1], @mean, NaN);
T = table(repmat(lbl, nb, 1), (edges(1:end-1)+edges(2:end))'/2, edges(1:end-1)', edges(2:end)', ...
    n, mC, sC./sqrt(max(n,1)), mD, mD./sqrt(max(aa.*bb,eps)), ...
    'VariableNames', {'class','rMid','rLo','rHi','n','meanCos','seCos','meanDot','corrNorm'});
T = T(n > 0, :);
end

function T = emptyPairs()
T = table(zeros(0,1), zeros(0,1), strings(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
          zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', {'trackA','trackB','class','n','rMedian','meanCos','seCos','z','meanDot','corrNorm'});
end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f), v = s.(f); else, v = d; end
end
