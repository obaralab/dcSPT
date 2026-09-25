function N = dc_comotion_null(R, opts)
%DC_COMOTION_NULL  How many steps before a co-motion measurement stops being noise.
%
%   N = dc_comotion_null(R)
%   N = dc_comotion_null(R, struct('rFarUm', 1.5, 'rNearUm', 0.5))
%
% R comes from dc_comotion. This builds the distribution a pair's mean cosine takes WHEN NOTHING IS
% GOING ON, as a function of how many steps went into it, and reads two things off it: the number of
% steps a pair needs before a given correlation is distinguishable from noise, and how much of the
% near-field signal survives that bar.
%
% TWO NULLS, BECAUSE THEY ANSWER DIFFERENT QUESTIONS.
%
%   FAR PAIRS (r > rFarUm). The same tracks, the same frames, the same drift, the same localization
%   error — only too far apart to be interacting. This is the honest baseline for "is the near field
%   above the background of this movie", and it is the one to quote. It does NOT centre on zero if
%   drift or a flow survived dc_drift, and that is the point: it measures whatever is left.
%
%   TIME-SHIFTED. The same pairs, but one track's steps are rotated in time by a large lag, so
%   simultaneity is destroyed while each track keeps its own step-size and direction distribution.
%
%   Read this one carefully, because it does NOT simply centre on zero. It removes only
%   TIME-VARYING common motion. A steady drift gives every step the same mean vector, so two steps
%   taken at different times still share it and the shifted null stays just as high as the far one —
%   measured on a fixture with a constant 21 nm/step drift, far = +0.127 and shifted = +0.128. What
%   the shifted null actually measures is the correlation explained by a persistent directional bias
%   in the steps, with no simultaneity required.
%
%   So the two nulls decompose the far field rather than one of them being "the clean one":
%       shifted            correlation from a persistent directional bias (a steady drift, a fixed flow)
%       far - shifted      correlation that genuinely needs the two to have moved AT THE SAME MOMENT
%                          (a wandering stage, a time-varying flow)
%   Reported as .driftShare. Near zero means any drift here is steady; large means the stage moved
%   around during the movie, which dc_drift is the fix for.
%
% WHY NOT THE TEXTBOOK NUMBER. For independent uniform directions in 2D, cos has mean 0 and variance
% 1/2, so the mean of n steps has SE = 1/sqrt(2n) and n = 2/c^2 steps would resolve a correlation c.
% That number is optimistic here, usually by a lot: consecutive steps of one track are not
% independent (confinement, drift, and the anticorrelation that localization error puts into
% neighbouring steps), so the effective count is smaller than n. This function measures the SE that
% the data actually has at each n and fits SE = k/sqrt(n); k > 1/sqrt(2) is the size of that penalty.
%
% OPTIONS
%   .rNearUm  0.5    pairs at or below this are "near" — the ones being tested
%   .rFarUm   1.5    pairs beyond this form the null
%   .nGrid    []     step counts to evaluate; default geometric from 5 to the largest pair
%   .nBoot    2000   resamples per n
%   .alpha    0.05   one-sided, since the hypothesis is positive correlation
%
% OUTPUT N
%   .curve    per n: seFar, seShift, thresh (the (1-alpha) quantile of the far null), and the power
%             to detect the near field's own median correlation
%   .k        the fitted SE = k/sqrt(n) for the far null, and kIdeal = 1/sqrt(2)
%   .nNeeded  steps needed to resolve a correlation of 0.1, 0.2, 0.3 and the observed near-field one
%   .far .near   the two pair populations, as used
%   .text     the answer in one paragraph

if nargin < 2 || ~isstruct(opts), opts = struct(); end
rNear = getf(opts,'rNearUm', 0.5);
rFar  = getf(opts,'rFarUm', 1.5);
nBoot = getf(opts,'nBoot', 2000);
alpha = getf(opts,'alpha', 0.05);

assert(~isempty(R.steps), 'dc_comotion_null:noSteps', 'nothing to build a null from');
St = R.steps(isfinite(R.steps.cos), :);

farS  = St(St.r >  rFar,  :);
nearS = St(St.r <= rNear, :);
assert(height(farS) >= 200, 'dc_comotion_null:thinFar', ...
    ['only %d step pairs beyond %.2f um. The null is built from these, so it needs to be the big ' ...
     'population; lower rFarUm or raise rMaxUm in dc_comotion.'], height(farS), rFar);

% ---- the far-field null, resampled at each n ----------------------------------------------------
nMaxPair = max([R.pairs.n; 50]);
nGrid = getf(opts,'nGrid', unique(round(logspace(log10(5), log10(min(nMaxPair, 2000)), 14))));

cFar   = farS.cos;
cShift = shiftedCos(St);                 % the time-rotated null, over all separations

rng(7);
seFar = nan(numel(nGrid),1); seShift = nan(numel(nGrid),1);
muFar = nan(numel(nGrid),1); muShift = nan(numel(nGrid),1); thr = nan(numel(nGrid),1);
for i = 1:numel(nGrid)
    n = nGrid(i);
    mf = mean(reshape(cFar(randi(numel(cFar), n, nBoot)), n, nBoot), 1);
    seFar(i) = std(mf);  muFar(i) = mean(mf);  thr(i) = quantile(mf, 1-alpha);
    if ~isempty(cShift)
        ms = mean(reshape(cShift(randi(numel(cShift), n, nBoot)), n, nBoot), 1);
        seShift(i) = std(ms);  muShift(i) = mean(ms);
    end
end

% SE = k / sqrt(n)
k = mean(seFar .* sqrt(nGrid(:)), 'omitnan');
kIdeal = 1/sqrt(2);

% ---- the near field, and the power to see it ----------------------------------------------------
cNear = nearS.cos;
obs = NaN; power = nan(numel(nGrid),1);
if ~isempty(cNear)
    obs = mean(cNear, 'omitnan');
    for i = 1:numel(nGrid)
        n = nGrid(i);
        mn = mean(reshape(cNear(randi(numel(cNear), n, nBoot)), n, nBoot), 1);
        power(i) = mean(mn > thr(i));
    end
end

N.curve = table(nGrid(:), seFar, seShift, muFar, muShift, thr, power, ...
    'VariableNames', {'n','seFar','seShift','meanFar','meanShift','thresh','power'});
% The gap between the two nulls is the simultaneous part of the far-field correlation: drift or a
% field-wide flow, as opposed to anything about each track's own shape.
N.driftShare = mean(muFar,'omitnan') - mean(muShift,'omitnan');
N.k = struct('fitted', k, 'ideal', kIdeal, 'penalty', k/kIdeal);
N.far  = struct('r', [rFar Inf], 'nStepPairs', height(farS), 'meanCos', mean(cFar,'omitnan'), ...
                'sdCos', std(cFar,'omitnan'));
N.near = struct('r', [0 rNear], 'nStepPairs', height(nearS), 'meanCos', obs, ...
                'sdCos', std(cNear,'omitnan'));
N.excess = N.near.meanCos - N.far.meanCos;

% steps needed to clear the far null at this alpha, for a few effect sizes
z = norminvApprox(1-alpha);
want = [0.1 0.2 0.3];
if isfinite(N.excess) && N.excess > 0, want = [want N.excess]; end
N.nNeeded = table(want(:), ceil((z*k ./ want(:)).^2), ...
    'VariableNames', {'correlation','steps'});

% Say which way the penalty runs rather than printing a ratio the reader has to interpret.
if k > kIdeal*1.05
    effWord = sprintf('each step is worth about %.2f of an independent one (consecutive steps are correlated)', (kIdeal/k)^2);
elseif k < kIdeal*0.95
    effWord = sprintf('the steps scatter slightly LESS than independent ones (%.2fx), which usually means the far field is not fully mixed', k/kIdeal);
else
    effWord = 'the steps behave as independent draws';
end

nEnough = NaN;
if any(isfinite(power)) && any(power >= 0.8), nEnough = nGrid(find(power >= 0.8, 1)); end
N.nEnough80 = nEnough;

N.text = sprintf([ ...
    'far field (r > %.2f um, %d step pairs) sits at cos = %+.4f; near field (r <= %.2f um, %d pairs) ' ...
    'at %+.4f, an excess of %+.4f. The measured scatter is SE = %.3f/sqrt(n) against the ideal ' ...
    '%.3f/sqrt(n) for independent steps, so %s. ' ...
    'Detecting a correlation of 0.1 needs %d steps, 0.2 needs %d, 0.3 needs %d%s.'], ...
    rFar, height(farS), N.far.meanCos, rNear, height(nearS), N.near.meanCos, N.excess, ...
    k, kIdeal, effWord, N.nNeeded.steps(1), N.nNeeded.steps(2), N.nNeeded.steps(3), ...
    tern(isfinite(nEnough), sprintf('; the observed near-field excess reaches 80%% power at n = %d', nEnough), ...
         '; the observed near-field excess does not reach 80% power at any n tested'));
end

% =================================================================================================
function c = shiftedCos(St)
% The same pairs of tracks, but not the same moment: A's step from row i against B's step from row
% i+lag. Every step vector and every track pairing survives; only the simultaneity is destroyed. So
% whatever correlation is left was never about the two moving at the same time — which is what makes
% this the null that centres on zero even when drift keeps the far field off it.
n = height(St);
if n < 20, c = []; return; end
lag = max(round(n/3), 7);
j = mod((1:n)' + lag - 1, n) + 1;
ua = [St.uax,    St.uay];
ub = [St.ubx(j), St.uby(j)];
la = hypot(ua(:,1), ua(:,2));  lb = hypot(ub(:,1), ub(:,2));
c = sum(ua .* ub, 2) ./ max(la .* lb, eps);
c(la <= 0 | lb <= 0) = NaN;
c = c(isfinite(c));
end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f), v = s.(f); else, v = d; end
end
function y = tern(c,a,b), if c, y=a; else, y=b; end, end

function z = norminvApprox(p)
% One-sided normal quantile without the Statistics toolbox (Acklam's rational approximation).
a = [-3.969683028665376e+01 2.209460984245205e+02 -2.759285104469687e+02 1.383577518672690e+02 ...
     -3.066479806614716e+01 2.506628277459239e+00];
b = [-5.447609879822406e+01 1.615858368580409e+02 -1.556989798598866e+02 6.680131188771972e+01 ...
     -1.328068155288572e+01];
c = [-7.784894002430293e-03 -3.223964580411365e-01 -2.400758277161838e+00 -2.549732539343734e+00 ...
      4.374664141464968e+00 2.938163982698783e+00];
d = [7.784695709041462e-03 3.224671290700398e-01 2.445134137142996e+00 3.754408661907416e+00];
pl = 0.02425;
if p < pl
    q = sqrt(-2*log(p));
    z = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1);
elseif p <= 1-pl
    q = p - 0.5; r = q*q;
    z = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / ...
        (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1);
else
    q = sqrt(-2*log(1-p));
    z = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1);
end
end
