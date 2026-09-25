function dc_comotion_smoke()
%DC_COMOTION_SMOKE  Co-motion between nearby tracks, against data whose answer is known in advance.
%
% WHAT IS ASSERTED:
%   1. INDEPENDENT TRACKS READ AS ZERO. Nothing coupled, no drift: the mean cosine over every
%      separation sits within its own error bar of zero. If this fails nothing else means anything.
%   2. INJECTED CO-MOTION IS RECOVERED AT ITS TRUE STRENGTH BY THE PAIR STATISTIC, and only between
%      the pairs it was injected into. A per-component correlation of 0.6 must read as a mean cosine
%      near (pi/4)*0.6 = 0.47 on the coupled pairs themselves. The BINNED curve reads far lower
%      (~0.13) because most pairs that happen to be close by are not coupled at all, and the test
%      pins that gap deliberately: the curve says whether short-range co-motion exists, the pairs say
%      how strong it is, and confusing the two underestimates the coupling.
%   3. DRIFT IS THE FALSE POSITIVE, AND IT IS CAUGHT AND REMOVED. A field-wide drift makes every
%      pair correlated at EVERY separation — the signature is a flat, raised curve rather than one
%      that decays with distance. dc_drift removes it and the curve returns to zero.
%   4. THE FAR-FIELD NULL ABSORBS WHAT DRIFT CORRECTION MISSES. Even with drift left in, the near
%      field minus the far field recovers the injected coupling, because both carry the same drift.
%      This is the whole reason for an empirical null rather than comparing against zero.
%   5. THE TWO NULLS DECOMPOSE THE FAR FIELD, and neither is "the clean one". Time-shifting removes
%      only TIME-VARYING common motion: against a STEADY drift the shifted null stays just as high as
%      the far one, because every step shares the same mean vector whatever moment it was taken. Give
%      the stage a random walk instead and the shifted null drops away while the far one does not.
%      Both directions are asserted, because assuming the shifted null always centres on zero is
%      exactly the mistake that would make a steady drift look like real co-motion.
%   6. MORE STEPS, LESS NOISE, AT THE RATE MEASURED NOT ASSUMED. The null SE falls as k/sqrt(n), and
%      the n the tool reports as sufficient really does separate coupled from uncoupled pairs.
%   7. CROSS-COLOUR ONLY IS THE DEFAULT, and it is applied before anything is pooled, so a count in
%      .bins counts what was used. Same-colour pairs are reachable but not by accident.
%   8. THE PAIR PANEL reads integrated intensity, counts the bleaching steps that were injected into
%      it, and refuses a same-colour pair — whose two intensity traces could have been swapped by the
%      linker, which is exactly what a step count must not be handed.
%   9. SPANS ARE NOT MIXED. A step over two timepoints is not simultaneous with one over a single
%      timepoint, and pairing them is refused rather than silently averaged.
%
% Synthetic; reads no dataset.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'core'), fullfile(fileparts(here),'drivers'), here);
rng(11);

dt = 0.012; Dc = 0.05; sd = sqrt(2*Dc*dt);      % ~35 nm per step
nT = 400;                                        % timepoints
nPair = 14;                                      % coupled pairs
nFree = 28;                                      % independent tracks

%% (1) nothing coupled, no drift -------------------------------------------------------------------
D0 = makeSet(nT, 0, nPair, nFree, sd, 0, 0);
S0 = dc_steps(D0);
R0 = dc_comotion(S0, struct('rMaxUm', 3, 'nMin', 20));
nearBin = R0.bins(R0.bins.class=="all" & R0.bins.rHi <= 0.5, :);
m = sum(nearBin.meanCos .* nearBin.n) / sum(nearBin.n);
se = sqrt(sum((nearBin.seCos.^2) .* nearBin.n.^2)) / sum(nearBin.n);
assert(abs(m) < max(4*se, 0.03), ...
    'uncoupled tracks must read as zero co-motion, got %+.4f +- %.4f', m, se);
fprintf('(1) independent tracks: near-field cos = %+.4f +- %.4f\n', m, se);

%% (2) injected co-motion --------------------------------------------------------------------------
cInj = 0.6;                                      % the shared fraction of each coupled step
[D1, truth] = makeSet(nT, cInj, nPair, nFree, sd, 0, 0);
S1 = dc_steps(D1);
R1 = dc_comotion(S1, struct('rMaxUm', 3, 'nMin', 20));
nb = R1.bins(R1.bins.class=="all" & R1.bins.rHi <= 0.5, :);
mNear = sum(nb.meanCos .* nb.n) / sum(nb.n);
fb = R1.bins(R1.bins.class=="all" & R1.bins.rLo >= 1.5, :);
mFar = sum(fb.meanCos .* fb.n) / sum(fb.n);
assert(mNear > 0.08, 'injected co-motion must show up near, got %+.4f', mNear);
assert(mNear - mFar > 0.07, ...
    'and it must be a NEAR-field effect: near %+.4f vs far %+.4f', mNear, mFar);
fprintf('(2) injected %.2f: near %+.4f, far %+.4f, excess %+.4f\n', cInj, mNear, mFar, mNear-mFar);

% the coupled pairs, and only those, come out on top
P = R1.pairs;
isC = ismember([P.trackA P.trackB], [truth.trackA truth.trackB], 'rows');
assert(any(isC) && any(~isC), ...
    'the fixture should contain both coupled (%d of %d found) and uncoupled pairs', nnz(isC), height(truth));
mC = median(P.meanCos(isC)); mU = median(P.meanCos(~isC));
assert(mC - mU > 0.2, ...
    'a coupled pair must score above an uncoupled one in the same movie (medians %+.3f vs %+.3f)', mC, mU);
% the pair statistic must recover the strength that was injected, not merely rank the pairs
assert(abs(mC - pi/4*cInj) < 0.12, ...
    ['a coupled pair must read the coupling it was given: expected cos near (pi/4)*%.2f = %.3f, ' ...
     'got %.3f'], cInj, pi/4*cInj, mC);
% ...and the binned curve must sit WELL BELOW it, because chance encounters dilute the bin
assert(mNear < 0.6*mC, ...
    ['the near-field bin must be diluted relative to the coupled pairs themselves (bin %+.3f vs ' ...
     'pairs %+.3f) — reading coupling strength off the curve underestimates it'], mNear, mC);
fprintf('    coupled pairs median cos %+.3f (injected rho %.2f -> expected %.3f), uncoupled %+.3f\n', ...
    mC, cInj, pi/4*cInj, mU);
fprintf('    the r<=0.5 bin reads %+.3f: about %.0f%% of near pairs are chance encounters\n', ...
    mNear, 100*(1 - mNear/mC));

%% (3) drift is a false positive, and it is removable ------------------------------------------------
drift = 0.6 * sd;                                 % a steady push, smaller than a typical step
D2 = makeSet(nT, 0, nPair, nFree, sd, drift, 0);  % NO coupling, drift only
S2 = dc_steps(D2);
Rd = dc_comotion(S2, struct('rMaxUm', 3, 'nMin', 20));
bAll = Rd.bins(Rd.bins.class=="all", :);
assert(mean(bAll.meanCos) > 0.08, ...
    'drift alone should make every pair look correlated, got mean %+.4f', mean(bAll.meanCos));
flat = max(bAll.meanCos) - min(bAll.meanCos);
assert(flat < 0.20, 'and it should be FLAT in separation, not decaying (range %.3f)', flat);

[S2c, Dr] = dc_drift(S2);
assert(Dr.medianSpeed > 0.3*drift, 'the drift estimate should find it (%.1f nm/step)', 1000*Dr.medianSpeed);
R2 = dc_comotion(S2c, struct('rMaxUm', 3, 'nMin', 20));
b2 = R2.bins(R2.bins.class=="all", :);
assert(abs(mean(b2.meanCos)) < 0.05, ...
    'after dc_drift the curve must come back to zero, got %+.4f', mean(b2.meanCos));
fprintf('(3) drift %.0f nm/step: raised every bin to %+.4f (flat, range %.3f); after correction %+.4f\n', ...
    1000*drift, mean(bAll.meanCos), flat, mean(b2.meanCos));

%% (4) the far-field null absorbs the drift that is left in -------------------------------------------
D3 = makeSet(nT, cInj, nPair, nFree, sd, drift, 0);   % coupling AND drift, uncorrected
S3 = dc_steps(D3);
R3 = dc_comotion(S3, struct('rMaxUm', 3, 'nMin', 20));
N3 = dc_comotion_null(R3, struct('rNearUm', 0.5, 'rFarUm', 1.5, 'nBoot', 600));
assert(N3.far.meanCos > 0.05, ...
    'with drift left in, the far field must NOT sit at zero — that is the point (%+.4f)', N3.far.meanCos);
assert(N3.excess > 0.10, ...
    'and near minus far must still recover the coupling (%+.4f)', N3.excess);
fprintf('(4) drift left in: far %+.4f (not zero), near %+.4f, excess %+.4f\n', ...
    N3.far.meanCos, N3.near.meanCos, N3.excess);

%% (5) the time-shifted null centres on zero even so --------------------------------------------------
muF = mean(N3.curve.meanFar, 'omitnan');
muS = mean(N3.curve.meanShift, 'omitnan');
assert(all(N3.curve.seShift > 0), 'the time-shifted null should have been computed with real scatter');
% STEADY drift: shifting in time cannot remove it, so the two nulls agree and driftShare is ~0
assert(muF > 0.05, 'the far-field null should carry the steady drift (%+.4f)', muF);
assert(abs(muF - muS) < 0.04, ...
    ['against a STEADY drift the shifted null must stay just as high (%+.4f vs %+.4f): every step ' ...
     'shares the same mean vector whatever moment it came from'], muF, muS);
assert(abs(N3.driftShare) < 0.04, ...
    'so driftShare should be near zero for a steady drift, got %+.4f', N3.driftShare);
fprintf('(5a) steady drift: far %+.4f, shifted %+.4f, driftShare %+.4f (shifting cannot remove it)\n', ...
    muF, muS, N3.driftShare);

% A WANDERING stage: now the common motion IS time-varying, so shifting destroys it
Dw = makeSet(nT, 0, nPair, nFree, sd, 1.2*sd, 0, 'walk');
Rw = dc_comotion(dc_steps(Dw), struct('rMaxUm', 3, 'nMin', 20));
Nw = dc_comotion_null(Rw, struct('rNearUm',0.5,'rFarUm',1.5,'nBoot',600));
wF = mean(Nw.curve.meanFar,'omitnan'); wS = mean(Nw.curve.meanShift,'omitnan');
assert(wF > 0.05, 'a wandering stage should still correlate the far field (%+.4f)', wF);
assert(wF - wS > 0.04, ...
    ['and NOW shifting should remove most of it (far %+.4f vs shifted %+.4f) — this is the case ' ...
     'driftShare is built to detect'], wF, wS);
fprintf('(5b) wandering stage: far %+.4f, shifted %+.4f, driftShare %+.4f (shifting removes it)\n', ...
    wF, wS, Nw.driftShare);

%% (6) the n curve --------------------------------------------------------------------------------
assert(issorted(N3.curve.seFar, 'descend'), 'the null SE must fall as n grows');
r = corr_(log(N3.curve.n), log(N3.curve.seFar));
assert(r < -0.95, 'and fall as a power law in n (log-log r = %.3f)', r);
assert(N3.k.fitted > 0.3 && N3.k.fitted < 1.5, 'the fitted k should be order 1/sqrt(2), got %.3f', N3.k.fitted);
assert(height(N3.nNeeded) >= 3 && all(N3.nNeeded.steps > 0), 'nNeeded must be filled');
big = N3.nNeeded.steps(N3.nNeeded.correlation == 0.1);
small = N3.nNeeded.steps(N3.nNeeded.correlation == 0.3);
assert(big > small, 'a weaker correlation must need MORE steps (%d vs %d)', big, small);
fprintf('(6) SE = %.3f/sqrt(n) (ideal %.3f, penalty %.2fx); 0.1 needs %d steps, 0.3 needs %d\n', ...
    N3.k.fitted, N3.k.ideal, N3.k.penalty, big, small);
fprintf('    %s\n', N3.text);

%% (7) cross-colour is the default -----------------------------------------------------------------
Rx = dc_comotion(S1, struct('rMaxUm', 3, 'nMin', 20));           % no classes given
assert(Rx.class == "cross" && all(Rx.steps.class == "cross"), ...
    'cross-colour must be the default class, got %s', Rx.class);
assert(all(Rx.steps.chA ~= Rx.steps.chB), 'and every kept pair must span the two colours');
Rall = dc_comotion(S1, struct('rMaxUm', 3, 'nMin', 20, 'classes', "all"));
assert(height(Rall.steps) > height(Rx.steps), ...
    'classes="all" must keep more (%d vs %d)', height(Rall.steps), height(Rx.steps));
nb2 = Rx.bins(Rx.bins.class=="all" & Rx.bins.rHi <= 0.5, :);
assert(sum(nb2.n) == nnz(Rx.steps.r <= 0.5), ...
    'the bins must count only what survived the class filter');
fprintf('(7) default kept %d cross pairs of %d total; near-field cross cos %+.4f\n', ...
    height(Rx.steps), height(Rall.steps), sum(nb2.meanCos.*nb2.n)/sum(nb2.n));

%% (8) the pair panel -------------------------------------------------------------------------------
% A dataset with intensity: each track's trace bleaches in a known number of steps.
[Db, truthB] = makeSet(nT, cInj, nPair, nFree, sd, 0, 0);
nStepsInj = 2;                                    % two emitters -> two drops
Db = addIntensity(Db, nStepsInj);
Sb = dc_steps(Db);
Rb = dc_comotion(Sb, struct('rMaxUm', 3, 'nMin', 20));
Nb = dc_comotion_null(Rb, struct('rNearUm',0.5,'rFarUm',1.5,'nBoot',400));

crossTruth = truthB(truthB.chA ~= truthB.chB, :);
assert(~isempty(crossTruth), 'the fixture must contain cross-colour coupled pairs');
k = find(ismember([Rb.pairs.trackA Rb.pairs.trackB], ...
                  [crossTruth.trackA crossTruth.trackB], 'rows'), 1);
assert(~isempty(k), 'a coupled cross-colour pair should be in R.pairs');

H = dc_pair_panel(Db, Rb, k, struct('visible','off', 'null', Nb, 'dtS', dt));
assert(numel(H.ax) == 3 && all(isgraphics(H.ax)), 'three panels');
assert(H.stepsA.k == nStepsInj && H.stepsB.k == nStepsInj, ...
    ['both intensity traces must yield the %d bleaching steps that were injected (got %d and %d) ' ...
     '— if this drifts, the panel is drawing a staircase nobody checked'], ...
    nStepsInj, H.stepsA.k, H.stepsB.k);
fprintf('(8) panel: %d shared steps, mean cos %+.3f, %d/%d bleaching steps recovered\n', ...
    H.n, H.meanCos, H.stepsA.k, H.stepsB.k);
close(H.fig);

% a same-colour pair must be refused rather than silently drawn
sameTruth = truthB(truthB.chA == truthB.chB, :);
if ~isempty(sameTruth)
    Rs = dc_comotion(Sb, struct('rMaxUm',3, 'nMin',20, 'classes',"all"));
    ks = find(ismember([Rs.pairs.trackA Rs.pairs.trackB], ...
                       [sameTruth.trackA sameTruth.trackB], 'rows'), 1);
    if ~isempty(ks)
        err = ''; try, dc_pair_panel(Db, Rs, ks, struct('visible','off')); catch ME, err = ME.identifier; end
        assert(strcmp(err,'dc_pair_panel:notCross'), ...
            'a same-colour pair must be refused by the panel, got "%s"', err);
        fprintf('    same-colour pair refused, as it must be\n');
    end
end

%% (9) spans are not mixed -----------------------------------------------------------------------
Dg = makeSet(nT, cInj, nPair, nFree, sd, 0, 0.12);   % 12%% of localizations dropped -> real gaps
Smix = dc_steps(Dg, struct('span', []));
if numel(unique(Smix.span)) > 1
    err = ''; try, dc_comotion(Smix); catch ME, err = ME.identifier; end
    assert(strcmp(err,'dc_comotion:mixedSpan'), ...
        'pairing steps of different spans must be refused, got "%s"', err);
    fprintf('(9) mixed spans refused: %s\n', mat2str(unique(Smix.span)'));
else
    fprintf('(9) the fixture produced no gaps; span guard not exercised\n');
end

fprintf('\nDC-COMOTION SMOKE PASSED.\n');
end

% =================================================================================================
function [D, truth] = makeSet(nT, coupling, nPair, nFree, sd, driftPerStep, gapProb, driftMode)
% A field of 2D random walks. A coupled pair shares `coupling` of each step's variance and starts
% close together; free tracks are independent and scattered.
%
% The fixture RETURNS ITS OWN GROUND TRUTH — `truth`, the global ids of the coupled pairs — rather
% than leaving the test to re-derive them from the numbering. dc_dataset renumbers track ids per
% colour as it adds them, so any rule the test guessed would be a second, silent copy of that
% renumbering, and the first thing to break when it changes.
%
% Half the coupled pairs are put in ONE colour and half are split ACROSS the two, so both the
% swap-prone 'same' class and the swap-immune 'cross' class carry real signal.
fov = 8;
tracks = {}; colour = [];                      % colour(k) = 1 or 2 for tracks{k}
pairIdx = zeros(nPair, 2);                     % the two entries of tracks{} for each coupled pair
for p = 1:nPair
    c0 = 1 + (fov-2)*rand(1,2);
    off = 0.15*(rand(1,2)-0.5);                % the pair starts within ~150 nm
    shared = sd*sqrt(coupling)*randn(nT-1, 2);
    tracks{end+1} = cumsum([c0;       shared + sd*sqrt(1-coupling)*randn(nT-1,2)], 1); %#ok<AGROW>
    pairIdx(p,1) = numel(tracks);
    tracks{end+1} = cumsum([c0 + off; shared + sd*sqrt(1-coupling)*randn(nT-1,2)], 1); %#ok<AGROW>
    pairIdx(p,2) = numel(tracks);
    if mod(p,2) == 1, colour(pairIdx(p,:)) = [1 1];     % same-colour pair
    else,             colour(pairIdx(p,:)) = [1 2];     % cross-colour pair
    end
end
for k = 1:nFree
    tracks{end+1} = cumsum([1 + (fov-2)*rand(1,2); sd*randn(nT-1,2)], 1); %#ok<AGROW>
    colour(numel(tracks)) = 1 + mod(k,2);
end
if nargin < 8 || isempty(driftMode), driftMode = 'steady'; end
if driftPerStep ~= 0
    switch driftMode
        case 'steady'                                   % the same push every frame
            step = repmat(driftPerStep*[0.8 0.6], nT-1, 1);
        case 'walk'                                     % a stage that wanders: common but time-varying
            step = driftPerStep * randn(nT-1, 2);
        otherwise, error('unknown driftMode %s', driftMode);
    end
    d = cumsum([0 0; step], 1);
    for k = 1:numel(tracks), tracks{k} = tracks{k} + d; end
end

keys = {'c1','c2'};
C = struct('key',{},'label',{},'pages',{},'tp',{},'dt_s',{},'nFrames',{});
for i = 1:2
    C(i) = struct('key',keys{i},'label',keys{i},'pages',(1:nT)','tp',(1:nT)','dt_s',0.012,'nFrames',nT);
end
D = dc_dataset('new','synthetic',C);

globalOf = zeros(1, numel(tracks));            % tracks{} index -> global track id
offset = 0;
for i = 1:2
    mine = find(colour == i);
    fr = []; x = []; y = []; tid = [];
    for q = 1:numel(mine)
        t = tracks{mine(q)};
        keep = true(nT,1);
        if gapProb > 0, keep = rand(nT,1) > gapProb; keep([1 end]) = true; end
        f = find(keep) - 1;
        fr = [fr; f]; x = [x; t(keep,1)]; y = [y; t(keep,2)]; %#ok<AGROW>
        tid = [tid; repmat(q-1, nnz(keep), 1)]; %#ok<AGROW>
        globalOf(mine(q)) = offset + q - 1;     % dc_dataset offsets by the ids already present
    end
    R = struct('key',keys{i}, 'frame',fr, 'x',x, 'y',y, 'q',100+0*fr, ...
        'trackId',tid, 'spotId',(0:numel(fr)-1)', 'tracks',{repmat({zeros(0,3)},1,numel(mine))}, ...
        'nFrames',nT, 'nTracks',numel(mine), 'nDets',numel(fr));
    D = dc_dataset('addChannel', D, R);
    offset = offset + numel(mine);
end

a = globalOf(pairIdx(:,1))';  b = globalOf(pairIdx(:,2))';
truth = table(min(a,b), max(a,b), string(keys(colour(pairIdx(:,1)))'), ...
              string(keys(colour(pairIdx(:,2)))'), ...
              'VariableNames', {'trackA','trackB','chA','chB'});
% dc_comotion keys a pair as (lower id, higher id); match that so a join just works
assert(all(truth.trackA < truth.trackB), 'coupled pairs must be keyed low-to-high');
end

function D = addIntensity(D, nSteps)
% Give every track an intensity trace that bleaches in nSteps discrete drops, so the panel's step
% count can be checked against a number that was put there on purpose.
S = D.spots;
lv = 1200;                                    % counts per emitter
S.iTot = nan(height(S),1);
for id = unique(S.trackId)'
    if ~isfinite(id), continue; end
    m = find(S.trackId == id);
    n = numel(m);
    when = round(linspace(0, n, nSteps+2)); when = when(2:end-1);   % where the drops fall
    emitters = nSteps + 1 - sum((1:n)' > when, 2);
    S.iTot(m) = lv*emitters + 18*randn(n,1);  % noise well under a step, so the fit is unambiguous
end
D.spots = S;
end

function r = corr_(x, y)
x = x(:) - mean(x); y = y(:) - mean(y);
r = sum(x.*y) / sqrt(sum(x.^2)*sum(y.^2));
end
