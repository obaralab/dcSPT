function H = dc_pair_panel(D, R, pair, opts)
%DC_PAIR_PANEL  One cross-colour pair, three panels: each colour's intensity, and their co-motion.
%
%   H = dc_pair_panel(D, R, k)                     the k-th row of R.pairs
%   H = dc_pair_panel(D, R, [idA idB], opts)       that pair by track id
%
% Three panels, stacked on a shared time axis so a bleaching step and a change in co-motion line up:
%
%   1  colour A's INTEGRATED intensity per localization, with the photobleaching steps fitted to it
%   2  colour B's, the same
%   3  the co-motion: the per-step cosine, its running mean, and the pair's separation
%
% WHY INTEGRATED AND NOT PEAK. A bleaching step is one fluorophore's photons ceasing. The summed
% intensity over the spot's disk drops by that whole amount; the peak drops by a fraction that
% depends on how well the PSF happened to be centred on a pixel, so a peak trace has extra scatter
% that is nothing to do with the photophysics. dc_pbsa_steps therefore reads iTot.
%
% WHAT THE BLEACHING PANELS ARE FOR HERE. A pair that co-moves is only interesting if both marks are
% single molecules for the stretch being measured: two fluorophores in one diffraction spot move as
% one object by construction, and would read as perfectly correlated with nothing biological behind
% it. The step count says how many emitters each mark had, and WHEN one of them went — so a
% correlation measured over a window where both traces sit at one step is the claim worth making.
%
% READING PANEL 3. The per-step cosine is scatter between -1 and +1 and always looks like noise; the
% running mean is what carries the signal. The shaded band is the far-field null at the same window
% length, from dc_comotion_null — a running mean inside the band is not distinguishable from two
% unrelated molecules in this movie. Separation is drawn on the right axis because co-motion that
% appears only while the two are close is the interesting shape, and co-motion that persists as they
% separate usually means drift that dc_drift did not take out.
%
% OPTIONS
%   .win        [] running-mean window in steps; default the n that dc_comotion_null says is needed
%   .null       [] an N from dc_comotion_null, for the band; omitted, no band is drawn
%   .minStep    'auto'  passed to dc_pbsa_steps
%   .dtS        []  seconds per frame, for a time axis; omitted, the axis is in timepoints
%   .visible    'on'
%
% OUTPUT H: .fig .ax(3) .stepsA .stepsB (the dc_pbsa_steps results) .n .meanCos .rMedian

if nargin < 4 || ~isstruct(opts), opts = struct(); end
vis   = getf(opts,'visible','on');
dtS   = getf(opts,'dtS',[]);
N     = getf(opts,'null',[]);
minSt = getf(opts,'minStep','auto');

% ---- which pair -----------------------------------------------------------------------------
if isscalar(pair)
    assert(pair >= 1 && pair <= height(R.pairs), 'dc_pair_panel:range', ...
        'R.pairs has %d rows, not %d', height(R.pairs), pair);
    idA = R.pairs.trackA(pair); idB = R.pairs.trackB(pair);
else
    idA = min(pair); idB = max(pair);
end
row = R.pairs(R.pairs.trackA == idA & R.pairs.trackB == idB, :);
assert(~isempty(row), 'dc_pair_panel:noPair', ...
    'tracks %g and %g are not a pair in R.pairs (too few shared steps, or a class that was dropped)', idA, idB);
assert(row.class == "cross", 'dc_pair_panel:notCross', ...
    ['tracks %g and %g are a %s-colour pair. This panel is for cross-colour pairs: two tracks of ' ...
     'ONE colour can have had their identities swapped by the linker, and their intensity traces ' ...
     'would then be swapped too, which is exactly what a bleaching step count must not be fed.'], ...
    idA, idB, row.class);

win = getf(opts,'win', []);
if isempty(win)
    win = 25;
    if ~isempty(N) && isfinite(N.nEnough80), win = min(N.nEnough80, max(10, floor(row.n/3))); end
    win = max(5, min(win, max(5, floor(row.n/2))));
end

% ---- the two tracks' localizations, in timepoint order ---------------------------------------
[tA, iA, chA] = trackTrace(D, idA);
[tB, iB, chB] = trackTrace(D, idB);

% ---- the pair's steps -------------------------------------------------------------------------
st = R.steps(R.steps.trackA == idA & R.steps.trackB == idB, :);
st = sortrows(st, 'tp');
assert(~isempty(st), 'dc_pair_panel:noSteps', 'no shared steps for this pair');

xs = @(tp) toX(tp, dtS);
xl = 'timepoint'; if ~isempty(dtS), xl = 'time (s)'; end

H.fig = figure('Color','w','Visible',vis,'Position',[60 60 1040 860], ...
    'Name', sprintf('pair %g-%g', idA, idB));
tl = tiledlayout(H.fig, 3, 1, 'TileSpacing','compact', 'Padding','compact');

% ---- panels 1 and 2: integrated intensity, with bleaching steps --------------------------------
[H.ax(1), H.stepsA] = intensityPanel(tl, xs(tA), iA, chA, minSt, xl);
[H.ax(2), H.stepsB] = intensityPanel(tl, xs(tB), iB, chB, minSt, xl);

% ---- panel 3: co-motion -------------------------------------------------------------------------
ax = nexttile(tl); H.ax(3) = ax; hold(ax,'on');
x = xs(st.tp);
plot(ax, x, st.cos, '.', 'Color', [0.62 0.67 0.72], 'MarkerSize', 5, ...
     'DisplayName', 'per-step cos\theta');
rm = movmean(st.cos, win, 'omitnan');
if ~isempty(N)
    % the far-field null at this window length: everything inside is indistinguishable from two
    % unrelated molecules in this same movie
    [~, j] = min(abs(N.curve.n - win));
    lo = N.curve.meanFar(j) - 2*N.curve.seFar(j);
    hi = N.curve.meanFar(j) + 2*N.curve.seFar(j);
    fill(ax, [x(1) x(end) x(end) x(1)], [lo lo hi hi], [0.85 0.88 0.92], ...
        'EdgeColor','none', 'FaceAlpha', 0.7, 'DisplayName', sprintf('far-field null, n=%d', N.curve.n(j)));
    yline(ax, N.curve.meanFar(j), ':', 'Color', [0.45 0.5 0.55], 'HandleVisibility','off');
end
plot(ax, x, rm, '-', 'Color', [0.80 0.20 0.20], 'LineWidth', 1.6, ...
    'DisplayName', sprintf('running mean, %d steps', win));
yline(ax, 0, '-', 'Color', [0.6 0.6 0.6], 'HandleVisibility','off');
ylabel(ax, 'cos\theta between steps'); ylim(ax, [-1.05 1.05]);

yyaxis(ax, 'right');
plot(ax, x, st.r*1000, '-', 'Color', [0.15 0.45 0.70], 'LineWidth', 1.0, 'DisplayName','separation');
ylabel(ax, 'separation (nm)');
ax.YAxis(2).Color = [0.15 0.45 0.70];
yyaxis(ax, 'left');
xlabel(ax, xl); grid(ax,'on');
% Below the axes, not on top of them: every part of this panel carries data, so 'best' has nowhere
% good to go.
legend(ax, 'Location','southoutside', 'Orientation','horizontal', 'Box','off', 'FontSize',8);
title(ax, sprintf('co-motion: %d shared steps, mean cos %+.3f%s, median separation %.0f nm', ...
    row.n, row.meanCos, nullNote(N, row), 1000*row.rMedian));

title(tl, sprintf('cross-colour pair: %s track %g  vs  %s track %g', chA, idA, chB, idB), ...
    'FontWeight','bold');

H.n = row.n; H.meanCos = row.meanCos; H.rMedian = row.rMedian; H.win = win;
H.trackA = idA; H.trackB = idB;
end

% =================================================================================================
function [ax, St] = intensityPanel(tl, x, y, ch, minSt, xl)
ax = nexttile(tl); hold(ax,'on');
St = struct('k',NaN,'fit',[],'text','no intensity recorded');
if all(isnan(y))
    plot(ax, x, zeros(size(x)), '-', 'Color',[0.8 0.8 0.8]);
    text(ax, 0.5, 0.5, 'no integrated intensity in this dataset', 'Units','normalized', ...
        'HorizontalAlignment','center', 'Color',[0.5 0.5 0.5]);
    ylabel(ax, sprintf('%s  \\Sigma I', ch)); xlabel(ax, xl); grid(ax,'on'); return
end
plot(ax, x, y, '-', 'Color', [0.55 0.60 0.66], 'LineWidth', 0.7);
plot(ax, x, y, '.', 'Color', [0.35 0.40 0.46], 'MarkerSize', 4);
try
    St = dc_pbsa_steps(y, 'MinStep', minSt);
    if isfield(St,'fit') && numel(St.fit) == numel(y)
        stairs(ax, x, St.fit, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6);
    end
catch ME
    St.text = ['step fit failed: ' ME.message];
end
ns = NaN; if isfield(St,'k'), ns = St.k; end
% Build the label before the title: MATLAB evaluates both branches of an inline conditional, so
% reaching for St.text when the fit succeeded (and left no such field) would error.
if isfinite(ns)
    lbl = sprintf('%d bleaching step(s)', ns);
elseif isfield(St,'text')
    lbl = St.text;
else
    lbl = 'no step fit';
end
ylabel(ax, sprintf('%s  \\Sigma I (counts)', ch));
xlabel(ax, xl); grid(ax,'on');
title(ax, sprintf('%s integrated intensity — %s', ch, lbl));
end

function [tp, iTot, ch] = trackTrace(D, id)
m = D.spots.trackId == id;
T = sortrows(D.spots(m,:), 'tp');
tp = T.tp; iTot = T.iTot; ch = char(string(T.ch(1)));
end

function x = toX(tp, dtS)
if isempty(dtS), x = double(tp); else, x = double(tp) * dtS; end
end

function s = nullNote(N, row)
s = '';
if isempty(N), return; end
if isfinite(N.far.meanCos)
    s = sprintf(' (far field %+.3f)', N.far.meanCos);
end
end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
function y = tern(c,a,b), if c, y=a; else, y=b; end, end
