function H = dc_draw(ax, what, varargin)
%DC_DRAW  The drawing every view shares: detection circles, track paths, current marks.
%
%   P = dc_draw([], 'prep', D)                    once, when the dataset changes
%   dc_draw(ax, 'spots',  xy, rPx, colour)
%   dc_draw(ax, 'tracks', D, tpNow, opts)
%   dc_draw(ax, 'pair',   D, idA, idB, tpNow, opts)
%
% ONE PLACE, BECAUSE A CIRCLE DRAWN TWO WAYS IS TWO BUGS. The Detect tab, the Track tab and the Pair
% tab all draw spots on a frame, and if each does it its own way they disagree about what a
% detection looks like the moment one of them is edited.
%
% CIRCLES ARE DRAWN IN DATA COORDINATES, NOT AS MARKERS. A 'o' marker has a size in POINTS: it stays
% the same size on screen when you zoom, so it stops covering the spot it marks and you cannot judge
% whether the detection is centred. A circle built from cos/sin at the spot's real radius is in the
% image's own units, so zooming in makes it bigger with the spot, the way SPTinMatlab draws it — and
% the radius is the true spot radius, which makes the circle a statement about size rather than a
% decoration.
%
% All of them return the graphics handles so a caller can delete just the overlay and redraw it
% without touching the image underneath.
%
% PREPARE ONCE, DRAW MANY TIMES. A player redraws several panels per frame, and each redraw was
% filtering the spot TABLE again: selecting a colour, dropping the untracked, taking a time window
% and sorting. On 34,000 rows that is 16 ms a call and four calls a frame — the whole frame budget,
% spent re-deriving something that only changes when the dataset does. dc_draw([], 'prep', D)
% returns the same columns as plain numeric arrays, pre-sorted; pass that in place of D and the same
% work is about a tenth of a millisecond. Passing D itself still works and prepares on the fly.

switch lower(char(what))
    case 'prep',   H = prep(varargin{:});
    case 'spots',  H = drawSpots(ax, varargin{:});
    case 'tracks', H = drawTracks(ax, varargin{:});
    case 'pair',   H = drawPair(ax, varargin{:});
    otherwise, error('dc_draw:what','unknown mode ''%s''', what);
end
end

% =================================================================================================
function H = drawSpots(ax, xy, rPx, col, style)
% One line object holding every circle, separated by NaN. N separate plot() calls on a frame with a
% few hundred spots is what makes a preview feel slow; this is one child no matter how many spots.
if nargin < 5 || isempty(style), style = '-'; end
H = gobjects(0);
if isempty(xy), return; end
r = max(rPx, 0.75);
th = linspace(0, 2*pi, 25)';
X = xy(:,1)' + r*cos(th);          % 25 x N
Y = xy(:,2)' + r*sin(th);
X(end+1,:) = NaN; Y(end+1,:) = NaN;
H = plot(ax, X(:), Y(:), style, 'Color', col, 'LineWidth', 1.0);
end

function P = prep(D)
% The spot table as plain arrays, sorted by track then timepoint. Numeric indexing on these is two
% orders of magnitude cheaper than the equivalent table operations, and nothing about the data
% changes — only its container.
S = D.spots;
keep = isfinite(S.trackId);
S = S(keep, :);
[~, ord] = sortrows([double(S.trackId), double(S.tp)]);
% Assigned field by field, NOT through struct('f', array): struct() treats a non-scalar value as one
% element per entry and hands back a 1xN struct ARRAY, so P.tid becomes a comma-separated list and
% every later size() on it gets N arguments instead of one.
P = struct('isPrep', true);
P.ch     = double(S.ch(ord));
P.chName = categories(S.ch);
P.tid    = double(S.trackId(ord));
P.tp     = double(S.tp(ord));
P.x      = double(S.x(ord));
P.y      = double(S.y(ord));
end

function H = drawTracks(ax, D, tpNow, opts)
% Every track of one colour, as a path up to tpNow, with the current position marked. A track that
% has no localization at tpNow is drawn faded: it exists in the movie but not in this frame, and
% hiding it would make the field look emptier than it is.
if nargin < 4 || ~isstruct(opts), opts = struct(); end
ch      = getf(opts,'ch','');
colPath = getf(opts,'colour',[1 0.3 1]);
tail    = getf(opts,'tail',Inf);        % how many timepoints of history to draw
rPx     = getf(opts,'rPx',3);
pxUm    = getf(opts,'pxUm',1);          % positions are in um; the image is in px
only    = getf(opts,'tracks',[]);       % restrict to these track ids
H = struct('paths',gobjects(0), 'now',gobjects(0), 'n',0);

if isstruct(D) && isfield(D,'isPrep'), P = D; else, P = prep(D); end
m = true(size(P.tid));
if ~isempty(ch)
    k = find(strcmp(P.chName, char(ch)), 1);
    if isempty(k), return; end
    m = m & (P.ch == k);
end
if ~isempty(only), m = m & ismember(P.tid, only); end
m = m & P.tp <= tpNow & P.tp > tpNow - tail;
if ~any(m), return; end

tid = P.tid(m); X = P.x(m)/pxUm; Y = P.y(m)/pxUm; tp = P.tp(m);

% One polyline for all paths, NaN between tracks. The arrays arrive sorted by track then timepoint,
% so a track boundary is just a change in tid and no sort is needed here.
brk = [find(diff(tid) ~= 0); numel(tid)];
n = numel(tid) + numel(brk);
xs = nan(n,1); ys = nan(n,1);
last = 0; o = 0;
for b = brk'
    len = b - last;
    xs(o+1:o+len) = X(last+1:b); ys(o+1:o+len) = Y(last+1:b);
    o = o + len + 1;                 % leave the NaN that separates this track from the next
    last = b;
end
H.paths = plot(ax, xs, ys, '-', 'Color', [colPath 0.75], 'LineWidth', 1.0);
H.n = numel(brk);

cur = tp == tpNow;
if any(cur)
    H.now = drawSpots(ax, [X(cur), Y(cur)], rPx, colPath);
end
end

function H = drawPair(ax, D, idA, idB, tpNow, opts)
% The two tracks of one cross-colour pair, each in its channel's colour, with a line joining them at
% the current timepoint. That connector is the whole point of the view: it is the separation the
% co-motion statistic is binned by, drawn where you can see whether it is real or a coincidence of
% two unrelated molecules passing.
if nargin < 6 || ~isstruct(opts), opts = struct(); end
pxUm = getf(opts,'pxUm',1);
rPx  = getf(opts,'rPx',3);
cA   = getf(opts,'colourA',[1 0.25 1]);
cB   = getf(opts,'colourB',[0.25 1 0.35]);
tail = getf(opts,'tail',Inf);
H = struct('a',[], 'b',[], 'link',gobjects(0), 'rNm',NaN);

H.a = drawTracks(ax, D, tpNow, struct('tracks',idA,'colour',cA,'pxUm',pxUm,'rPx',rPx,'tail',tail));
H.b = drawTracks(ax, D, tpNow, struct('tracks',idB,'colour',cB,'pxUm',pxUm,'rPx',rPx,'tail',tail));

if isstruct(D) && isfield(D,'isPrep')
    pa = prepPos(D, idA, tpNow); pb = prepPos(D, idB, tpNow);
else
    pa = posAt(D, idA, tpNow); pb = posAt(D, idB, tpNow);
end
if ~isempty(pa) && ~isempty(pb)
    H.link = plot(ax, [pa(1) pb(1)]/pxUm, [pa(2) pb(2)]/pxUm, '-', ...
        'Color', [1 1 1 0.85], 'LineWidth', 1.2);
    H.rNm = 1000*hypot(pa(1)-pb(1), pa(2)-pb(2));
end
end

function p = prepPos(P, id, tp)
m = P.tid == id & P.tp == tp;
p = []; if any(m), k = find(m,1); p = [P.x(k) P.y(k)]; end
end

function p = posAt(D, id, tp)
m = D.spots.trackId == id & D.spots.tp == tp;
p = [];
if any(m), r = find(m,1); p = [D.spots.x(r) D.spots.y(r)]; end
end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
