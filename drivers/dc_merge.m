function [S, info] = dc_merge(D, opts)
%DC_MERGE  Both colours in one list, at one instant — the view the merged panel draws.
%
%   S = dc_merge(D)                              every spot of both colours, flagged
%   [S, info] = dc_merge(D, struct('t_s', t))    what was visible at that instant
%   ... struct('t_s', t, 'tol_s', tol)           how near in time counts as "then"
%   ... struct('frame', k, 'ref', 'a')           that frame of colour a, and whatever of b lines up
%
% TIME IS THE INDEX, NOT THE FRAME NUMBER. Frame k of one colour is frame k of the other only when
% the two share a clock (see dc_align); in every other case stepping by index silently shows two
% instants at once. So a query by frame is resolved through the reference colour's clock into a
% time, and the other colour is taken from that time.
%
% WHAT IT REFUSES TO DO. When a colour has no frame within tol_s of the instant asked for, its spots
% are ABSENT and info says so — it does not fall back to that colour's last known positions. Drawing
% a stale position beside a current one is how a merged view invents colocalization: the two marks
% sit together on screen because one of them is old, not because the molecules met.
%
% OUTPUT
%   S     the spot table, filtered: .ch .frame .t_s .x .y .q .trackId .spotId
%   info  .t_s        the instant resolved to
%         .perChannel one row per colour: .key .frame .t_s .gap_s .present .n
%         .missing    the colours with nothing near that instant
%         .text       one sentence for the panel

if nargin < 2 || ~isstruct(opts), opts = struct(); end
S = D.spots;
K = arrayfun(@(c) char(c.key), D.channels(:)', 'uni', 0);
info = struct('t_s',NaN, 'perChannel',struct('key',{},'frame',{},'t_s',{},'gap_s',{},'present',{},'n',{}), ...
              'missing',{{}}, 'text','all spots, both colours');
if isempty(fieldnames(opts)) || (~isfield(opts,'t_s') && ~isfield(opts,'frame')), return; end

% ---- resolve the instant ----
if isfield(opts,'frame') && ~isempty(opts.frame)
    ref = K{1}; if isfield(opts,'ref') && ~isempty(opts.ref), ref = char(opts.ref); end
    c = chanOf(D, ref);
    t = c.t0_s + double(opts.frame)*c.dt_s;
else
    t = double(opts.t_s);
end
info.t_s = t;

% Default tolerance: half the COARSEST frame interval. Anything nearer than that was the nearest
% thing imaged; anything further away is a different moment.
tol = max(arrayfun(@(c) c.dt_s, D.channels))/2;
if isfield(opts,'tol_s') && ~isempty(opts.tol_s), tol = double(opts.tol_s); end

keep = false(height(S),1);
for i = 1:numel(D.channels)
    c = D.channels(i);
    f = round((t - c.t0_s)/c.dt_s);                 % this colour's nearest frame to that instant
    f = max(f, 0);
    tf = c.t0_s + f*c.dt_s;
    gap = tf - t;
    present = abs(gap) <= tol*(1+1e-9) + 1e-12;
    m = present & (S.ch == c.key) & (S.frame == f);
    keep = keep | m;
    info.perChannel(end+1) = struct('key',c.key, 'frame',f, 't_s',tf, 'gap_s',gap, ...
                                    'present',present, 'n',nnz(m)); %#ok<AGROW>
    if ~present, info.missing{end+1} = c.key; end %#ok<AGROW>
end
S = S(keep,:);

parts = arrayfun(@(p) sprintf('%s frame %d (%+.1f ms, %d spots)', p.key, p.frame, 1000*p.gap_s, p.n), ...
                 info.perChannel, 'uni', 0);
info.text = sprintf('t = %.4f s — %s', t, strjoin(parts, ' · '));
if ~isempty(info.missing)
    info.text = sprintf('%s — %s not imaged within %.1f ms of this instant', ...
        info.text, strjoin(info.missing, ' and '), 1000*tol);
end
end

% =================================================================================================
function c = chanOf(D, key)
k = find(strcmp(arrayfun(@(x) char(x.key), D.channels(:)', 'uni', 0), key), 1);
assert(~isempty(k), 'dc_merge:noChannel', 'no colour ''%s''', key);
c = D.channels(k);
end
