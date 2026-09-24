function [S, info] = dc_merge(D, opts)
%DC_MERGE  Both colours at one TIMEPOINT — the view the merged panel draws.
%
%   S = dc_merge(D)                                 every spot of both colours, flagged
%   [S, info] = dc_merge(D, struct('tp', k))        what each colour has at timepoint k
%   ... struct('frame', f, 'ref', 'c2')             that frame of c2, and whatever shares its timepoint
%
% THE TIMEPOINT IS THE INDEX. A frame query is resolved through the reference colour's own map into
% a timepoint, and every other colour is taken at THAT timepoint — not at its own frame of the same
% number, which is a different moment whenever the colours are not matched (see dc_align).
%
% WHAT IT REFUSES TO DO. A colour with no frame at that timepoint contributes NOTHING, and info says
% so. It does not fall back to that colour's nearest frame: a mark drawn from a neighbouring
% timepoint sits beside a current one and looks like a coincidence that was never observed. Where
% the colours are subsampled this is the normal case, not an error — half the timepoints simply have
% one colour in them.
%
% OUTPUT
%   S     the spot table, filtered: .ch .frame .tp .page .x .y .q .trackId .trackLocal .spotId
%   info  .tp .perChannel(.key .frame .page .present .n) .missing .text

if nargin < 2 || ~isstruct(opts), opts = struct(); end
S = D.spots;
info = struct('tp',NaN, 'perChannel',struct('key',{},'frame',{},'page',{},'present',{},'n',{}), ...
              'missing',{{}}, 'text','all spots, both colours');
if ~isfield(opts,'tp') && ~isfield(opts,'frame'), return; end

if isfield(opts,'frame') && ~isempty(opts.frame)
    ref = char(D.channels(1).key);
    if isfield(opts,'ref') && ~isempty(opts.ref), ref = char(opts.ref); end
    c = chanOf(D, ref);
    f = double(opts.frame);
    assert(f >= 0 && f < c.nFrames, 'dc_merge:frameRange', ...
        'colour ''%s'' has frames 0..%d, not %d', ref, c.nFrames-1, f);
    tp = c.tp(f+1);
else
    tp = double(opts.tp);
end
info.tp = tp;

keep = false(height(S),1);
for i = 1:numel(D.channels)
    c = D.channels(i);
    f = find(c.tp == tp, 1);                 % this colour's frame at that timepoint, if it has one
    present = ~isempty(f);
    fr = NaN; pg = NaN; m = false(height(S),1);
    if present
        fr = f - 1;                          % 0-based
        pg = c.pages(f);
        m = (S.ch == c.key) & (S.frame == fr);
    end
    keep = keep | m;
    info.perChannel(end+1) = struct('key',c.key, 'frame',fr, 'page',pg, ...
                                    'present',present, 'n',nnz(m)); %#ok<AGROW>
    if ~present, info.missing{end+1} = c.key; end %#ok<AGROW>
end
S = S(keep,:);

parts = arrayfun(@(p) sprintf('%s frame %g (page %g, %d spots)', p.key, p.frame, p.page, p.n), ...
                 info.perChannel, 'uni', 0);
info.text = sprintf('timepoint %d — %s', tp, strjoin(parts, ' · '));
if ~isempty(info.missing)
    info.text = sprintf('timepoint %d — %s · %s has no frame at this timepoint', tp, ...
        strjoin(parts(~ismember({info.perChannel.key}, info.missing)), ' · '), ...
        strjoin(info.missing, ' and '));
end
end

% =================================================================================================
function c = chanOf(D, key)
k = find(strcmp(arrayfun(@(x) char(x.key), D.channels(:)', 'uni', 0), key), 1);
assert(~isempty(k), 'dc_merge:noChannel', 'no colour ''%s''', key);
c = D.channels(k);
end
