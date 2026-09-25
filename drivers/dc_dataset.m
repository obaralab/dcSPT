function D = dc_dataset(mode, varargin)
%DC_DATASET  One cell, both colours — indexed by integers the acquisition itself supplied.
%
%   D = dc_dataset('new', cellKey, channels)
%   D = dc_dataset('addChannel', D, R)
%   D = dc_dataset('validate', D)
%
% See docs/DATA_MODEL.md. Per-channel arrays stay per channel, because the two colours do not share
% a frame axis. The merged view is a flat spot list in which the shared index is a COLUMN — and that
% index is the acquisition TIMEPOINT, not a time in seconds.
%
% WHY NOT SECONDS. A second is a frame number multiplied by a nominal interval; these stacks carry
% one `finterval` for the whole file and no per-page timestamps at all, so every second is derived
% and every derivation accumulates — the 1.1% between 0.054 and 0.0534237 is half a second across a
% 1,791-frame track. A timepoint is what the microscope recorded. It is exact, and it is also the
% unit a tracking gap is naturally counted in: frame k+2 after frame k is a gap of one, with nothing
% to compare against a tolerance.
%
% THE INVARIANTS, all checked by 'validate':
%   1. every spot's channel is a key in the registry
%   2. a spot's page and timepoint are EXACTLY what the channel map says for its frame — integer
%      equality, not a tolerance, because there is nothing approximate about an index
%   3. frames are 0-based and within that channel's range
%   4. a track id belongs to exactly one channel: tracking is per colour
%   5. spot ids are unique within a channel
%   6. every spot's local track index points at a track that colour actually has
%
% TRACK IDS ARE MADE GLOBAL ON THE WAY IN. Each colour is tracked on its own, so each colour's
% tracker numbers from zero and 'track 0' exists in both — the ids collide for the most ordinary
% reason there is. Renumbering on add means one number identifies one track anywhere in the cell, so
% a list, an export or a selection never has to carry the colour beside it to be unambiguous. The
% colour's OWN index is kept as .trackLocal, because a track has to be findable in the .tracks cell
% the tracker returned; .idOffset on the per-colour record is the mapping between them, recorded
% rather than left to be re-derived.

switch lower(char(mode))
    case 'new',        D = newD(varargin{:});
    case 'addchannel', D = addChannel(varargin{:});
    case 'validate',   D = validate(varargin{:});
    otherwise, error('dc_dataset:mode','unknown mode ''%s''', mode);
end
end

% =================================================================================================
function D = newD(cellKey, channels)
D = struct('cell', char(cellKey), 'channels', channels, 'tracks', struct([]), ...
           'spots', emptySpots(keysOf(channels)), 'masks', struct());
end

function D = addChannel(D, R)
K = keysOf(D.channels);
k = find(strcmp(K, R.key), 1);
assert(~isempty(k), 'dc_dataset:unknownChannel', ...
    'colour ''%s'' is not in the registry (%s)', R.key, strjoin(K, ', '));
ch = D.channels(k);

% The colour's tracks keep their own 0-based numbering; the offset carries them into the cell's.
off = 0;
if ~isempty(D.spots) && any(isfinite(D.spots.trackId))
    off = max(D.spots.trackId(isfinite(D.spots.trackId))) + 1;
end
T = struct('key', R.key, 'tracks', {R.tracks}, 'nFrames', R.nFrames, 'nTracks', R.nTracks, ...
           'idOffset', off);
if isempty(D.tracks), D.tracks = T; else, D.tracks(end+1) = T; end

n = numel(R.frame);
if n > 0
    f = double(R.frame(:));
    assert(all(f >= 0 & f < ch.nFrames), 'dc_dataset:frameRange', ...
        'colour ''%s'': a frame is outside 0..%d', R.key, ch.nFrames-1);
    local = double(R.trackId(:));
    assert(all(~isfinite(local) | (local >= 0 & local < max(R.nTracks,1))), 'dc_dataset:badLocalTrack', ...
        'colour ''%s'': a spot names a track index this colour does not have', R.key);
    % Intensity is optional: a caller that only has positions (a synthetic fixture, an import from
    % a tracker that did not measure) gets NaN rather than a missing column, so every consumer can
    % assume the column exists and check isfinite.
    new = table(repmat(categorical({R.key}, K), n, 1), f, ch.tp(f+1), ch.pages(f+1), ...
                double(R.x(:)), double(R.y(:)), double(R.q(:)), ...
                col(R,'iMean',n), col(R,'iMax',n), col(R,'iTot',n), ...
                local + off, local, double(R.spotId(:)), ...
                'VariableNames', {'ch','frame','tp','page','x','y','q', ...
                                  'iMean','iMax','iTot','trackId','trackLocal','spotId'});
    D.spots = [D.spots; new];
end
D = validate(D);
end

function D = validate(D)
K = keysOf(D.channels);
assert(numel(unique(K)) == numel(K), 'dc_dataset:dupKey', 'two colours share a key');
S = D.spots;
if isempty(S), return; end
assert(all(ismember(cellstr(S.ch), K)), 'dc_dataset:badFlag', ...
    'a spot carries a channel flag that is not in the registry');
for i = 1:numel(D.channels)
    ch = D.channels(i);
    m = S.ch == ch.key;
    if ~any(m), continue; end
    f = S.frame(m);
    assert(all(f >= 0 & f == round(f) & f < ch.nFrames), 'dc_dataset:badFrame', ...
        'colour ''%s'': frames must be 0-based integers below %d', ch.key, ch.nFrames);
    assert(isequal(S.tp(m), ch.tp(f+1)) && isequal(S.page(m), ch.pages(f+1)), ...
        'dc_dataset:indexDrift', ...
        ['colour ''%s'': a spot''s timepoint or page is not what the channel map says for its ' ...
         'frame. These are indices the acquisition supplied; a second copy that disagrees means ' ...
         'one of them was computed instead of read.'], ch.key);
    sid = S.spotId(m);
    assert(numel(unique(sid)) == numel(sid), 'dc_dataset:dupSpot', ...
        'colour ''%s'': spot ids repeat within the colour', ch.key);
    t = find(strcmp(keysOfTracks(D), ch.key), 1);
    if ~isempty(t)
        loc = S.trackLocal(m); gl = S.trackId(m); ok = isfinite(loc);
        assert(all(~ok | (loc(ok) >= 0 & loc(ok) < D.tracks(t).nTracks)), 'dc_dataset:badLocalTrack', ...
            ['colour ''%s'': a spot names a track index this colour does not have. The local index ' ...
             'is what finds a track in the cell the tracker returned, so it has to stay valid.'], ch.key);
        assert(all(~ok | gl(ok) == loc(ok) + D.tracks(t).idOffset), 'dc_dataset:trackIdDrift', ...
            ['colour ''%s'': a global track id is not its local index plus the recorded offset. ' ...
             'The two numbering schemes have to agree, or a selection and an export mean different ' ...
             'tracks.'], ch.key);
    end
end
tid = S.trackId; ok = isfinite(tid);
if any(ok)
    pairs = unique([double(S.ch(ok)), tid(ok)], 'rows');
    [~, ia] = unique(pairs(:,2));
    assert(numel(ia) == size(pairs,1), 'dc_dataset:trackSpansChannels', ...
        ['a track id appears in more than one colour. Tracking is per colour — a track that spans ' ...
         'them would be a link the tracker never made.']);
end
end

function v = col(R, f, n)
if isfield(R, f) && numel(R.(f)) == n, v = double(R.(f)(:)); else, v = nan(n,1); end
end

function K = keysOf(C), K = arrayfun(@(c) char(c.key), C(:)', 'uni', 0); end

function K = keysOfTracks(D)
if isempty(D.tracks), K = {}; return, end
K = arrayfun(@(t) char(t.key), D.tracks(:)', 'uni', 0);
end

function T = emptySpots(K)
T = table(categorical(cell(0,1), K), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
          zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
          zeros(0,1), zeros(0,1), zeros(0,1), ...
          'VariableNames', {'ch','frame','tp','page','x','y','q', ...
                            'iMean','iMax','iTot','trackId','trackLocal','spotId'});
end
