function D = dc_dataset(mode, varargin)
%DC_DATASET  One cell, both colours: build it, and check it is still true.
%
%   D = dc_dataset('new', cellKey, channels)
%   D = dc_dataset('addChannel', D, R)        % R from dc_process_cell, one colour
%   D = dc_dataset('validate', D)             % errors on anything that cannot be true
%
% See docs/DATA_MODEL.md. The short version: the per-channel arrays are canonical and stay on their
% own frame axis, because the two colours do not share one; the flat spot list is derived from them
% and is the only thing that can be merged, because in it time is a column rather than an index.
%
% THE INVARIANTS, all checked by 'validate':
%   1. every spot's channel is a key in the registry — a flag that names nothing is worse than none
%   2. t_s == t0_s + frame*dt_s, for that spot's OWN channel. This is the one that catches a second
%      copy of dt drifting from the first, which is how the previous attempt lost 1.1% of its clock
%   3. frames are 0-based within a channel and never compared across channels
%   4. a track id belongs to exactly one channel: tracking is per colour, so a track cannot span them
%   5. spot ids are unique within a channel
%
% The structure is deliberately dumb data — no handles, no closures — so it saves, loads, and can be
% compared field by field in a test.

switch lower(char(mode))
    case 'new'
        D = newD(varargin{:});
    case 'addchannel'
        D = addChannel(varargin{:});
    case 'validate'
        D = validate(varargin{:});
    otherwise
        error('dc_dataset:mode','unknown mode ''%s''', mode);
end
end

% =================================================================================================
function D = newD(cellKey, channels)
D = struct('cell', char(cellKey), 'channels', channels, 'tracks', struct([]), ...
           'spots', emptySpots(keysOf(channels)), 'masks', struct('er','','mito','','dtPage_s',NaN));
end

function D = addChannel(D, R)
% Fold one colour's processed result into the dataset: its tracks stay on their own axis, its spots
% join the shared list with the flag and the cell clock.
k = find(strcmp(keysOf(D.channels), R.key), 1);
assert(~isempty(k), 'dc_dataset:unknownChannel', ...
    'colour ''%s'' is not in the registry (%s)', R.key, strjoin(keysOf(D.channels), ', '));
ch = D.channels(k);
assert(abs(ch.dt_s - R.dt_s) < 1e-12 && abs(ch.t0_s - R.t0_s) < 1e-12, ...
    'dc_dataset:clockMismatch', ...
    ['colour ''%s'' was processed at dt %.6g / t0 %.6g but the registry says %.6g / %.6g. One of ' ...
     'them is a second copy of the clock, which is exactly what must not exist.'], ...
    R.key, R.dt_s, R.t0_s, ch.dt_s, ch.t0_s);

T = struct('key', R.key, 'tracks', {R.tracks}, 'nFrames', R.nFrames, 'nTracks', R.nTracks);
if isempty(D.tracks), D.tracks = T; else, D.tracks(end+1) = T; end

n = numel(R.frame);
if n > 0
    new = table(repmat(categorical({R.key}, keysOf(D.channels)), n, 1), ...
                double(R.frame(:)), ch.t0_s + double(R.frame(:))*ch.dt_s, ...
                double(R.x(:)), double(R.y(:)), double(R.q(:)), ...
                double(R.trackId(:)), double(R.spotId(:)), ...
                'VariableNames', {'ch','frame','t_s','x','y','q','trackId','spotId'});
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
    want = ch.t0_s + S.frame(m)*ch.dt_s;
    err = max(abs(S.t_s(m) - want));
    assert(err < 1e-9, 'dc_dataset:clockDrift', ...
        ['colour ''%s'': a spot''s time is %.3g s from t0 + frame*dt. The clock is in the registry ' ...
         'and nowhere else; a time that disagrees with it means a second copy has drifted.'], ch.key, err);
    assert(all(S.frame(m) >= 0) && all(S.frame(m) == round(S.frame(m))), ...
        'dc_dataset:badFrame', 'colour ''%s'': frames must be 0-based integers', ch.key);
    sid = S.spotId(m);
    assert(numel(unique(sid)) == numel(sid), 'dc_dataset:dupSpot', ...
        'colour ''%s'': spot ids repeat within the colour', ch.key);
end
% a track belongs to one colour only
tid = S.trackId; ok = isfinite(tid);
if any(ok)
    pairs = unique([double(S.ch(ok)), tid(ok)], 'rows');
    [~, ia] = unique(pairs(:,2));
    assert(numel(ia) == size(pairs,1), 'dc_dataset:trackSpansChannels', ...
        ['a track id appears in more than one colour. Tracking is per colour — a track that spans ' ...
         'them would be a link the tracker never made.']);
end
end

function K = keysOf(C)
K = arrayfun(@(c) char(c.key), C(:)', 'uni', 0);
end

function T = emptySpots(K)
T = table(categorical(cell(0,1), K), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
          zeros(0,1), zeros(0,1), zeros(0,1), ...
          'VariableNames', {'ch','frame','t_s','x','y','q','trackId','spotId'});
end
