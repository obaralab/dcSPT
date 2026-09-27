function R = dc_import_spots(csvPath, ch, opts)
%DC_IMPORT_SPOTS  One colour's tracks read from a spots CSV instead of tracked here.
%
%   R = dc_import_spots(path, ch)
%   R = dc_import_spots(path, ch, struct('pxUm',0.0968,'alreadyUm',true))
%
% ch : the channel entry this file belongs to (from dc_channels), so the frames it names can be
%      checked against the ones the acquisition says that colour has.
%
% WHY IMPORT AT ALL. Tracking is the slow, fiddly, already-solved part. A project that has been
% through SPTinMatlab — or TrackMate, or anything that writes one row per localization — has tracks
% that were tuned and curated by hand, and re-deriving them here would be both slower and different.
% The co-motion analysis only needs positions, frames and track ids; where they came from is not its
% business.
%
% THE COLUMNS IT LOOKS FOR, case-insensitively, with the alternatives each accepts:
%
%   track   TRACK_ID | TRACKID | track | id          0-based or 1-based; blank/NaN = not tracked
%   frame   FRAME | frame | t | POSITION_T           0-based
%   x, y    X_um | X | POSITION_X                    microns unless 'alreadyUm' is false
%   iTot    TOTAL_INTENSITY | TOTAL_INTENSITY_CH1 | INTENSITY   optional; bleaching needs it
%   q       QUALITY                                  optional
%
% UNITS ARE THE ONE THING IT CANNOT GUESS. A column called X might be microns or pixels, and the
% difference is silent: every co-motion radius would be wrong by the pixel size and nothing would
% look broken. When the header names a unit (X_um) it is believed; otherwise 'alreadyUm' decides,
% and it defaults to TRUE with a warning, because a file written by this pipeline's own exporter is
% in microns. Pass it explicitly when importing from anything else.
%
% Frames are checked against the channel's own frame count, and a file naming frames the colour does
% not have is refused rather than silently truncated — that is the signature of the two colours'
% files having been swapped.

if nargin < 3 || ~isstruct(opts), opts = struct(); end
assert(isfile(csvPath), 'dc_import_spots:noFile', 'no file at %s', csvPath);
T = readtable(csvPath, 'VariableNamingRule','preserve');
V = string(T.Properties.VariableNames);

fTrack = pick(V, ["TRACK_ID","TRACKID","track_id","track","TrackID","id"]);
fFrame = pick(V, ["FRAME","frame","POSITION_T","t","Frame"]);
fX     = pick(V, ["X_um","x_um","X","x","POSITION_X"]);
fY     = pick(V, ["Y_um","y_um","Y","y","POSITION_Y"]);
assert(~isempty(fFrame) && ~isempty(fX) && ~isempty(fY), 'dc_import_spots:columns', ...
    ['%s needs a frame column and x/y columns. Found: %s'], csvPath, strjoin(V, ', '));

x = double(T.(fX)); y = double(T.(fY)); fr = double(T.(fFrame));
nameSaysUm = contains(lower(fX), '_um');
alreadyUm = getf(opts,'alreadyUm', true);
if ~nameSaysUm && ~isfield(opts,'alreadyUm')
    warning('dc_import_spots:units', ...
        ['%s has no unit in its x column name ("%s"), so microns is assumed. If it is in pixels, ' ...
         'every separation and every co-motion radius will be wrong by the pixel size and nothing ' ...
         'will look broken — pass alreadyUm=false.'], csvPath, fX);
end
if ~alreadyUm
    px = getf(opts,'pxUm',[]);
    assert(~isempty(px), 'dc_import_spots:noPixel', ...
        'alreadyUm is false, so a pixel size is needed to convert %s to microns', csvPath);
    x = x*px; y = y*px;
end

tid = nan(numel(fr),1);
if ~isempty(fTrack)
    tid = double(T.(fTrack));
    if iscell(tid), tid = cellfun(@(v) str2double(string(v)), T.(fTrack)); end
    ok = isfinite(tid);
    if any(ok) && min(tid(ok)) >= 1
        tid(ok) = tid(ok) - 1;         % 0-based inside this toolkit; a 1-based file is shifted once
    end
end

% frames the colour does not have are a swapped pair of files, not a rounding problem
assert(min(fr) >= 0, 'dc_import_spots:frameRange', '%s has negative frame numbers', csvPath);
if ~isempty(ch)
    assert(max(fr) < ch.nFrames, 'dc_import_spots:frameRange', ...
        ['%s names frame %d, but colour %s has only %d frames (0..%d). Two colours'' files swapped ' ...
         'looks exactly like this.'], csvPath, max(fr), ch.key, ch.nFrames, ch.nFrames-1);
end

nT = 0; okT = isfinite(tid);
if any(okT), nT = max(tid(okT)) + 1; end
R = struct('key', ch.key, 'label', ch.label, ...
    'frame', fr, 'x', x, 'y', y, ...
    'q',     col(T, pick(V, ["QUALITY","quality","q"]), numel(fr)), ...
    'iMean', col(T, pick(V, ["MEAN_INTENSITY","mean_intensity"]), numel(fr)), ...
    'iMax',  col(T, pick(V, ["MAX_INTENSITY","max_intensity"]), numel(fr)), ...
    'iTot',  col(T, pick(V, ["TOTAL_INTENSITY","total_intensity","INTENSITY","intensity"]), numel(fr)), ...
    'trackId', tid, 'spotId', (0:numel(fr)-1)', ...
    'tracks', {repmat({zeros(0,3)}, 1, nT)}, ...
    'nFrames', ch.nFrames, 'nTracks', nT, 'nDets', numel(fr), ...
    'source', csvPath);
end

% =================================================================================================
function f = pick(V, names)
f = '';
for n = names
    k = find(strcmpi(V, n), 1);
    if ~isempty(k), f = char(V(k)); return; end
end
end

function v = col(T, f, n)
if isempty(f), v = nan(n,1); return; end
v = double(T.(f));
if numel(v) ~= n, v = nan(n,1); end
end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
