function varargout = dc_project(mode, varargin)
%DC_PROJECT  The project on disk: what cells there are, how they were matched, and what was written.
%
%   P = dc_project('new',   folder)
%   P = dc_project('scan',  folder, matchOpts)      match the pairs and build the manifest
%   P = dc_project('load',  folder)
%       dc_project('save',  P)
%   p = dc_project('paths', P, cellKey)             where that cell's files go
%       dc_project('writeTracks', P, cellKey, D)    one spots CSV per colour
%   D = dc_project('readTracks',  P, cellKey, C)    and back again
%
% WHY A PROJECT AND NOT A SESSION. Everything this tool computed used to live in memory and die with
% the window: the pairs it matched, the thresholds that were tuned, which tracks were rejected by
% hand. That is fine for a demonstration and useless for work, because the expensive parts —
% tracking, and the human judgement in curation — are exactly the parts that were being thrown away.
%
% THE LAYOUT, which follows SPTinMatlab's so the two are navigable by the same habits:
%
%   <project>/
%     experiment_dc.mat          the manifest: cells, the matching rule, calibration, conditions
%     tracks/
%       <cell>_<colour>_spots.csv    one row per localization, per colour. The curatable artefact:
%                                    edit or filter it and read it back, the way the single-colour
%                                    pipeline's spots CSVs work.
%     analysis/
%       <cell>_dcspt.mat             the dataset and the co-motion result for that cell
%
% SPOTS CSVs ARE THE HAND-OFF, not the .mat. A CSV can be opened, sorted, filtered and understood
% without this toolkit, and a track rejected by deleting rows in a spreadsheet comes back in as a
% rejected track. The .mat is a cache of what those CSVs mean; delete it and it rebuilds.
%
% THE MATCHING RULE IS SAVED WITH THE PROJECT. Which two files are one cell is a decision, not a
% fact, and a project re-scanned next month with a different default would silently pair different
% files. P.match holds the tokens or pattern that were used, and 'scan' reuses them unless told
% otherwise.
%
% P — the manifest:
%   .folder .created .version
%   .match      the rule: .tokens or .pattern, as passed to dc_match
%   .cells      1xN struct: key, stacks{1x2}, chTokens, interleaved, condition, day, exclude,
%               notes, pxUm, dtS, pxSrc, dtSrc, status, nTracks, nPairs
%   .params     the shared settings, so a project reopens where it was left

switch lower(char(mode))
    case 'new',        varargout{1} = newP(varargin{:});
    case 'scan',       varargout{1} = scanP(varargin{:});
    case 'load',       varargout{1} = loadP(varargin{:});
    case 'save',       saveP(varargin{:});
    case 'paths',      varargout{1} = pathsFor(varargin{:});
    case 'writetracks',varargout{1} = writeTracks(varargin{:});
    case 'readtracks', varargout{1} = readTracks(varargin{:});
    otherwise, error('dc_project:mode','unknown mode ''%s''', mode);
end
end

% =================================================================================================
function P = newP(folder)
folder = char(folder);
P = struct('folder',folder, 'created',datetime('now'), 'version',1, ...
           'match', struct('tokens',{{'Ch1','Ch2'}}, 'pattern',''), ...
           'cells', emptyCells(), 'params', struct());
ensure(folder);
end

function P = scanP(folder, mopts)
% Match the pairs, and keep whatever a previous manifest already knew about the same cells. A
% re-scan after adding one movie must not discard the conditions and calibrations already entered.
folder = char(folder);
if nargin < 2 || ~isstruct(mopts), mopts = struct(); end
old = emptyCells(); P = newP(folder);
if isfile(manifestPath(folder))
    P = loadP(folder);
    old = P.cells;
    if ~isfield(mopts,'tokens') && ~isfield(mopts,'pattern'), mopts = P.match; end
end
P.match = struct('tokens', {getf(mopts,'tokens',{'Ch1','Ch2'})}, 'pattern', getf(mopts,'pattern',''));

[found, info] = dc_match(folder, P.match);
P.scanInfo = info;
C = emptyCells();
for k = 1:numel(found)
    c = blankCell();
    c.key = found(k).key; c.stacks = found(k).stacks;
    c.chTokens = found(k).chTokens; c.interleaved = found(k).interleaved;
    j = find(strcmp({old.key}, c.key), 1);
    if ~isempty(j)                      % keep what a human typed; refresh only what was matched
        keep = {'condition','day','exclude','notes','pxUm','dtS','pxSrc','dtSrc','status', ...
                'nTracks','nPairs','rejected'};
        for f = keep, if isfield(old(j),f{1}), c.(f{1}) = old(j).(f{1}); end, end
    end
    C(end+1) = c; %#ok<AGROW>
end
P.cells = C;
end

function P = loadP(folder)
p = manifestPath(char(folder));
assert(isfile(p), 'dc_project:noManifest', 'no experiment_dc.mat in %s', folder);
L = load(p);
assert(isfield(L,'P'), 'dc_project:badManifest', '%s does not hold a dcSPT manifest', p);
P = L.P;
P.folder = char(folder);                % a project that moved still knows where it is now
end

function saveP(P)
ensure(P.folder);
p = manifestPath(P.folder);
save(p, 'P', '-v7.3');
end

function q = pathsFor(P, key)
q = struct();
q.tracks   = fullfile(P.folder, 'tracks');
q.analysis = fullfile(P.folder, 'analysis');
q.spots    = @(ch) fullfile(q.tracks, sprintf('%s_%s_spots.csv', key, ch));
q.dataset  = fullfile(q.analysis, sprintf('%s_dcspt.mat', key));
end

function files = writeTracks(P, key, D)
% One CSV per colour, in the column shape the single-colour pipeline writes, so the same habits and
% the same spreadsheet work on both. TRACK_ID blank where a localization was never linked.
q = pathsFor(P, key); ensure(q.tracks);
files = {};
for ch = categories(D.spots.ch)'
    S = D.spots(D.spots.ch == ch{1}, :);
    if isempty(S), continue; end
    S = sortrows(S, {'trackId','tp'});
    T = table(S.trackId, S.spotId, S.frame, S.tp, S.page, S.x, S.y, S.q, ...
              S.iMean, S.iMax, S.iTot, ...
        'VariableNames', {'TRACK_ID','SPOT_ID','FRAME','TIMEPOINT','PAGE','X_um','Y_um', ...
                          'QUALITY','MEAN_INTENSITY','MAX_INTENSITY','TOTAL_INTENSITY'});
    f = q.spots(ch{1});
    writetable(T, f);
    files{end+1} = f; %#ok<AGROW>
end
end

function D = readTracks(P, key, C)
q = pathsFor(P, key);
D = dc_dataset('new', key, C);
for i = 1:numel(C)
    f = q.spots(C(i).key);
    assert(isfile(f), 'dc_project:noTracks', 'no spots file for %s at %s', C(i).key, f);
    R = dc_import_spots(f, C(i));
    D = dc_dataset('addChannel', D, R);
end
end

% =================================================================================================
function p = manifestPath(folder), p = fullfile(folder, 'experiment_dc.mat'); end
function ensure(d), if ~isfolder(d), mkdir(d); end, end

function C = emptyCells()
C = struct('key',{},'stacks',{},'chTokens',{},'interleaved',{},'condition',{},'day',{}, ...
           'exclude',{},'notes',{},'pxUm',{},'dtS',{},'pxSrc',{},'dtSrc',{},'status',{}, ...
           'nTracks',{},'nPairs',{},'rejected',{});
end

function c = blankCell()
c = struct('key','','stacks',{{'',''}},'chTokens',{{'',''}},'interleaved',false, ...
           'condition','','day','','exclude',false,'notes','', ...
           'pxUm',[],'dtS',[],'pxSrc','','dtSrc','','status','matched', ...
           'nTracks',[],'nPairs',[],'rejected',[]);
end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
