function dc_project_smoke()
%DC_PROJECT_SMOKE  A project on disk: matched pairs, a manifest that survives, curatable tracks.
%
% WHAT IS ASSERTED:
%   1. PAIRS ARE MATCHED BY A RULE THAT IS VISIBLE. Two colours of one cell are found by their
%      token, and a stack holding both colours is one cell on its own without any token at all.
%   2. A KEY WITH THE WRONG NUMBER OF FILES IS REPORTED, NOT PAIRED WITH WHATEVER IS NEAREST. Two
%      colours of DIFFERENT cells would produce tracks, separations and a co-motion number that all
%      look entirely normal, so a bad pairing has to be refused loudly or not at all.
%   3. THE MANIFEST SURVIVES A ROUND TRIP, and the matching rule goes with it — which two files are
%      one cell is a decision, and a re-scan next month with a different default would silently pair
%      different files.
%   4. A RE-SCAN KEEPS WHAT A HUMAN TYPED. Adding a movie must not discard the conditions and
%      calibrations already entered against the other cells.
%   5. TRACKS ROUND-TRIP THROUGH CSV, so the artefact a person curates is the artefact the tool
%      reads back. Deleting rows in a spreadsheet is a valid way to reject a track.
%
% Synthetic; writes only under tempdir.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root,'core'), fullfile(root,'drivers'), here);

proj = fullfile(tempdir, sprintf('dc_proj_%d', feature('getpid')));
if isfolder(proj), rmdir(proj,'s'); end
mkdir(proj);
cleanup = onCleanup(@() rmdir(proj,'s'));

nT = 24; px = 0.1;
% two cells as separate colour files, one cell as a single interleaved stack, one orphan
writePair(proj, 'cellA', nT, px);
writePair(proj, 'cellB', nT, px);
writeDual(proj, 'cellD', nT, px);
writeOne (proj, 'cellC_Ch1', nT, px);        % no partner: must be reported, not paired

%% (1)(2) matching ---------------------------------------------------------------------------------
[found, info] = dc_match(proj, struct('tokens',{{'Ch1','Ch2'}}));
keys = sort(string({found.key}));
assert(numel(found) == 3, 'cellA, cellB and the interleaved cellD should pair; got %d (%s)', ...
    numel(found), strjoin(keys, ', '));
dual = found([found.interleaved]);
assert(numel(dual) == 1 && contains(dual.key,'cellD'), ...
    'a stack whose labels name two channels is one cell with no token needed');
assert(~isempty(info.unpaired), 'the orphan must be reported');
orph = [info.unpaired{:}];
assert(any(contains(string([orph.files]), 'cellC')), ...
    'and it must be cellC that is named: %s', strjoin(string([orph.files]), ', '));
fprintf('(1) %s\n', info.text);
fprintf('(2) unpaired reported: %s\n', strjoin(string([orph.files]), ', '));

%% (3) the manifest round-trips ---------------------------------------------------------------------
P = dc_project('scan', proj, struct('tokens',{{'Ch1','Ch2'}}));
assert(numel(P.cells) == 3, 'the manifest should hold the three matched cells');
P.cells(1).condition = 'Baseline';
P.cells(1).pxUm = px; P.cells(1).dtS = 0.012; P.cells(1).notes = 'typed by hand';
dc_project('save', P);
assert(isfile(fullfile(proj,'experiment_dc.mat')), 'experiment_dc.mat should have been written');

P2 = dc_project('load', proj);
assert(numel(P2.cells) == numel(P.cells), 'the cells should come back');
assert(strcmp(P2.cells(1).condition,'Baseline') && strcmp(P2.cells(1).notes,'typed by hand'), ...
    'and what was typed with them');
assert(isequal(P2.match.tokens, {'Ch1','Ch2'}), ...
    'the matching rule must be saved with the project, or a re-scan could pair different files');
fprintf('(3) manifest round-trips with the rule %s\n', strjoin(P2.match.tokens, '/'));

%% (4) a re-scan keeps what was typed ----------------------------------------------------------------
writePair(proj, 'cellE', nT, px);             % a movie arrives later
P3 = dc_project('scan', proj);                % no options: the saved rule is reused
assert(numel(P3.cells) == 4, 'the new cell should be picked up, got %d', numel(P3.cells));
j = find(strcmp({P3.cells.key}, P2.cells(1).key), 1);
assert(~isempty(j) && strcmp(P3.cells(j).condition,'Baseline'), ...
    'and the condition already entered must survive the re-scan');
assert(isequal(P3.match.tokens, {'Ch1','Ch2'}), 'with the rule carried forward');
fprintf('(4) re-scan: %d cells, conditions kept\n', numel(P3.cells));

%% (5) tracks round-trip through CSV -------------------------------------------------------------
C = dc_channels('manual', [ ...
    struct('key','c1','label','c1','pages',(1:nT)','tp',(1:nT)','dt_s',0.012), ...
    struct('key','c2','label','c2','pages',(1:nT)','tp',(1:nT)','dt_s',0.012)]);
D = fakeDataset(C, nT);
nTr0 = numel(unique(D.spots.trackId(isfinite(D.spots.trackId))));
files = dc_project('writeTracks', P3, 'cellA', D);
assert(numel(files) == 2 && all(cellfun(@isfile, files)), 'one spots CSV per colour');

D2 = dc_project('readTracks', P3, 'cellA', C);
nTr1 = numel(unique(D2.spots.trackId(isfinite(D2.spots.trackId))));
assert(nTr1 == nTr0, 'the track count must survive the round trip (%d -> %d)', nTr0, nTr1);
assert(max(abs(sort(D2.spots.x) - sort(D.spots.x))) < 1e-4, 'and the positions');

% curation by deleting rows in the CSV is a valid way to reject a track
T = readtable(files{1}, 'VariableNamingRule','preserve');
victim = T.TRACK_ID(find(isfinite(T.TRACK_ID),1));
writetable(T(T.TRACK_ID ~= victim, :), files{1});
D3 = dc_project('readTracks', P3, 'cellA', C);
n1 = numel(unique(D3.spots.trackId(D3.spots.ch=='c1' & isfinite(D3.spots.trackId))));
n0 = numel(unique(D2.spots.trackId(D2.spots.ch=='c1' & isfinite(D2.spots.trackId))));
assert(n1 == n0 - 1, ...
    ['deleting a track''s rows in the CSV must remove exactly that track (%d -> %d) — the file is ' ...
     'the curatable artefact, not a dump of one'], n0, n1);
fprintf('(5) %d tracks written and read back; a track deleted in the CSV stays deleted\n', nTr0);

fprintf('\nDC-PROJECT SMOKE PASSED.\n');
end

% =================================================================================================
function writePair(folder, key, nT, px)
for c = 1:2
    lab = arrayfun(@(t) sprintf('c:%d/2 t:%d/%d - %s #1', c, t, nT, key), 1:nT, 'uni', 0);
    ijStack(fullfile(folder, sprintf('%s_Ch%d.tif', key, c)), lab, nT, px);
end
end

function writeOne(folder, name, nT, px)
lab = arrayfun(@(t) sprintf('c:1/2 t:%d/%d - %s #1', t, nT, name), 1:nT, 'uni', 0);
ijStack(fullfile(folder, [name '.tif']), lab, nT, px);
end

function writeDual(folder, key, nT, px)
% one file, both colours interleaved: the case no token rule can see
lab = {};
for t = 1:nT
    lab{end+1} = sprintf('c:1/2 t:%d/%d - %s #1', t, nT, key); %#ok<AGROW>
    lab{end+1} = sprintf('c:2/2 t:%d/%d - %s #1', t, nT, key); %#ok<AGROW>
end
ijStack(fullfile(folder, [key '_dual.tif']), lab, 2*nT, px);
end

function ijStack(path, labels, n, px)
pages = arrayfun(@(k) uint16(100 + 20*rand(16)), 1:n, 'uni', 0);
dcWriteIJ(path, labels, pages, px);
end

function D = fakeDataset(C, nT)
rng(5);
D = dc_dataset('new','cellA',C);
for i = 1:2
    nTk = 4; fr = []; x = []; y = []; tid = [];
    for k = 1:nTk
        f = (0:nT-1)';
        fr = [fr; f]; %#ok<AGROW>
        x = [x; 2 + k + cumsum(0.02*randn(nT,1))]; %#ok<AGROW>
        y = [y; 3 + k + cumsum(0.02*randn(nT,1))]; %#ok<AGROW>
        tid = [tid; repmat(k-1, nT, 1)]; %#ok<AGROW>
    end
    R = struct('key',C(i).key, 'frame',fr, 'x',x, 'y',y, 'q',100+0*fr, ...
        'iMean',50+0*fr, 'iMax',80+0*fr, 'iTot',1000+0*fr, ...
        'trackId',tid, 'spotId',(0:numel(fr)-1)', 'tracks',{repmat({zeros(0,3)},1,nTk)}, ...
        'nFrames',nT, 'nTracks',nTk, 'nDets',numel(fr));
    D = dc_dataset('addChannel', D, R);
end
end
