function dc_dataset_smoke()
%DC_DATASET_SMOKE  The data model: one flag per spot, one integer clock, and a merge driven by
%TIMEPOINT rather than by frame number.
%
% WHAT IS ASSERTED:
%   1. THE FLAG NAMES SOMETHING. Every spot's channel is a key in the registry, and a spot whose
%      flag is not is rejected rather than carried.
%   2. THE INDICES LIVE IN ONE PLACE. A spot's timepoint and page are EXACTLY what the channel map
%      says for its frame, and a copy that has drifted is caught. This is integer equality: there is
%      nothing approximate about an index, so a tolerance here would only hide a bug.
%   3. A TRACK BELONGS TO ONE COLOUR, and its id says which without being asked. Each colour is
%      tracked on its own and so numbers from zero — 'track 0' exists in both, for the most ordinary
%      reason there is — so the ids are made global on the way in and the colour's own index is kept
%      beside them. A track spanning both colours would be a link the tracker never made.
%   4. THE FOUR RELATIONS are told apart, INCLUDING the dangerous one: the same frame COUNT and the
%      same timepoints, but the two exposures are consecutive pages rather than one instant. That is
%      the case a panel stepping by frame index gets silently wrong.
%   5. THE MERGE IS BY TIMEPOINT. Asking for frame k of one colour returns the other colour's frame
%      at that TIMEPOINT, not its own frame k.
%   6. AN UNIMAGED COLOUR IS ABSENT, NOT STALE. Where one colour has no frame at that timepoint its
%      spots are left out and the caller is told — a stale position drawn beside a current one is how
%      a merged view invents colocalization.
%
% Synthetic; reads no dataset.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'core'), fullfile(fileparts(here),'drivers'), here);

dt = 0.0267118;

%% (1)(2)(3) the invariants -------------------------------------------------------------------------
% 'a' on the odd pages, 'b' on every fourth page — so 'b' shares only every other timepoint.
a = mkCh('a','fast', (1:2:20)', (1:10)',  2*dt);
b = mkCh('b','slow', (2:4:20)', (1:2:9)', 4*dt);
D = dc_dataset('new','cellA',[a b]);
% Both colours' trackers numbered from 0: 'a' has tracks 0 and 1, 'b' has tracks 0 and 1 too.
D = dc_dataset('addChannel', D, fakeR(a, 0:9, [0 0 0 0 0 1 1 1 1 1]));
D = dc_dataset('addChannel', D, fakeR(b, 0:4, [0 0 1 1 1]));
assert(height(D.spots) == 15, 'both colours'' spots should be in one list, got %d', height(D.spots));
assert(all(ismember(cellstr(unique(D.spots.ch)), {'a','b'})), 'and each flagged');
tpa = D.spots.tp(D.spots.ch=='a'); tpb = D.spots.tp(D.spots.ch=='b');
assert(isequal(tpa, (1:10)'), 'the fast colour is at every timepoint, got %s', mat2str(tpa'));
assert(isequal(tpb, (1:2:9)'), 'the slow colour at every other one, got %s', mat2str(tpb'));
assert(isequal(D.spots.page(D.spots.ch=='b'), (2:4:20)'), 'and each spot remembers its page');

% the ids collided at 0 and 1 in both colours; globally they must not
ga = unique(D.spots.trackId(D.spots.ch=='a')); gb = unique(D.spots.trackId(D.spots.ch=='b'));
assert(isempty(intersect(ga,gb)), ...
    ['both trackers numbered from zero, so the ids collided; globally they must be distinct ' ...
     '(%s vs %s)'], mat2str(ga'), mat2str(gb'));
assert(isequal(ga',[0 1]) && isequal(gb',[2 3]), ...
    'the second colour''s tracks continue the numbering, got %s and %s', mat2str(ga'), mat2str(gb'));
assert(isequal(unique(D.spots.trackLocal(D.spots.ch=='b'))',[0 1]), ...
    'while each colour keeps its own index, which is what finds a track in the tracker''s output');
assert(D.tracks(2).idOffset == 2, 'and the mapping between them is recorded, got %g', D.tracks(2).idOffset);

bad3 = D; bad3.spots.trackId(end) = 99;
err = ''; try, dc_dataset('validate', bad3); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_dataset:trackIdDrift'), ...
    'a global id that no longer matches its local index must be caught, got "%s"', err);
bad4 = D; bad4.spots.trackLocal(end) = 7; bad4.spots.trackId(end) = 9;
err = ''; try, dc_dataset('validate', bad4); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_dataset:badLocalTrack'), ...
    'and so must a local index naming a track the colour does not have, got "%s"', err);

bad = D; bad.spots.tp(3) = bad.spots.tp(3) + 1;          % an index that drifted
err = ''; try, dc_dataset('validate', bad); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_dataset:indexDrift'), ...
    'a timepoint that disagrees with the channel map must be caught, got "%s"', err);
bad = D; bad.spots.page(3) = bad.spots.page(3) + 1;
err = ''; try, dc_dataset('validate', bad); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_dataset:indexDrift'), 'and so must a page, got "%s"', err);

% What a WRONG offset would look like: the second colour renumbered from zero as well, so both
% colours' local and global ids agree but the ids themselves collide. This is the failure the global
% numbering exists to prevent, so it is worth having a backstop that catches it.
bad2 = D; bad2.tracks(2).idOffset = 0;
mb = D.spots.ch=='b'; bad2.spots.trackId(mb) = bad2.spots.trackLocal(mb);
err = ''; try, dc_dataset('validate', bad2); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_dataset:trackSpansChannels'), 'a track cannot span colours, got "%s"', err);

cch = mkCh('c','other',(1:4)',(1:4)',dt);
err = ''; try, dc_dataset('addChannel', D, fakeR(cch, 0:3)); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_dataset:unknownChannel'), 'a flag naming nothing must be refused, got "%s"', err);

%% (4) the four relations ---------------------------------------------------------------------------
% Interleaved: BOTH colours have a frame at every one of the 10 timepoints, and the same frame count
% — but their exposures are consecutive pages. 'matched' must say that rather than imply simultaneity.
A1 = dc_align(mkCh('a','a',(1:2:20)',(1:10)',2*dt), mkCh('b','b',(2:2:20)',(1:10)',2*dt));
assert(strcmp(A1.relation,'matched'), 'the same timepoints -> matched (got %s)', A1.relation);
assert(isequal(A1.nFrames,[10 10]), 'and the same frame count — which is what makes it dangerous');
assert(A1.pageGap == 1, 'the exposures are one page apart, got %g', A1.pageGap);
assert(contains(A1.text,'not simultaneous'), ...
    'a panel must be told the two exposures are consecutive, not one instant: "%s"', A1.text);

A2 = dc_align(mkCh('a','a',(1:2:20)',(1:10)',2*dt), mkCh('b','b',(2:4:20)',(1:2:9)',4*dt));
assert(strcmp(A2.relation,'subsampled'), 'half the timepoints -> subsampled (got %s)', A2.relation);
assert(contains(A2.text,'1 in 2'), 'and said how sparse: "%s"', A2.text);

A3 = dc_align(mkCh('a','a',(1:10)',(1:10)',dt), mkCh('b','b',(1:10)',(6:15)',dt));
assert(strcmp(A3.relation,'partial'), 'an overlap that contains neither -> partial (got %s)', A3.relation);
assert(A3.nShared == 5, 'sharing 5 timepoints, got %d', A3.nShared);

A4 = dc_align(mkCh('a','a',(1:5)',(1:5)',dt), mkCh('b','b',(6:10)',(20:24)',dt));
assert(strcmp(A4.relation,'disjoint'), 'no shared timepoint -> disjoint (got %s)', A4.relation);
assert(contains(A4.text,'Nothing can be'), 'and said plainly: "%s"', A4.text);

%% (5) the merge is by timepoint --------------------------------------------------------------------
% Frame 2 of 'a' is timepoint 3. 'b' HAS a frame there — its frame 1, not its frame 2.
[S5, i5] = dc_merge(D, struct('frame', 2, 'ref', 'a'));
assert(i5.tp == 3, 'a frame query resolves through the reference colour''s own map, got tp %g', i5.tp);
bRow = i5.perChannel(strcmp({i5.perChannel.key},'b'));
assert(bRow.present && bRow.frame == 1, ...
    ['the slow colour''s frame at timepoint 3 is 1, and a panel stepping by index would have drawn ' ...
     'its frame 2 — a different moment (got %g)'], bRow.frame);
assert(all(S5.frame(S5.ch=='a') == 2) && all(S5.frame(S5.ch=='b') == 1), ...
    'and the table holds each colour''s own frame at that timepoint');
assert(all(S5.tp == 3), 'every row of a merged view is the one timepoint');

%% (6) an unimaged colour is absent, not stale ------------------------------------------------------
% Timepoint 4: 'a' has a frame, 'b' does not (it is on the odd timepoints only).
[S6, i6] = dc_merge(D, struct('tp', 4));
assert(any(S6.ch=='a'), 'the colour that was imaged then still shows');
assert(~any(S6.ch=='b'), 'the colour that was not must contribute nothing');
bIn = i6.perChannel(strcmp({i6.perChannel.key},'b'));
assert(~bIn.present && bIn.n == 0 && isnan(bIn.frame), ...
    'and its row says it had no frame, rather than a stale one');
assert(isequal(i6.missing,{'b'}) && contains(i6.text,'no frame at this timepoint'), ...
    'the panel line must name it: "%s"', i6.text);

fprintf('data model: %d spots over 2 colours, one integer clock; relations %s / %s / %s / %s told apart\n', ...
    height(D.spots), A1.relation, A2.relation, A3.relation, A4.relation);
fprintf('merge: %s\n', i5.text);
fprintf('\nDC-DATASET SMOKE PASSED.\n');
end

% =================================================================================================
function c = mkCh(key, label, pages, tp, dt)
c = struct('key',key, 'label',label, 'pages',pages(:), 'tp',tp(:), 'dt_s',dt, 'nFrames',numel(pages));
end

function R = fakeR(ch, frames, tid)
f = double(frames(:));
if nargin < 3, tid = nan(numel(f),1); end
tid = double(tid(:));
nT = 0; if any(isfinite(tid)), nT = max(tid(isfinite(tid))) + 1; end
R = struct('key',ch.key, 'frame',f, 'x',5+0.01*f, 'y',5-0.01*f, 'q',100+0*f, ...
    'trackId',tid, 'spotId',(0:numel(f)-1)', 'tracks',{repmat({zeros(0,3)},1,nT)}, ...
    'nFrames',numel(f), 'nTracks',nT, 'nDets',numel(f));
end
