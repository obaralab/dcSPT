function dc_dataset_smoke()
%DC_DATASET_SMOKE  The data model: one flag per spot, one clock per colour, and a merge driven by
%time rather than by frame number.
%
% WHAT IS ASSERTED:
%   1. THE FLAG NAMES SOMETHING. Every spot's channel is a key in the registry, and a spot whose
%      flag is not is rejected rather than carried.
%   2. THE CLOCK LIVES IN ONE PLACE. t_s == t0 + frame*dt for that spot's OWN colour, and a dt that
%      has drifted from the registry is caught — the failure the first attempt actually had.
%   3. A TRACK BELONGS TO ONE COLOUR. Tracking is per colour, so a track spanning both would be a
%      link the tracker never made.
%   4. THE FOUR CLOCK RELATIONS are told apart, INCLUDING the dangerous one: same frame COUNT and
%      same interval, but never simultaneous. That is the case a panel stepping by frame index gets
%      silently wrong.
%   5. THE MERGE IS BY TIME. Asking for frame k of one colour returns the other colour's frame at
%      that INSTANT, not its frame k.
%   6. AN UNIMAGED COLOUR IS ABSENT, NOT STALE. Where one colour has no frame near the instant, its
%      spots are left out and the caller is told — a stale position drawn beside a current one is
%      how a merged view invents colocalization.
%
% Synthetic; reads no dataset.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'core'), fullfile(fileparts(here),'drivers'), here);

dt = 0.0267118;
mk = @(key,label,dtv,t0) struct('key',key,'label',label,'stride',2,'offset',0, ...
        'dt_s',dtv,'t0_s',t0,'file','');

%% (1)(2)(3) the invariants -------------------------------------------------------------------------
C = [mk('a','fast',dt,0), mk('b','slow',2*dt,dt)];
D = dc_dataset('new','cellA',C);
D = dc_dataset('addChannel', D, fakeR('a', 0:9, dt, 0));
D = dc_dataset('addChannel', D, fakeR('b', 0:4, 2*dt, dt));
assert(height(D.spots) == 15, 'both colours'' spots should be in one list, got %d', height(D.spots));
assert(all(ismember(cellstr(unique(D.spots.ch)), {'a','b'})), 'and each flagged');
ta = D.spots.t_s(D.spots.ch=='a'); tb = D.spots.t_s(D.spots.ch=='b');
assert(abs(ta(1)) < 1e-12 && abs(tb(1) - dt) < 1e-12, 'each colour on its own clock');
assert(abs(tb(2) - tb(1) - 2*dt) < 1e-12, 'the slow colour''s spots are 2 dt apart');

bad = D; bad.spots.t_s(3) = bad.spots.t_s(3) + 0.001;     % a clock that drifted
err = ''; try, dc_dataset('validate', bad); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_dataset:clockDrift'), ...
    'a time that disagrees with the registry must be caught, got "%s"', err);

bad2 = D; bad2.spots.trackId(D.spots.ch=='a') = 1; bad2.spots.trackId(D.spots.ch=='b') = 1;
err = ''; try, dc_dataset('validate', bad2); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_dataset:trackSpansChannels'), 'a track cannot span colours, got "%s"', err);

err = ''; try, dc_dataset('addChannel', D, fakeR('c', 0:3, dt, 0)); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_dataset:unknownChannel'), 'a flag naming nothing must be refused, got "%s"', err);
Rw = fakeR('a', 0:3, dt*1.5, 0);                           % processed on a different clock
err = ''; try, dc_dataset('addChannel', D, Rw); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_dataset:clockMismatch'), 'a second copy of the clock must be refused, got "%s"', err);

%% (4) the four relations ---------------------------------------------------------------------------
A1 = dc_align(mkN('a',dt,0,10),      mkN('b',dt,0,10));
assert(strcmp(A1.relation,'simultaneous') && A1.sameCount, 'same clock -> simultaneous');
A2 = dc_align(mkN('a',dt,0,10),      mkN('b',dt,dt/2,10));
assert(strcmp(A2.relation,'interleaved'), 'same dt, offset t0 -> interleaved (got %s)', A2.relation);
assert(A2.sameCount, 'and CRUCIALLY the same frame count — which is what makes it dangerous');
assert(contains(A2.text,'never imaged at the same instant'), ...
    'the panel must be told they are not simultaneous despite the matching counts: "%s"', A2.text);
A3 = dc_align(mkN('a',dt,0,10),      mkN('b',2*dt,dt,5));
assert(strcmp(A3.relation,'subsampled'), 'a 2x slower colour -> subsampled (got %s)', A3.relation);
assert(~A3.sameCount && contains(A3.text,'more slowly'), 'and said so: "%s"', A3.text);
A4 = dc_align(mkN('a',dt,0,10),      mkN('b',dt*1.37,0,7));
assert(strcmp(A4.relation,'unequal'), 'unrelated rates -> unequal (got %s)', A4.relation);
assert(contains(A4.text,'Frame numbers mean nothing'), 'and said plainly: "%s"', A4.text);
% a rounded dt (0.054 for 0.0534237, as the first attempt carried) must still read as subsampled
A5 = dc_align(mkN('a',dt,0,10), mkN('b',0.054,dt,5));
assert(strcmp(A5.relation,'subsampled'), 'a 1%% rounded dt should not become "unequal" (got %s)', A5.relation);

%% (5) the merge is by time -------------------------------------------------------------------------
[S3, i3] = dc_merge(D, struct('frame', 2, 'ref', 'a'));       % t = 2*dt; b's nearest frame is 0 or 1
assert(abs(i3.t_s - 2*dt) < 1e-12, 'a frame query resolves through the reference colour''s clock');
bRow = i3.perChannel(strcmp({i3.perChannel.key},'b'));
assert(bRow.frame ~= 2, ...
    ['the slow colour''s frame at this instant is %d, and a panel stepping by index would have ' ...
     'drawn its frame 2 — a different moment'], bRow.frame);
assert(abs(bRow.t_s - (dt + bRow.frame*2*dt)) < 1e-12, 'and it is that colour''s own clock that placed it');
assert(all(S3.frame(S3.ch=='a') == 2), 'the reference colour shows the frame asked for');

%% (6) an unimaged colour is absent, not stale --------------------------------------------------------
% A colour that stops early: ask for an instant past its last frame.
D2 = dc_dataset('new','cellB',[mk('a','fast',dt,0), mk('b','slow',2*dt,dt)]);
D2 = dc_dataset('addChannel', D2, fakeR('a', 0:19, dt, 0));
D2 = dc_dataset('addChannel', D2, fakeR('b', 0:1,  2*dt, dt));   % b stops after 2 frames
[S6, i6] = dc_merge(D2, struct('t_s', 15*dt));
assert(~any(S6.ch=='b'), 'the colour that was not imaged then must contribute nothing');
assert(any(S6.ch=='a'), 'while the one that was still shows');
bIn = i6.perChannel(strcmp({i6.perChannel.key},'b'));
assert(bIn.n == 0, 'and its row says it had nothing, rather than a stale position');
assert(contains(i6.text,'spots') , 'the panel line should describe what it drew: "%s"', i6.text);

fprintf('data model: %d spots over 2 colours, one clock each; relations %s / %s / %s / %s told apart\n', ...
    height(D.spots), A1.relation, A2.relation, A3.relation, A4.relation);
fprintf('merge: %s\n', i3.text);
fprintf('\nDC-DATASET SMOKE PASSED.\n');
end

% =================================================================================================
function R = fakeR(key, frames, dt, t0)
f = double(frames(:));
R = struct('key',key, 'frame',f, 'x',5+0.01*f, 'y',5-0.01*f, 'q',100+0*f, ...
    'trackId',nan(numel(f),1), 'spotId',(0:numel(f)-1)', 'tracks',{{}}, ...
    'dt_s',dt, 't0_s',t0, 'nFrames',numel(f), 'nTracks',0, 'nDets',numel(f));
end

function c = mkN(key, dt, t0, n)
c = struct('key',key,'label',key,'stride',1,'offset',0,'dt_s',dt,'t0_s',t0,'file','','nFrames',n);
end
