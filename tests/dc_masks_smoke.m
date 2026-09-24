function dc_masks_smoke()
%DC_MASKS_SMOKE  One organelle mask, shared by both colours, indexed by time — and the assumption
%that makes sharing safe, checked rather than assumed.
%
% WHAT IS ASSERTED:
%   1. THE PAGE IS A FUNCTION OF TIME, not of frame number. Both colours at the same instant get the
%      SAME page, even though that instant is a different frame number in each — which is the whole
%      reason the lookup is by time.
%   2. THE MASK OUTLASTING THE MOVIE IS NORMAL: a time past the last page clamps rather than failing.
%   3. THE SHARING ASSUMPTION IS CHECKED. When the organelle is imaged far more slowly than the
%      particles — the case here — sharing is declared safe and says by how much. When it is not,
%      that is said plainly instead of being assumed.
%   4. A MASK WITH NO DECLARED RATE IS REFUSED, because deriving one from page counts needs an exact
%      integer ratio that real acquisitions do not give.
%
% Synthetic; reads no dataset.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'core'), fullfile(fileparts(here),'drivers'), here);
root = fullfile(tempdir, sprintf('dc_masks_%d', feature('getpid')));
if isfolder(root), rmdir(root,'s'); end
mkdir(root);
cleanup = onCleanup(@() rmdir(root,'s'));

dt = 0.0267118;                       % the fast colour
mk = @(key,dtv,t0) struct('key',key,'label',key,'stride',2,'offset',0,'dt_s',dtv,'t0_s',t0,'file','');
D = dc_dataset('new','cellA',[mk('a',dt,0), mk('b',2*dt,dt)]);

% a mito stack that changes every 2.7 s — 100 frames of the fast colour
S = 32; nPg = 6; maskDt = 100*dt;
p = fullfile(root,'mito.tif');
for k = 1:nPg
    m = false(S); m(6+k:20+k, 8:24) = true;       % it moves a little, page to page
    if k == 1, imwrite(uint8(m)*255, p); else, imwrite(uint8(m)*255, p, 'WriteMode','append'); end
end
D = dc_masks('set', D, struct('mito', p, 'dt_s', maskDt, 't0_s', 0));

%% (1) both colours at one instant get one page -----------------------------------------------------
t = 250*dt;                                        % some instant mid-movie
fa = round((t - D.channels(1).t0_s)/D.channels(1).dt_s);
fb = round((t - D.channels(2).t0_s)/D.channels(2).dt_s);
assert(fa ~= fb, 'the fixture should put this instant at different FRAME numbers in the two colours (%d, %d)', fa, fb);
[pgA, tPg] = dc_masks('page', D, D.channels(1).t0_s + fa*D.channels(1).dt_s);
[pgB, ~]   = dc_masks('page', D, D.channels(2).t0_s + fb*D.channels(2).dt_s);
assert(pgA == pgB, ...
    ['the same instant must give the same mask page in both colours (%d vs %d) — indexing by frame ' ...
     'number would have given two'], pgA, pgB);
assert(abs(tPg - t) <= maskDt/2 + 1e-9, 'and the page returned should be the one nearest in time');
m1 = dc_masks('at', D, 'mito', t);
assert(islogical(m1) && any(m1(:)), 'the page should come back as a mask');

%% (2) past the end of the mask stack ----------------------------------------------------------------
mLate = dc_masks('at', D, 'mito', 10000*dt);
assert(islogical(mLate) && any(mLate(:)), ...
    'a mask stack shorter than the movie is normal and must clamp, not fail');

%% (3) the assumption, checked ------------------------------------------------------------------------
C = dc_masks('check', D);
assert(C.ok, 'a mask 100x slower than the frames should be safe to share: "%s"', C.text);
assert(contains(C.text,'has not moved between them'), 'and should say why: "%s"', C.text);
% the ratio is against the COARSEST colour (2*dt here), which is the one most at risk of being
% given a stale mask: 100*dt / 2*dt = 50
assert(abs(C.ratio - 50) < 1e-9, 'the ratio should be against the slowest colour (%.2f, wanted 50)', C.ratio);
% now imaged nearly as fast as the particles: the assumption no longer holds
D2 = dc_masks('set', D, struct('mito', p, 'dt_s', 2*dt, 't0_s', 0));
C2 = dc_masks('check', D2);
assert(~C2.ok, 'a mask imaged at nearly the frame rate must NOT be declared safe to share');
assert(contains(C2.text,'no longer safe'), 'and must say so plainly: "%s"', C2.text);

%% (4) no declared rate -------------------------------------------------------------------------------
err = ''; try, dc_masks('set', D, struct('mito', p)); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_masks:noRate'), 'a mask with no page interval must be refused, got "%s"', err);

fprintf('masks: one page (%d) for both colours at t=%.3f s though their frame numbers there are %d and %d\n', ...
    pgA, t, fa, fb);
fprintf('check: %s\n', C.text);
fprintf('\nDC-MASKS SMOKE PASSED.\n');
end
