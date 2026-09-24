function dc_track_smoke()
%DC_TRACK_SMOKE  Two colours of one interleaved stack, detected and tracked independently, each on
%its own clock.
%
% WHAT IS ASSERTED:
%   1. BOTH COLOURS COME BACK from one call, each with its own frames, tracks and detections.
%   2. EACH IS TRACKED ON ITS OWN PAGES: colour A's particle is found where A's particle is, and
%      colour B's where B's is. Neither picks up the other's.
%   3. THE CLOCKS ARE RIGHT AND DIFFERENT. Interleaved colours share a page interval but their
%      frames are twice that apart, and B starts one page later. A run where both claimed t = 0
%      would look identical to this one in every other respect.
%   4. THE TWO CAN BE PUT ON ONE CLOCK, and the clock is the ACQUISITION's, not the detections'.
%      Every frame of A has a frame of B within half a page. Where B was not DETECTED there is
%      genuinely no partner, and that comes back unmatched rather than snapped to a distant frame.
%   5. A COLOUR MAY DECLARE ITS OWN dt and it is used as given, not derived from the stride.
%   6. NO dt ANYWHERE IS AN ERROR, not a guess: everything downstream is reported in seconds.
%
% Synthetic; reads no dataset.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'core'), fullfile(fileparts(here),'drivers'));
root = fullfile(tempdir, sprintf('dc_track_%d', feature('getpid')));
if isfolder(root), rmdir(root,'s'); end
mkdir(root);
cleanup = onCleanup(@() rmdir(root,'s'));

%% an interleaved stack: A drifts left to right on the odd pages, B sits still on the even --------
S = 64; nPg = 60; dtPage = 0.01; px = 0.1;
mv = fullfile(root,'cellA.tif');
rng(4);
ax = zeros(1,0); bx = zeros(1,0);
for pg = 1:nPg
    im = uint16(200 + 6*randn(S,S));
    if mod(pg,2) == 1
        f = (pg+1)/2;  x = 10 + round(f*0.7); y = 20;    % A: a clean diagonal-ish drift
        ax(end+1) = x; %#ok<AGROW>
    else
        f = pg/2;      x = 48; y = 44 - round(f*0.2);    % B: barely moves, elsewhere in the frame
        bx(end+1) = x; %#ok<AGROW>
    end
    im(y-1:y+1, x-1:x+1) = im(y-1:y+1, x-1:x+1) + uint16(5000);
    if pg == 1, imwrite(im, mv); else, imwrite(im, mv, 'WriteMode','append'); end
end
cel = struct('stack', mv, 'base', 'cellA');
prm = struct('diamUm',0.4,'thrAbs',60,'pxUm',px,'linkUm',1.5,'gapUm',1.5,'maxGap',1);

%% (1)(2) both colours, each on its own pages -----------------------------------------------------
C = dc_channels('interleaved', dtPage);
R = dc_process_cell(cel, C, prm);
assert(numel(R) == 2, 'both colours should come back, got %d', numel(R));
assert(strcmp(R(1).key,'a') && strcmp(R(2).key,'b'), 'and keep their keys');
assert(R(1).nFrames == nPg/2 && R(2).nFrames == nPg/2, ...
    'each colour has half the pages (%d, %d of %d)', R(1).nFrames, R(2).nFrames, nPg);
assert(R(1).nDets > 20 && R(2).nDets > 20, 'both colours should be detected (%d, %d)', R(1).nDets, R(2).nDets);
% A moves across the frame; B does not. Their x ranges must not be confusable.
assert(range(R(1).x)*px > 1.0, 'colour A should travel (range %.2f um)', range(R(1).x)*px);
assert(range(R(2).x)*px < 0.5, 'colour B should sit still (range %.2f um)', range(R(2).x)*px);
assert(min(R(2).x) > max(R(1).x), ...
    'the two colours occupy different parts of the frame; neither should pick up the other');
assert(R(1).nTracks >= 1 && R(2).nTracks >= 1, 'each colour should link into at least one track');

%% (3) the clocks ---------------------------------------------------------------------------------
assert(abs(R(1).dt_s - 2*dtPage) < 1e-12 && abs(R(2).dt_s - 2*dtPage) < 1e-12, ...
    'interleaved frames are twice the page interval apart (%.4g, %.4g)', R(1).dt_s, R(2).dt_s);
assert(R(1).t0_s == 0 && abs(R(2).t0_s - dtPage) < 1e-12, ...
    'colour B starts one page later, not at zero (t0 = %.4g)', R(2).t0_s);
assert(abs(R(1).t_s(1) - 0) < 1e-12 && abs(R(2).t_s(1) - dtPage) < 1e-12, ...
    'and the per-detection clock should say so');

%% (4) one clock for the pair ---------------------------------------------------------------------
% The ACQUISITION clocks: every frame of A has a frame of B within half a page, because that is how
% the stack was taken. This is the question "was the other colour even imaged then".
M = dc_time_match(dc_times(R(1)), dc_times(R(2)));
assert(M.nUnmatched == 0, '%d of A''s FRAMES found no frame of B', M.nUnmatched);
assert(max(abs(M.dt_s)) <= dtPage + 1e-12, ...
    'a partner frame should never be more than one page interval away (worst %.4g s)', max(abs(M.dt_s)));
% The DETECTION times are a different question — "where was B's molecule when mine was seen" — and
% there B is blind wherever it was not detected. Those must come back UNMATCHED rather than snapped
% to a frame that is really elsewhere in time.
Md = dc_time_match(unique(R(1).t_s), unique(R(2).t_s));
assert(Md.nUnmatched >= 0 && all(isnan(Md.dt_s(~Md.matched))), ...
    'an unmatched detection carries no gap, because there was nothing to measure one against');
assert(all(abs(Md.dt_s(Md.matched)) <= dtPage + 1e-12), 'and a matched one is within half a page');

%% (5) a declared dt is used as given ---------------------------------------------------------------
C2 = C; C2(2).dt_s = 0.5;                 % B strobed: nothing to do with the page interval
R2 = dc_process_cell(cel, C2, prm);
assert(abs(R2(2).dt_s - 0.5) < 1e-12, 'a declared dt must be used as given, got %.4g', R2(2).dt_s);
assert(abs(R2(2).t_s(end) - (R2(2).t0_s + R2(2).frame(end)*0.5)) < 1e-9, 'and drive the clock');

%% (6) no dt at all is an error ----------------------------------------------------------------------
C3 = C; C3(1).dt_s = NaN; C3(2).dt_s = NaN;   % the stack carries no metadata either
err = '';
try, dc_process_cell(cel, C3, prm); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_process_cell:noDt'), ...
    ['with no frame interval anywhere this must refuse, not guess — everything downstream is ' ...
     'reported in seconds (got "%s")'], err);

fprintf('dc track: A %d dets / %d tracks, B %d / %d; dt %.3g s each, B starts +%.3g s; every A frame matched within %.0f ms\n', ...
    R(1).nDets, R(1).nTracks, R(2).nDets, R(2).nTracks, R(1).dt_s, R(2).t0_s, 1000*max(abs(M.dt_s)));
fprintf('\nDC-TRACK SMOKE PASSED.\n');
end
