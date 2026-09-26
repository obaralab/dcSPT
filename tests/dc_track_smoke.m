function dc_track_smoke()
%DC_TRACK_SMOKE  Two colours of one interleaved stack, detected and tracked independently, each with
%its own frame numbering and its own page and timepoint maps.
%
% WHAT IS ASSERTED:
%   1. BOTH COLOURS COME BACK from one call, each with its own frames, tracks and detections.
%   2. EACH IS TRACKED ON ITS OWN PAGES: colour 1's particle is found where colour 1's particle is,
%      and colour 2's where colour 2's is. Neither picks up the other's.
%   3. EVERY DETECTION CARRIES THE ACQUISITION'S OWN INDICES — which page it came from and which
%      timepoint it belongs to — and they are the channel map's, not recomputed. A run in which both
%      colours claimed the same pages would look identical to this one in every other respect.
%   4. THE TWO CAN BE PUT ON ONE CLOCK, and it is the ACQUISITION's: their timepoints, compared as
%      sets. Interleaved colours share every timepoint while sitting on consecutive PAGES, and the
%      page gap is what says the two exposures were not simultaneous.
%   5. A COLOUR MAY DECLARE ITS OWN dt and it is carried as given — for reporting a result in
%      seconds, never for indexing.
%   6. NO dt ANYWHERE IS NOT AN ERROR. Nothing in detection, linking or the merged view divides by
%      seconds, so a colour whose interval is unknown still tracks; the seconds are needed only when
%      a diffusion coefficient is finally asked for.
%
% Synthetic; reads no dataset.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'core'), fullfile(fileparts(here),'drivers'));
root = fullfile(tempdir, sprintf('dc_track_%d', feature('getpid')));
if isfolder(root), rmdir(root,'s'); end
mkdir(root);
cleanup = onCleanup(@() rmdir(root,'s'));

%% an interleaved stack: colour 1 drifts left to right on the odd pages, colour 2 sits still on the even
S = 64; nPg = 60; dtPage = 0.01; px = 0.1;
mv = fullfile(root,'cellA.tif');
rng(4);
for pg = 1:nPg
    im = uint16(200 + 6*randn(S,S));
    if mod(pg,2) == 1
        f = (pg+1)/2;  x = 10 + round(f*0.7); y = 20;    % a clean drift
    else
        f = pg/2;      x = 48; y = 44 - round(f*0.2);    % barely moves, elsewhere in the frame
    end
    im(y-1:y+1, x-1:x+1) = im(y-1:y+1, x-1:x+1) + uint16(5000);
    if pg == 1, imwrite(im, mv); else, imwrite(im, mv, 'WriteMode','append'); end
end
cel = struct('stack', mv, 'base', 'cellA');
prm = struct('diamUm',0.4,'thrAbs',60,'pxUm',px,'linkUm',1.5,'gapUm',1.5,'maxGap',1);

%% (1)(2) both colours, each on its own pages -------------------------------------------------------
C = dc_channels('interleaved', nPg, dtPage);
R = dc_process_cell(cel, C, prm);
assert(numel(R) == 2, 'both colours should come back, got %d', numel(R));
assert(strcmp(R(1).key,'c1') && strcmp(R(2).key,'c2'), 'and keep their keys');
assert(R(1).nFrames == nPg/2 && R(2).nFrames == nPg/2, ...
    'each colour has half the pages (%d, %d of %d)', R(1).nFrames, R(2).nFrames, nPg);
assert(R(1).nDets > 20 && R(2).nDets > 20, 'both colours should be detected (%d, %d)', R(1).nDets, R(2).nDets);
% Positions come back in MICRONS. The fixture puts colour 1's particle between x = 10 and 31 px at
% 0.1 um/px, so anything near 10-31 would mean the detector's pixels leaked through unscaled.
assert(max(R(1).x) < 0.9*size(imread(mv,1),2)*px, ...
    'x must be in um, not pixels (max %.2f against a %.1f um field)', max(R(1).x), size(imread(mv,1),2)*px);
assert(range(R(1).x) > 1.0, 'colour 1 should travel (range %.2f um)', range(R(1).x));
assert(range(R(2).x) < 0.5, 'colour 2 should sit still (range %.2f um)', range(R(2).x));
assert(min(R(2).x) > max(R(1).x), ...
    'the two colours occupy different parts of the frame; neither should pick up the other');
assert(R(1).nTracks >= 1 && R(2).nTracks >= 1, 'each colour should link into at least one track');

%% (3) the indices ----------------------------------------------------------------------------------
assert(isequal(R(1).pages', 1:2:nPg) && isequal(R(2).pages', 2:2:nPg), ...
    'colour 1 is the odd pages and colour 2 the even ones, got %s / %s', ...
    mat2str(R(1).pages(1:3)'), mat2str(R(2).pages(1:3)'));
assert(isequal(R(1).tp', 1:nPg/2) && isequal(R(2).tp', 1:nPg/2), ...
    'and both are at every timepoint of the acquisition');
for c = 1:2
    assert(isequal(R(c).page, C(c).pages(R(c).frame+1)) && isequal(R(c).tp, C(c).tp(R(c).frame+1)), ...
        ['colour %d: a detection''s page and timepoint must be the channel map''s for its frame, ' ...
         'not a second copy computed here'], c);
    assert(all(R(c).frame >= 0 & R(c).frame < R(c).nFrames), 'colour %d: frames are 0-based within the colour', c);
end
% the pixel a detection came from is the page it is labelled with — read the raw page back and check
[readPage, closeStack] = dc_tiff_pages(mv);
cl2 = onCleanup(closeStack);
i1 = find(R(1).frame == 5, 1);
pg = R(1).page(i1);
raw = double(readPage(pg));
% R.x/R.y are MICRONS, so they divide by the pixel size to index the raw frame. If this ever reads
% the array directly again, positions have silently gone back to detector pixels.
assert(raw(round(R(1).y(i1)/px), round(R(1).x(i1)/px)) > 1000, ...
    'the detection at frame 5 should be bright on page %d, the page it says it came from', pg);
clear cl2

%% (4) one clock for the pair -----------------------------------------------------------------------
D = dc_dataset('new','cellA',C);
D = dc_dataset('addChannel', D, R(1));
D = dc_dataset('addChannel', D, R(2));
A = dc_align(D);
assert(strcmp(A.relation,'matched'), ...
    'interleaved colours share every timepoint (got %s)', A.relation);
assert(A.nShared == nPg/2, 'all %d of them, got %d', nPg/2, A.nShared);
assert(A.pageGap == 1, ...
    ['and the page gap is what says they were not simultaneous — one page apart, got %g'], A.pageGap);
% the merged view at one timepoint holds both colours' own frames there, which are the same NUMBER
% here but different PAGES — and it is the pages that were actually exposed.
[S4, i4] = dc_merge(D, struct('tp', 6));
assert(all(S4.tp == 6) && isempty(i4.missing), 'both colours are present at timepoint 6: %s', i4.text);
pgs = unique(S4.page);
assert(numel(pgs) == 2 && abs(diff(pgs)) == 1, ...
    'from two consecutive pages, got %s', mat2str(pgs'));

%% (5) a declared dt is carried as given ------------------------------------------------------------
C2 = C; C2(2).dt_s = 0.5;                 % colour 2 strobed: nothing to do with the page interval
R2 = dc_process_cell(cel, C2, prm);
assert(abs(R2(2).dt_s - 0.5) < 1e-12, 'a declared dt must be carried as given, got %.4g', R2(2).dt_s);
assert(isequal(R2(2).page, R(2).page), ...
    'and must change NOTHING about the indices — the seconds are for reporting, not for finding pages');

%% (6) no dt at all still tracks --------------------------------------------------------------------
C3 = C; C3(1).dt_s = NaN; C3(2).dt_s = NaN;
R3 = dc_process_cell(cel, C3, prm);
assert(R3(1).nDets == R(1).nDets && isequal(R3(1).page, R(1).page), ...
    ['a colour with no frame interval must still detect and track identically: nothing in this ' ...
     'pipeline indexes by seconds']);

fprintf('dc track: c1 %d dets / %d tracks on pages %d..%d, c2 %d / %d on pages %d..%d\n', ...
    R(1).nDets, R(1).nTracks, R(1).pages(1), R(1).pages(end), ...
    R(2).nDets, R(2).nTracks, R(2).pages(1), R(2).pages(end));
fprintf('align: %s\n', A.text);
fprintf('\nDC-TRACK SMOKE PASSED.\n');
end
