function dc_masks_smoke()
%DC_MASKS_SMOKE  One organelle mask, shared by both colours, indexed by TIMEPOINT — and the
%assumption that makes sharing safe, checked rather than assumed.
%
% WHAT IS ASSERTED:
%   1. THE PAGE IS A FUNCTION OF THE TIMEPOINT, not of frame number. Both colours at one timepoint
%      get the SAME mask page even though that timepoint is a different frame number in each — which
%      is the whole reason the lookup takes a timepoint.
%   2. THE PAGE IS THE ONE CONTAINING THE TIMEPOINT. A window is a half-open range, so the last
%      timepoint of a page and the first of the next fall either side of the boundary. Nearest-centre
%      would put a boundary timepoint on the wrong page by half a window.
%   3. THE MASK OUTLASTING THE MOVIE IS NORMAL: a timepoint past the last page clamps rather than
%      failing.
%   4. THE SHARING ASSUMPTION IS CHECKED. Where the organelle is averaged over many timepoints,
%      sharing is declared safe and says over how many. Where it is not, that is said plainly.
%   5. A MASK WITH NO DECLARED WINDOW IS REFUSED, because deriving one from page counts needs an
%      exact integer ratio real acquisitions do not give, and the fallback is silently 1.
%
% Synthetic; reads no dataset.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'core'), fullfile(fileparts(here),'drivers'), here);
root = fullfile(tempdir, sprintf('dc_masks_%d', feature('getpid')));
if isfolder(root), rmdir(root,'s'); end
mkdir(root);
cleanup = onCleanup(@() rmdir(root,'s'));

dt = 0.0267118;
% 'a' on the odd pages at every timepoint; 'b' on every fourth page, at the odd timepoints only.
a = mkCh('a','fast', (1:2:1000)', (1:500)',   2*dt);
b = mkCh('b','slow', (2:4:1000)', (1:2:499)', 4*dt);
D = dc_dataset('new','cellA',[a b]);

% a mito stack of 6 pages, each averaging 100 timepoints
S = 32; nPg = 6; tpPerPage = 100;
p = fullfile(root,'mito.tif');
for k = 1:nPg
    m = false(S); m(6+k:20+k, 8:24) = true;       % it moves a little, page to page
    if k == 1, imwrite(uint8(m)*255, p); else, imwrite(uint8(m)*255, p, 'WriteMode','append'); end
end
D = dc_masks('set', D, struct('mito', p, 'tpPerPage', tpPerPage, 'tp0', 1));

%% (1) both colours at one timepoint get one page ---------------------------------------------------
tp = 251;
fa = find(D.channels(1).tp == tp) - 1;             % 0-based frame of each colour at that timepoint
fb = find(D.channels(2).tp == tp) - 1;
assert(~isempty(fa) && ~isempty(fb) && fa ~= fb, ...
    'the fixture should put this timepoint at different FRAME numbers in the two colours (%d, %d)', fa, fb);
[pgA, span] = dc_masks('page', D, D.channels(1).tp(fa+1));
pgB         = dc_masks('page', D, D.channels(2).tp(fb+1));
assert(pgA == pgB, ...
    ['one timepoint must give one mask page in both colours (%d vs %d) — indexing by frame number ' ...
     'would have given two'], pgA, pgB);
assert(pgA == 3 && isequal(span,[201 300]), ...
    'timepoint 251 is on page 3, which averages 201..300 (got page %d, %s)', pgA, mat2str(span));
m1 = dc_masks('at', D, 'mito', tp);
assert(islogical(m1) && any(m1(:)), 'the page should come back as a mask');

%% (2) the window contains, it does not round -------------------------------------------------------
assert(dc_masks('page', D, 200) == 2 && dc_masks('page', D, 201) == 3, ...
    'the boundary between windows is a step, not a rounding: 200 -> page 2, 201 -> page 3');
assert(dc_masks('page', D, 1) == 1, 'the first timepoint is on the first page');
assert(dc_masks('page', D, 100) == 1, 'and so is the last one that page averages');

%% (3) past the end of the mask stack ---------------------------------------------------------------
mLate = dc_masks('at', D, 'mito', 20000);
assert(islogical(mLate) && any(mLate(:)), ...
    'a mask stack shorter than the movie is normal and must clamp, not fail');

%% (4) the assumption, checked ----------------------------------------------------------------------
C = dc_masks('check', D);
assert(C.ok, 'a mask averaging 100 timepoints should be safe to share: "%s"', C.text);
assert(C.tpPerPage == tpPerPage, 'and report the window (%g)', C.tpPerPage);
assert(contains(C.text,'has not moved between them'), 'and say why: "%s"', C.text);
% averaged over a single timepoint: no slower than the particles, so the assumption no longer holds
D2 = dc_masks('set', D, struct('mito', p, 'tpPerPage', 1, 'tp0', 1));
C2 = dc_masks('check', D2);
assert(~C2.ok, 'a mask averaging one timepoint must NOT be declared safe to share');
assert(contains(C2.text,'no longer safe'), 'and must say so plainly: "%s"', C2.text);

%% (5) no declared window --------------------------------------------------------------------------
err = ''; try, dc_masks('set', D, struct('mito', p)); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_masks:noWindow'), 'a mask with no averaging window must be refused, got "%s"', err);

fprintf('masks: one page (%d, averaging %d..%d) for both colours at timepoint %d, though their frame numbers there are %d and %d\n', ...
    pgA, span(1), span(2), tp, fa, fb);
fprintf('check: %s\n', C.text);
fprintf('\nDC-MASKS SMOKE PASSED.\n');
end

% =================================================================================================
function c = mkCh(key, label, pages, tp, dt)
c = struct('key',key, 'label',label, 'pages',pages(:), 'tp',tp(:), 'dt_s',dt, 'nFrames',numel(pages));
end
