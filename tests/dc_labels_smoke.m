function dc_labels_smoke()
%DC_LABELS_SMOKE  Reading the acquisition's own account of its pages, and building the colour maps
%from it — for two colours and for one.
%
% WHAT IS ASSERTED:
%   1. THE LABELS ARE READ FROM THE FILE, per page, as integers: which channel and which timepoint.
%      Not inferred from page parity — read.
%   2. A SINGLE-COLOUR STACK READS TOO. One channel is not a special case: c:1/1 with a timepoint per
%      page gives one colour whose frames are its pages, which is exactly what the single-colour
%      pipeline needs to stop calling frame numbers seconds.
%   3. A DROPPED PAGE IS SURVIVED, and this is the whole point. When an acquisition misses one page
%      the pages stop alternating, and a stride assumption silently hands half of one colour's frames
%      to the other. The labels still say what each page is, so fromStack stays right where
%      'interleaved' goes wrong — and the test shows 'interleaved' going wrong on the same file.
%   4. NO LABELS IS NOT AN ERROR: .ok is false, with a reason, and the caller falls back.
%   5. THE MAPS ARE INTEGERS THROUGHOUT: pages 1-based into the file, timepoints as the acquisition
%      numbered them, and dt_s derived from the page spacing only to be reported.
%
% Synthetic; reads no dataset.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'core'), fullfile(fileparts(here),'drivers'), here);
root = fullfile(tempdir, sprintf('dc_labels_%d', feature('getpid')));
if isfolder(root), rmdir(root,'s'); end
mkdir(root);
cleanup = onCleanup(@() rmdir(root,'s'));

dtPage = 0.02671184204518795;

%% (1) two colours, interleaved, as the microscope wrote them ---------------------------------------
nTp = 12; lab = {};
for t = 1:nTp
    lab{end+1} = sprintf('c:2/4 t:%d/%d - HVK-3C-Plate1-296-011 #1', t, nTp); %#ok<AGROW>
    lab{end+1} = sprintf('c:4/4 t:%d/%d - HVK-3C-Plate1-296-011 #1', t, nTp); %#ok<AGROW>
end
p2 = fullfile(root,'two.tif');
dc_ijtiff_write(p2, lab);

L = dc_tiff_labels(p2);
assert(L.ok && strcmp(L.source,'labels'), 'the labels should be read: %s', L.text);
assert(L.nPages == 2*nTp, '%d pages, got %d', 2*nTp, L.nPages);
assert(isequal(L.channels, [2 4]), 'the channels are the acquisition''s own numbers, got %s', mat2str(L.channels));
assert(isequal(L.ch(1:4)', [2 4 2 4]) && isequal(L.tp(1:4)', [1 1 2 2]), ...
    'channel and timepoint per page, got ch %s tp %s', mat2str(L.ch(1:4)'), mat2str(L.tp(1:4)'));
assert(L.nCh == 4 && L.nTp == nTp, 'and what the acquisition HAD: c:_/%g t:_/%g', L.nCh, L.nTp);

C = dc_channels('fromStack', p2, dtPage);
assert(numel(C) == 2 && strcmp(C(1).key,'c2') && strcmp(C(2).key,'c4'), ...
    'two colours keyed by the acquisition''s channel numbers');
assert(isequal(C(1).pages(1:3)', [1 3 5]) && isequal(C(2).pages(1:3)', [2 4 6]), ...
    'each colour''s frames are its own pages, got %s / %s', mat2str(C(1).pages(1:3)'), mat2str(C(2).pages(1:3)'));
assert(isequal(C(1).tp, (1:nTp)') && isequal(C(2).tp, (1:nTp)'), ...
    'and both are at every timepoint — one moment imaged twice');
assert(abs(C(1).dt_s - 2*dtPage) < 1e-12, ...
    'a frame of one colour is two pages on, so its dt is twice the page interval (%.6g)', C(1).dt_s);

A = dc_align(dc_dataset('new','cellA',C));
assert(strcmp(A.relation,'matched') && A.pageGap == 1, ...
    'the same timepoints, exposures one page apart (%s, gap %g)', A.relation, A.pageGap);

%% (2) one colour ------------------------------------------------------------------------------------
lab1 = arrayfun(@(t) sprintf('c:1/1 t:%d/8 - singlecolour #1', t), 1:8, 'uni', 0);
p1 = fullfile(root,'one.tif');
dc_ijtiff_write(p1, lab1);
L1 = dc_tiff_labels(p1);
assert(L1.ok && isequal(L1.channels, 1) && L1.nCh == 1, ...
    'one channel is not a special case: %s', L1.text);
assert(isequal(L1.tp', 1:8), 'its timepoints are its pages, got %s', mat2str(L1.tp'));
C1 = dc_channels('fromStack', p1, dtPage);
assert(numel(C1) == 1 && isequal(C1.pages', 1:8) && isequal(C1.tp', 1:8), ...
    'and the single colour maps frame k to page k, timepoint k');
assert(abs(C1.dt_s - dtPage) < 1e-12, 'with dt the page interval itself (%.6g)', C1.dt_s);

%% (3) a dropped page -------------------------------------------------------------------------------
% The acquisition misses the c:4 exposure of timepoint 3. Pages now go 2,4,2,4,2,2,4,... — still
% alternating plausibly for a while, which is exactly why a stride is dangerous.
labD = lab([1:4, 5, 7:end]);     % drop page 6, the c:4 of t:3
pD = fullfile(root,'dropped.tif');
dc_ijtiff_write(pD, labD);
LD = dc_tiff_labels(pD);
assert(LD.ok && LD.nPages == 2*nTp - 1, 'the short stack should still read (%d pages)', LD.nPages);
CD = dc_channels('fromStack', pD, dtPage);
c2 = CD(strcmp({CD.key},'c2')); c4 = CD(strcmp({CD.key},'c4'));
assert(c2.nFrames == nTp && c4.nFrames == nTp-1, ...
    'the colour that lost a page has one frame fewer (%d, %d)', c2.nFrames, c4.nFrames);
assert(~ismember(3, c4.tp), 'and no frame at the timepoint it missed');
assert(isequal(c4.tp', [1 2 4:nTp]), 'while keeping every other one, got %s', mat2str(c4.tp'));
assert(isequal(c2.pages', [1 3 5 6 8:2:(2*nTp-1)]), ...
    'the pages stop alternating after the gap, and the labels say so: %s', mat2str(c2.pages'));

AD = dc_align(dc_dataset('new','cellD',CD));
assert(strcmp(AD.relation,'partial') || strcmp(AD.relation,'subsampled'), ...
    'a colour missing a timepoint no longer matches the other exactly (got %s)', AD.relation);

% the same file, read by the assumption instead: it gets the pages wrong
CI = dc_channels('interleaved', LD.nPages, dtPage);
assert(~isequal(sort(CI(1).pages), sort(c2.pages)) || ~isequal(sort(CI(2).pages), sort(c4.pages)), ...
    ['a stride assumption should be shown to disagree with the file on this stack — if it agrees, ' ...
     'the fixture no longer exercises the dropped page']);

%% (4) no labels at all -----------------------------------------------------------------------------
pN = fullfile(root,'plain.tif');
imwrite(uint8(zeros(4)), pN); imwrite(uint8(zeros(4)), pN, 'WriteMode','append');
LN = dc_tiff_labels(pN);
assert(~LN.ok && strcmp(LN.source,'none'), 'a stack without labels is not an error');
assert(~isempty(LN.text), 'but says what it found: "%s"', LN.text);
err = ''; try, dc_channels('fromStack', pN); catch ME, err = ME.identifier; end
assert(strcmp(err,'dc_channels:noLabels'), ...
    'and fromStack refuses rather than guessing, pointing at the fallback (got "%s")', err);

fprintf('labels: %s\n', L.text);
fprintf('single colour: %s\n', L1.text);
fprintf('dropped page: %s -> c2 %d frames, c4 %d frames, %s\n', LD.text, c2.nFrames, c4.nFrames, AD.relation);
fprintf('\nDC-LABELS SMOKE PASSED.\n');
end
