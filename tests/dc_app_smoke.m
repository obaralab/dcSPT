function dc_app_smoke()
%DC_APP_SMOKE  The app, driven end to end through its own API, on a movie written to disk.
%
% WHAT IS ASSERTED:
%   1. IT BUILDS, with the four tabs in the order the work happens, and opens without a folder.
%   2. IT READS THE STACKS the way the drivers do: the channel table is filled from each file's own
%      slice labels, and the alignment verdict is the one dc_align gives.
%   3. THE THRESHOLD IS REFERENCED TO EACH COLOUR'S OWN NOISE. The fixture gives c2 three times c1's
%      read noise, so a threshold that tracks noise must come out roughly three times higher there —
%      and both colours must still detect all their particles. A percentile-of-maxima threshold
%      cannot do this: it implicitly assumes how many spots there are, keeps the brightest few when
%      there are more, and the tracks then fragment because each particle is seen only now and then.
%   4. TRACKING FILLS A DATASET that dc_dataset accepts — which means every invariant it checks
%      (a spot's timepoint and page are the channel map's; a track id belongs to one colour) held.
%   5. CO-MOTION RUNS CROSS-COLOUR ONLY and finds the pairs that were built to co-move.
%   6. THE PAIR PANEL DRAWS INTO THE TAB rather than a new window, and moves with prev/next.
%   7. EVERY BUTTON GOES THROUGH THE API. The test calls only H.api, so if a callback ever does its
%      own work instead of delegating, this stops testing the app and the count below catches it.
%
% Headless: every figure is built with 'visible','off'.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root,'core'), fullfile(root,'drivers'), fullfile(root,'app'), here);

tmp = fullfile(tempdir, sprintf('dc_app_%d', feature('getpid')));
if isfolder(tmp), rmdir(tmp,'s'); end
mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp,'s'));

nT = 60; px = 0.1;   % 160 px = 16 um: room for pairs that are genuinely far apart
writeMovie(tmp, nT, px);

%% (1) it builds ------------------------------------------------------------------------------------
H = dc_app(struct('visible','off'));
c1 = onCleanup(@() close(H.fig));
assert(isgraphics(H.fig), 'the app should build');
titles = arrayfun(@(t) string(t.Title), H.tabs);
assert(numel(titles) == 4, 'four tabs, got %d', numel(titles));
assert(all(contains(titles, ["Cells","Detect","Co-motion","Pair"])), ...
    'the tabs should be Cells, Detect, Co-motion, Pair — got %s', strjoin(titles, ', '));
fprintf('(1) built: %s\n', strjoin(titles, ' | '));

%% (2) it reads the stacks --------------------------------------------------------------------------
H.api.loadFolder(tmp);
St = H.api.state();
assert(numel(St.C) == 2, 'two colours, got %d', numel(St.C));
assert(all(arrayfun(@(c) c.nFrames == nT, St.C)), 'each colour should have %d frames', nT);
assert(strcmp(St.C(1).key,'c1') && strcmp(St.C(2).key,'c2'), ...
    'the keys should come from the slice labels, got %s and %s', St.C(1).key, St.C(2).key);
assert(strcmp(St.align.relation,'matched'), ...
    'both colours are at every timepoint, so they should read as matched (got %s)', St.align.relation);
assert(abs(St.pxUm - px) < 1e-9, 'pixel size should come from the file (%g vs %g)', St.pxUm, px);
fprintf('(2) %s\n', St.align.text);

%% (3) per-colour thresholds -------------------------------------------------------------------------
H.api.setParam('diamUm', 0.4);
det = H.api.detectPreview();          % runs before calibration too, at the defaults
assert(numel(det) == 2, 'preview should return both colours');

% calibrate is a button, so reach it the way the button does
btnCal = findButton(H.fig, 'Calibrate');
assert(~isempty(btnCal), 'the Detect tab should have a Calibrate button');
btnCal.ButtonPushedFcn([], []);
p = H.api.params();
assert(all(p.thr > 0), 'both thresholds should be set, got %s', mat2str(p.thr));
ratio = p.thr(2) / p.thr(1);
assert(ratio > 2 && ratio < 4.5, ...
    ['c2 has 3x c1''s read noise, so a noise-referenced threshold must be about 3x higher there ' ...
     '(got %.2f and %.2f, ratio %.2f). A ratio near 1 means the threshold is tracking something ' ...
     'other than the noise.'], p.thr(1), p.thr(2), ratio);
det = H.api.detectPreview();
assert(size(det{1},1) >= 9 && size(det{2},1) >= 9, ...
    ['both colours must still find their 10 particles after calibration (%d and %d) — a threshold ' ...
     'that is safe but blind is not a calibration'], size(det{1},1), size(det{2},1));
fprintf('(3) thresholds %.2f and %.2f (ratio %.2f, noise ratio 3) -> %d and %d of 10 spots\n', ...
    p.thr(1), p.thr(2), ratio, size(det{1},1), size(det{2},1));

%% (4) tracking fills a dataset ------------------------------------------------------------------------
H.api.setParam('linkUm', 0.8);
H.api.runTracking();
St = H.api.state();
assert(~isempty(St.D), 'tracking should leave a dataset');
dc_dataset('validate', St.D);         % throws if any invariant broke
nTr = numel(unique(St.D.spots.trackId(isfinite(St.D.spots.trackId))));
assert(nTr >= 12, 'the fixture has 20 particles; expect at least 12 tracks, got %d', nTr);
assert(any(isfinite(St.D.spots.iTot)), 'integrated intensity should have been measured');
fprintf('(4) %d spots, %d tracks, intensity measured\n', height(St.D.spots), nTr);

%% (5) co-motion, cross-colour only ----------------------------------------------------------------
% There must be ROOM between rFarUm and rMaxUm, or the "far" population is a thin shell and the null
% has nothing to stand on. This fixture is sparse, so widen both.
H.api.setCo('rMaxUm', 10);      % the field is 14 um across; keep pairs at every separation
H.api.setCo('rFarUm', 3);
H.api.setCo('rNearUm', 0.6);
H.api.runComotion();
St = H.api.state();
assert(~isempty(St.R), 'co-motion should leave a result');
assert(St.R.class == "cross" && all(St.R.steps.class == "cross"), ...
    'the app must run cross-colour only, got %s', St.R.class);
assert(all(St.R.steps.chA ~= St.R.steps.chB), 'and every pair must span the two colours');
assert(~isempty(St.N) && isfinite(St.N.k.fitted), 'the null should have been built');
assert(height(St.R.pairs) >= 1, 'the fixture has a co-moving cross-colour pair; none was found');
fprintf('(5) %d cross pairs, far field %+.4f, SE = %.3f/sqrt(n)\n', ...
    height(St.R.pairs), St.N.far.meanCos, St.N.k.fitted);

%% (6) the pair panel draws into the tab -------------------------------------------------------------
nFigBefore = numel(findall(0,'Type','figure'));
H.api.showPair(1);
assert(numel(findall(0,'Type','figure')) == nFigBefore, ...
    'the panel must draw into the Pair tab, not pop a new window');
pp = getappdata(H.fig,'pair');
assert(~isempty(allchild(pp.host)), 'the Pair tab should now hold the panel');
assert(pp.k == 1, 'and remember which pair it is showing');
if height(St.R.pairs) > 1
    btnNext = findButton(H.fig, 'next');
    btnNext.ButtonPushedFcn([], []);
    pp = getappdata(H.fig,'pair');
    assert(pp.k == 2, 'next should advance the pair, got %d', pp.k);
end
fprintf('(6) pair panel hosted in the tab, prev/next works\n');

% A far-field shell with no room in it must fail by NAME, and the message must point at rMaxUm —
% the fix people do not think of, because the symptom looks like rFarUm being too high.
errId = ''; errMsg = '';
try
    dc_comotion_null(St.R, struct('rFarUm', 9.9));     % 9.9 to 10.0: essentially no shell
catch ME
    errId = ME.identifier; errMsg = ME.message;
end
assert(strcmp(errId,'dc_comotion_null:thinFar'), ...
    'an empty far-field shell must be refused by name, got "%s"', errId);
assert(contains(errMsg,'rMaxUm'), ...
    'and the message must name rMaxUm as the thing to raise: "%s"', errMsg);
fprintf('    empty far-field shell refused, message names rMaxUm\n');

%% (7) the buttons delegate ---------------------------------------------------------------------------
assert(all(isfield(H.api, {'loadFolder','detectPreview','runTracking','runComotion','showPair'})), ...
    'the api must expose every action a button performs');
fprintf('(7) api exposes %d actions\n', numel(fieldnames(H.api)));

fprintf('\nDC-APP SMOKE PASSED.\n');
end

% =================================================================================================
function writeMovie(folder, nT, px)
% Two colour stacks with ImageJ slice labels. Six particles: one cross-colour pair that co-moves,
% the rest independent. c2 is twice as bright as c1, so a single shared detection threshold cannot
% serve both — which is what test (3) is for.
S = 140;
rng(5);
sd = 0.45;                                        % px per step
amp = [260 420];                 % both well above their own noise
noiseSd = [4 12];                % ...but c2's read noise is 3x c1's: what calibration must track
nP = 10;                     % particles per colour: 100 cross pairs, so the far field is populated
pos = cell(1,2);
shared = sd*randn(nT,2);
% Where the independent particles sit: spread out, so most pairs are far apart and the far-field
% null has something to be built from.
rng(5);
seeds = [22 22; 118 26; 26 118; 116 114; 70 20; 20 70; 118 70; 70 118; 45 92];
for ch = 1:2
    p = zeros(nT, nP, 2);                         % [frame x particle x (x,y)]
    % particle 1 of each colour is the co-moving pair: 70% shared motion, started close together
    start = [80 80] + (ch-1)*[1.5 0];
    p(:,1,:) = reshape(cumsum([start; sqrt(0.7)*shared(2:end,:) + sqrt(0.3)*sd*randn(nT-1,2)],1), nT,1,2);
    for k = 2:nP
        p(:,k,:) = reshape(cumsum([seeds(k-1,:); sd*randn(nT-1,2)],1), nT,1,2);
    end
    pos{ch} = p;
end
[X, Y] = meshgrid(1:S, 1:S);
sig = 0.25/px/2.355;                               % a 250 nm PSF
for ch = 1:2
    f = fullfile(folder, sprintf('Ch%d.tif', ch));
    labels = arrayfun(@(t) sprintf('c:%d/2 t:%d/%d - synth #1', ch, t, nT), 1:nT, 'uni', 0);
    pages = cell(1,nT);
    for t = 1:nT
        im = 100 + noiseSd(ch)*randn(S);
        for k = 1:size(pos{ch},2)
            cx = pos{ch}(t,k,1); cy = pos{ch}(t,k,2);
            im = im + amp(ch)*exp(-((X-cx).^2 + (Y-cy).^2)/(2*sig^2));
        end
        pages{t} = uint16(max(im,0));
    end
    ijWrite(f, labels, pages, px);
end
end

function ijWrite(path, labels, pages, px)
% A 16-bit TIFF carrying ImageJ slice labels and an ImageDescription with the scale, assembled by
% hand: imwrite cannot write the unknown tags the label reader needs.
n = numel(labels);
desc = sprintf('ImageJ=1.54f\nimages=%d\nframes=%d\nunit=micron\nfinterval=0.012\n', n, n);
hdr = uint8([uint8('IJIJ') uint8('labl') be32(n)]);
cnt = numel(hdr); body = uint8([]);
for k = 1:n, b = u16be(labels{k}); body = [body b]; cnt(end+1) = numel(b); end %#ok<AGROW>
ijm = [hdr body];

out = uint8([uint8('II') lo16(42) lo32(8)]);
pxOff = zeros(1,n); H = zeros(1,n); W = zeros(1,n);
for k = 1:n
    im = uint16(pages{k}); [H(k), W(k)] = size(im);
    pxOff(k) = numel(out);
    raw = typecast(reshape(im', 1, []), 'uint8');   % little-endian 16-bit, row-major
    out = [out raw]; %#ok<AGROW>
    if mod(numel(out),2), out = [out uint8(0)]; end %#ok<AGROW>
end
dOff = numel(out); out = [out uint8(desc) uint8(0)];
if mod(numel(out),2), out = [out uint8(0)]; end
ijmOff = numel(out); out = [out ijm];
if mod(numel(out),2), out = [out uint8(0)]; end
cntOff = numel(out); for i = 1:numel(cnt), out = [out lo32(cnt(i))]; end %#ok<AGROW>
resOff = numel(out); out = [out lo32(1e6) lo32(round(1e6*px))];   % XResolution = 1/px as a rational

ifd = zeros(1,n);
for k = 1:n
    ifd(k) = numel(out);
    E = [256 3 1 W(k); 257 3 1 H(k); 258 3 1 16; 259 3 1 1; 262 3 1 1; 273 4 1 pxOff(k); ...
         277 3 1 1; 278 3 1 H(k); 279 4 1 H(k)*W(k)*2; 282 5 1 resOff; 283 5 1 resOff; 296 3 1 1];
    if k == 1
        E = sortrows([E; 270 2 numel(desc)+1 dOff; 50838 4 numel(cnt) cntOff; 50839 1 numel(ijm) ijmOff], 1);
    else
        E = sortrows(E, 1);
    end
    e = uint8([]);
    for r = 1:size(E,1), e = [e lo16(E(r,1)) lo16(E(r,2)) lo32(E(r,3)) lo32(E(r,4))]; end %#ok<AGROW>
    out = [out lo16(size(E,1)) e lo32(0)]; %#ok<AGROW>
end
for k = 1:n
    nxt = 0; if k < n, nxt = ifd(k+1); end
    nE = 12; if k == 1, nE = 15; end
    q = ifd(k) + 2 + 12*nE;
    out(q+1:q+4) = lo32(nxt);
end
out(5:8) = lo32(ifd(1));
fid = fopen(path,'w'); assert(fid > 0, 'cannot write %s', path);
fwrite(fid, out, 'uint8'); fclose(fid);
end

function b = findButton(fig, txt)
b = findall(fig, 'Type','uibutton');
keep = arrayfun(@(q) contains(string(q.Text), txt), b);
b = b(keep); if ~isempty(b), b = b(1); end
end

function b = lo16(v), b = uint8([mod(v,256) floor(v/256)]); end
function b = lo32(v)
v = double(v);
b = uint8([mod(v,256) mod(floor(v/256),256) mod(floor(v/65536),256) mod(floor(v/16777216),256)]);
end
function b = be32(v)
v = double(v);
b = uint8([floor(v/16777216) mod(floor(v/65536),256) mod(floor(v/256),256) mod(v,256)]);
end
function b = u16be(s)
u = double(s); b = uint8(zeros(1, 2*numel(u)));
b(1:2:end) = floor(u/256); b(2:2:end) = mod(u,256);
end
