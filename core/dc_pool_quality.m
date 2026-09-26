function q = dc_pool_quality(sptPath, diamUm, pxUm, nSample, opts)
%DC_POOL_QUALITY  Pool DoG candidate qualities over evenly-spaced frames (drives the % tuner).
%
%   q = dc_pool_quality(sptPath, diamUm, pxUm, nSample)
%
% For each of nSample evenly-spaced frames, take the 5x5 local maxima above a LIGHT per-frame MAD
% floor (k=2, matching ERAware candidate_qualities) and collect their DoG values. The pooled
% vector `q` is the movie-level quality distribution: its histogram + a top-percentile cut give a
% single absolute threshold for the whole file (percentile detection is the proven, per-file way
% to set the threshold). Absolute threshold for keep_pct% = prctile(q, 100 - keep_pct).
%
%   sptPath : path to a multi-page single-particle TIFF
%   diamUm  : spot diameter in µm (default 0.5)
%   pxUm    : µm per pixel (default 0.10785)
%   nSample : number of frames to pool (default 40)
% OUTPUT q : pooled candidate DoG qualities (column vector)
if nargin<2 || isempty(diamUm),  diamUm  = 0.5; end
if nargin<3 || isempty(pxUm),    pxUm    = 0.10785; end
if nargin<4 || isempty(nSample), nSample = 40; end
% The tuner has to pool the SAME candidates detection will keep. With ridge rejection on, a
% filament contributes a long chain of maxima; leaving them in the pooled distribution would push
% the top-percentile threshold up and cost real spots — the threshold would be set by the very
% detections that are about to be thrown away.
if nargin<5 || ~isstruct(opts), opts = struct(); end
ridgeMax = []; if isfield(opts,'ridgeMax'), ridgeMax = opts.ridgeMax; end
sizeMax  = []; if isfield(opts,'sizeMax'),  sizeMax  = opts.sizeMax;  end
% Which frames the gate applies to, so the pooled distribution is the mixture detection will
% actually produce. Gating every sampled frame when the run gates only half of them sets the
% percentile off a distribution the run never sees.
bleedOn = 'all'; if isfield(opts,'bleedFrames') && ~isempty(opts.bleedFrames), bleedOn = lower(char(opts.bleedFrames)); end

% Sample only the pages belonging to THIS colour. A channel in dcSPT is a page list read from the
% acquisition's own slice labels, so it is passed in whole rather than reconstructed from a stride —
% a stride cannot express an acquisition that dropped a page, and pooling the other colour's frames
% into this one's distribution is exactly how a percentile threshold ends up set by the wrong signal.
info = imfinfo(sptPath);
if isfield(opts,'pages') && ~isempty(opts.pages)
    pages = round(opts.pages(:))';
    pages = pages(pages >= 1 & pages <= numel(info));
else
    stride = 1; if isfield(opts,'stride') && ~isempty(opts.stride), stride = max(1,round(opts.stride)); end
    offset = 0; if isfield(opts,'offset') && ~isempty(opts.offset), offset = max(0,round(opts.offset)); end
    pages = (offset+1) : stride : numel(info);
end
if isempty(pages), q = zeros(0,1); return; end
idx  = pages(unique(round(linspace(1, numel(pages), min(numel(pages), nSample)))));
q = zeros(0,1);
for k = idx
    f   = double(imread(sptPath, k));
    [dog, ~, dscale] = dc_dog(f, diamUm, pxUm);
    md  = median(dog(:)); madN = 1.4826*median(abs(dog(:)-md)) + 1e-6;
    thr = md + 2.0*madN;                          % light floor (k=2) = plausible candidates
    mx  = imdilate(dog, ones(5,5));
    cand = (dog >= mx) & (dog > thr);
    frameIdx = find(pages == k, 1);                  % page -> this colour's frame number
    gateThis = strcmp(bleedOn,'all') || ...
        (strcmp(bleedOn,'odd') && mod(round(frameIdx),2)==1) || ...
        (strcmp(bleedOn,'even') && mod(round(frameIdx),2)==0);
    if gateThis && ((~isempty(ridgeMax) && ridgeMax > 0) || (~isempty(sizeMax) && sizeMax > 0))
        [H_, W_] = size(dog);
        [cyq, cxq] = find(cand);
        ok = cxq > 2 & cxq <= W_-2 & cyq > 2 & cyq <= H_-2;   % spt_ridge needs a 3x3 neighbourhood
        cxq = cxq(ok); cyq = cyq(ok);
        [~, kr] = dc_ridge(dog, cxq, cyq, ridgeMax, dscale, sizeMax);
        q = [q; dog(sub2ind([H_ W_], cyq(kr), cxq(kr)))]; %#ok<AGROW>
    else
        q = [q; dog(cand)]; %#ok<AGROW>
    end
end
end
