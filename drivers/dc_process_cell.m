function R = dc_process_cell(cel, C, prm)
%DC_PROCESS_CELL  Detect and track BOTH colours of one cell, each on its own clock.
%
%   R = dc_process_cell(cel, C, prm)
%
% cel : struct with .stack (the movie; interleaved) or .stacks (one path per colour), and .base
% C   : the colours, from dc_channels
% prm : .diamUm .thrAbs .pxUm .linkUm .gapUm .maxGap, optional .maxFrames .detOpts
%
% ONE CALL, BOTH COLOURS. The single-colour pipeline this toolkit grew out of tracks one thing per
% run, so two colours meant running it twice and hoping the two runs agreed about pixel size, frame
% interval and which pages were whose — they were written to the same file names and the second
% quietly replaced the first. Here the pair is the unit: the colours are read from one description,
% tracked independently, and returned together with the clock each one is on.
%
% WHAT "INDEPENDENTLY" MEANS. Each colour gets its own detection, its own linking and its own frame
% numbering; neither is used to help find the other. That matters: a tracker that let one colour
% inform the other would manufacture exactly the correlation this experiment is trying to measure.
% The colours meet later, in seconds, and only to be compared.
%
% OUTPUT R : 1xN struct, one per colour —
%   .key .label            which colour
%   .frame .x .y .q        every detection: frame 0-based WITHIN this colour, x/y in MICRONS
%                          (detector pixels times pxUm — the unit every downstream ...Um option means)
%   .iMean .iMax .iTot     intensity on the RAW frame over a disk of half the detection diameter.
%                          iTot is what a bleaching step count reads: a step there is one
%                          fluorophore's photons leaving, while a mean also tracks the background.
%   .tp .page              the acquisition's own indices for each detection: which timepoint of the
%                          experiment it belongs to, and which page of the stack it came from
%   .trackId .spotId       linkage (NaN trackId = detected, not tracked)
%   .tracks                the raw chains, as dc_track returned them
%   .dt_s                  seconds per frame, carried for REPORTING a result — never an index
%   .pages .tp_all         this colour's full frame -> page and frame -> timepoint maps
%   .nFrames .nTracks .nDets

assert(~isempty(C), 'dc_process_cell:noChannels', 'no colours to process');
px = prm.pxUm;
R = repmat(emptyR(), 1, numel(C));

for c = 1:numel(C)
    ch = C(c);
    stack = stackFor(cel, ch);
    assert(isfile(stack), 'dc_process_cell:noStack', 'colour %s: no stack at %s', ch.key, stack);
    % The pages are the channel map's, read from the acquisition rather than computed from a stride.
    nfr = ch.nFrames;
    if isfield(prm,'maxFrames') && ~isempty(prm.maxFrames)
        nfr = min(nfr, prm.maxFrames);     % both maps truncate together, or they would disagree
    end
    pages = ch.pages(1:nfr)';
    tps   = ch.tp(1:nfr);
    assert(nfr > 0, 'dc_process_cell:noFrames', 'colour %s: no pages', ch.key);

    % dt is carried for REPORTING only — nothing here indexes by it. A colour with none is fine
    % until something asks for a result in seconds.
    dt = ch.dt_s;
    if ~isfinite(dt) || dt <= 0
        dtPage = NaN;
        try, tc = dc_tiff_calib(stack); if isstruct(tc) && isfield(tc,'dt_s'), dtPage = tc.dt_s; end, catch, end
        if isfinite(dtPage) && dtPage > 0 && numel(pages) > 1
            dt = median(diff(pages)) * dtPage;
        end
    end

    [readPage, closeStack] = dc_tiff_pages(stack);
    cleanup = onCleanup(closeStack);
    dets = cell(1, nfr);
    fr = {}; xs = {}; ys = {}; qs = {}; im = {};
    % The disk the intensity is summed over. Tied to the detection diameter so it follows the PSF
    % rather than being a second, independently wrong number.
    rPx = max(1.5, 0.5 * prm.diamUm / px);
    dopts = struct(); if isfield(prm,'detOpts') && isstruct(prm.detOpts), dopts = prm.detOpts; end
    for t = 1:nfr
        raw = double(readPage(pages(t)));
        xy = dc_detect(raw, prm.diamUm, px, prm.thrAbs, dopts);
        dets{t} = xy;
        n = size(xy,1);
        if n > 0
            fr{end+1} = repmat(t-1, n, 1); xs{end+1} = xy(:,1); ys{end+1} = xy(:,2); qs{end+1} = xy(:,3); %#ok<AGROW>
            % MEAN/MAX/TOTAL on the RAW frame. Bleaching is counted on TOTAL: a step there is one
            % fluorophore's worth of photons leaving, whereas a mean moves with the disk's
            % background as well.
            % dc_measure indexes the raw frame, so it takes PIXEL centres — the un-scaled xy.
            im{end+1} = dc_measure(raw, xy(:,1:2), rPx); %#ok<AGROW>
        end
    end
    clear cleanup

    [tracks, tinfo] = dc_track(dets, prm.linkUm, prm.gapUm, prm.maxGap, px);

    frame = cat(1, fr{:}); x = cat(1, xs{:}); y = cat(1, ys{:}); q = cat(1, qs{:});
    I = cat(1, im{:});
    if isempty(frame), frame = zeros(0,1); x = frame; y = frame; q = frame; I = zeros(0,3); end
    spotId = (0:numel(frame)-1)';
    trackId = nan(numel(frame),1);
    % x/y are still DETECTOR PIXELS at this point, and must stay that way until the tracks have been
    % matched below: dc_track works in pixels and returns pixel positions, and the match is made by
    % formatting the coordinates into a key. Scaling one side and not the other silently matches
    % nothing, and every track comes back with no id.
    % Label each detection with the track it ended up in, by (frame, x, y) — the tracker returns
    % positions, not indices back into the detection list.
    keyOf = containers.Map('KeyType','char','ValueType','double');
    for i = 1:numel(frame)
        keyOf(sprintf('%d|%.6f|%.6f', frame(i), x(i), y(i))) = i;
    end
    for k = 1:numel(tracks)
        tr = tracks{k};
        for r = 1:size(tr,1)
            kk = sprintf('%d|%.6f|%.6f', tr(r,1)-1, tr(r,2), tr(r,3));
            if isKey(keyOf, kk), trackId(keyOf(kk)) = k - 1; end
        end
    end

    % NOW to microns, once the pixel-keyed matching is done. Everything downstream is named ...Um
    % and documented in um — dc_steps, dc_comotion's radii, dc_masks — so leaving detector pixels
    % here would make all of those silently mean pixels: on a 0.1 um/px camera a radius of 5 would
    % be half a micron rather than five.
    x = x * px;  y = y * px;

    R(c) = struct('key',ch.key, 'label',ch.label, ...
        'frame',frame, 'x',x, 'y',y, 'q',q, ...
        'iMean',I(:,1), 'iMax',I(:,2), 'iTot',I(:,3), ...
        'tp', tps(frame+1), 'page', pages(frame+1)', ...
        'trackId',trackId, 'spotId',spotId, 'tracks',{tracks}, ...
        'dt_s',dt, 'pages',pages(:), 'tp_all',tps(:), ...
        'stack',stack, 'base',baseOf(cel), 'pxUm',px, ...
        'nFrames',nfr, 'nTracks',numel(tracks), 'nDets',tinfo.nDets);
end
end

% =================================================================================================
function s = stackFor(cel, ch)
% Interleaved colours share the cell's stack; paired colours name their own by suffix.
if isfield(cel,'stacks') && isstruct(cel.stacks) && isfield(cel.stacks, ch.key)
    s = cel.stacks.(ch.key); return
end
s = cel.stack;
end

function b = baseOf(cel)
if isfield(cel,'base') && ~isempty(cel.base), b = char(cel.base); return, end
[~, b] = fileparts(cel.stack);
end

function R = emptyR()
R = struct('key','', 'label','', 'frame',[], 'x',[], 'y',[], 'q',[], ...
    'iMean',[], 'iMax',[], 'iTot',[], 'tp',[], 'page',[], ...
    'trackId',[], 'spotId',[], 'tracks',{{}}, 'dt_s',NaN, 'pages',[], 'tp_all',[], ...
    'stack','', 'base','', 'pxUm',NaN, ...
    'nFrames',0, 'nTracks',0, 'nDets',0);
end
