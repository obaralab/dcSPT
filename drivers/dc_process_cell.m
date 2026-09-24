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
%   .frame .x .y .q        every detection (frame is 0-based WITHIN this colour)
%   .t_s                   the experiment clock: t0_s + frame * dt_s
%   .trackId .spotId       linkage (NaN trackId = detected, not tracked)
%   .tracks                the raw chains, as dc_track returned them
%   .dt_s .t0_s            this colour's clock
%   .stride .offset .pages how its frames sit in the stack it came from
%   .nFrames .nTracks .nDets

assert(~isempty(C), 'dc_process_cell:noChannels', 'no colours to process');
px = prm.pxUm;
R = repmat(emptyR(), 1, numel(C));

for c = 1:numel(C)
    ch = C(c);
    stack = stackFor(cel, ch);
    assert(isfile(stack), 'dc_process_cell:noStack', 'colour %s: no stack at %s', ch.key, stack);
    info = imfinfo(stack); nPages = numel(info);
    pages = (ch.offset+1) : ch.stride : nPages;
    if isfield(prm,'maxFrames') && ~isempty(prm.maxFrames)
        pages = pages(1:min(numel(pages), prm.maxFrames));
    end
    nfr = numel(pages);
    assert(nfr > 0, 'dc_process_cell:noFrames', ...
        'colour %s: offset %d and stride %d select no pages from %d', ch.key, ch.offset, ch.stride, nPages);

    dt = ch.dt_s;
    if ~isfinite(dt) || dt <= 0
        % Fall back to the stack's own page interval x stride, and say that is what happened — a
        % derived dt is a guess about the acquisition, not a reading of it.
        dtPage = NaN;
        try, tc = dc_tiff_calib(stack); if isstruct(tc) && isfield(tc,'dt_s'), dtPage = tc.dt_s; end, catch, end
        assert(isfinite(dtPage) && dtPage > 0, 'dc_process_cell:noDt', ...
            ['colour %s declares no frame interval and %s carries none in its metadata. Every ' ...
             'diffusion coefficient and dwell time is reported in seconds, so this cannot be guessed.'], ...
            ch.key, stack);
        dt = dtPage * ch.stride;
    end

    [readPage, closeStack] = dc_tiff_pages(stack);
    cleanup = onCleanup(closeStack);
    dets = cell(1, nfr);
    fr = {}; xs = {}; ys = {}; qs = {};
    dopts = struct(); if isfield(prm,'detOpts') && isstruct(prm.detOpts), dopts = prm.detOpts; end
    for t = 1:nfr
        raw = double(readPage(pages(t)));
        xy = dc_detect(raw, prm.diamUm, px, prm.thrAbs, dopts);
        dets{t} = xy;
        n = size(xy,1);
        if n > 0
            fr{end+1} = repmat(t-1, n, 1); xs{end+1} = xy(:,1); ys{end+1} = xy(:,2); qs{end+1} = xy(:,3); %#ok<AGROW>
        end
    end
    clear cleanup

    [tracks, tinfo] = dc_track(dets, prm.linkUm, prm.gapUm, prm.maxGap, px);

    frame = cat(1, fr{:}); x = cat(1, xs{:}); y = cat(1, ys{:}); q = cat(1, qs{:});
    if isempty(frame), frame = zeros(0,1); x = frame; y = frame; q = frame; end
    spotId = (0:numel(frame)-1)';
    trackId = nan(numel(frame),1);
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

    R(c) = struct('key',ch.key, 'label',ch.label, ...
        'frame',frame, 'x',x, 'y',y, 'q',q, 't_s', ch.t0_s + frame*dt, ...
        'trackId',trackId, 'spotId',spotId, 'tracks',{tracks}, ...
        'dt_s',dt, 't0_s',ch.t0_s, 'stride',ch.stride, 'offset',ch.offset, 'pages',pages, ...
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
if ~isempty(ch.file)
    [d, b, e] = fileparts(s);
    cand = fullfile(d, [b ch.file e]);
    if isfile(cand), s = cand; end
end
end

function b = baseOf(cel)
if isfield(cel,'base') && ~isempty(cel.base), b = char(cel.base); return, end
[~, b] = fileparts(cel.stack);
end

function R = emptyR()
R = struct('key','', 'label','', 'frame',[], 'x',[], 'y',[], 'q',[], 't_s',[], ...
    'trackId',[], 'spotId',[], 'tracks',{{}}, 'dt_s',NaN, 't0_s',0, ...
    'stride',1, 'offset',0, 'pages',[], 'stack','', 'base','', 'pxUm',NaN, ...
    'nFrames',0, 'nTracks',0, 'nDets',0);
end
