function [tracks, info] = dc_track(dets, linkUm, gapUm, maxGap, pxUm)
%DC_TRACK  LAP tracking of per-frame detections, for one colour.
%
%   [tracks, info] = dc_track(dets, linkUm, gapUm, maxGap, pxUm)
%
% Frame-to-frame linking is a linear assignment (matchpairs) with a cost-of-no-link above any
% feasible cost, so a within-range partner always beats a birth/death. Chains are then assembled and
% gap-closed: a track end is stitched to a later start within gapUm, across up to maxGap missing
% frames, nearest start at the smallest gap winning.
%
% DERIVED FROM SPTinMatlab's spt_track, WITH THE ER MODES REMOVED. That tracker can bias or forbid a
% link by how much of it lies off a segmented ER — a statement about that experiment's biology,
% resting on a mask this toolkit does not have. What is left is the plain LAP, which is the part that
% is about tracking rather than about ER. The two colours are tracked INDEPENDENTLY, each with its
% own frames and its own link radius; nothing here knows about the other colour, and that is
% deliberate — relating them is a separate stage with its own rules about time.
%
% INPUT
%   dets   : 1xT cell; dets{t} = Nx3 [x y quality] detections for frame t (1-based px)
%   linkUm : frame-to-frame max link distance (µm)
%   gapUm  : gap-closing max distance (µm)
%   maxGap : max missing frames to bridge
%   pxUm   : µm per pixel
% OUTPUT
%   tracks : 1xM cell; each Kx4 = [frame x y quality], sorted by frame, K>=2
%   info   : struct — nDets, nLinked

R = linkUm/pxUm; G = gapUm/pxUm;
nT = numel(dets);
d0 = R + 1;                                         % cost of leaving a spot unlinked (birth/death)
info = struct('nDets',0,'nLinked',0);
if nT > 0, info.nDets = sum(cellfun(@(d) size(d,1), dets)); end

% ---- frame-to-frame links: nxt{t}(i) = index in dets{t+1}, or 0 ----
nxt = cell(nT,1);
for t = 1:nT-1
    P = dets{t}; Q = dets{t+1};
    nxt{t} = zeros(size(P,1),1);
    if isempty(P) || isempty(Q), continue; end
    C = dc_link_cost(P, Q, R);
    Mm = matchpairs(C, d0);
    for r = 1:size(Mm,1)
        i = Mm(r,1); j = Mm(r,2);
        if isfinite(C(i,j)), nxt{t}(i) = j; end
    end
end

% ---- chain assembly: walk forward from every spot with no incoming link ----
started = cell(nT,1);
for t = 1:nT, started{t} = false(size(dets{t},1),1); end
for t = 1:nT-1
    j = nxt{t}; on = j > 0;
    started{t+1}(j(on)) = true;
end
tracks = {};
for t = 1:nT
    for i = 1:size(dets{t},1)
        if started{t}(i), continue; end
        chain = zeros(0,4); tt = t; ii = i;
        while true
            d = dets{tt}(ii,:);
            chain(end+1,:) = [tt, d(1), d(2), d(3)]; %#ok<AGROW>
            if tt < nT && nxt{tt}(ii) > 0, ii = nxt{tt}(ii); tt = tt+1; else, break; end
        end
        if size(chain,1) >= 2, tracks{end+1} = chain; end %#ok<AGROW>
    end
end

if maxGap >= 1 && numel(tracks) > 1
    tracks = gap_close(tracks, G, maxGap);
end
info.nLinked = sum(cellfun(@(c) size(c,1), tracks));
end

% =================================================================================================
function tracks = gap_close(tracks, G, maxGap)
n = numel(tracks);
endF = zeros(n,1); endXY = zeros(n,2); startF = zeros(n,1); startXY = zeros(n,2);
for k = 1:n
    tr = tracks{k};
    endF(k) = tr(end,1); endXY(k,:) = tr(end,2:3);
    startF(k) = tr(1,1); startXY(k,:) = tr(1,2:3);
end
usedStart = false(n,1); mergedInto = zeros(n,1);
[~, order] = sort(endF);
for oi = 1:n
    a = order(oi);
    te = endF(a); best = -1; bestd = inf;
    for dt = 2:maxGap+1
        for b = 1:n
            if usedStart(b) || b == a || startF(b) ~= te+dt, continue; end
            dd = hypot(endXY(a,1)-startXY(b,1), endXY(a,2)-startXY(b,2));
            if dd <= G && dd < bestd, bestd = dd; best = b; end
        end
        if best > 0, break; end                     % nearest start at the smallest gap
    end
    if best > 0, usedStart(best) = true; mergedInto(best) = a; end
end
childOf = zeros(n,1); isChild = false(n,1);
for b = 1:n
    if mergedInto(b) > 0, childOf(mergedInto(b)) = b; isChild(b) = true; end
end
out = {};
for k = 1:n
    if isChild(k), continue; end
    chain = tracks{k}; c = childOf(k);
    while c > 0, chain = [chain; tracks{c}]; c = childOf(c); end %#ok<AGROW>
    out{end+1} = chain; %#ok<AGROW>
end
tracks = out;
end
