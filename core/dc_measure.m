function m = dc_measure(raw, xy, rPx)
%DC_MEASURE  MEAN/MAX/TOTAL intensity over a disk on the raw frame. Port of ERAware measure_spots().
%
%   m = dc_measure(raw, xy, rPx)   ->  m = Nx3 = [MEAN, MAX, TOTAL(sum)]
%
% Disk of radius rPx around each spot (rounded centre), measured on the RAW (non-background-
% subtracted) frame. Partial disk at the borders — out-of-frame offsets are dropped (NOT clamped),
% so edge spots are not inflated.
[H, W] = size(raw);
raw = double(raw);
rr = ceil(rPx);
[ox, oy] = meshgrid(-rr:rr, -rr:rr);
in = (ox.^2 + oy.^2) <= rPx^2;               % disk offsets
ox = ox(in); oy = oy(in);
N = size(xy,1); m = zeros(N,3);
for k = 1:N
    cx = round(xy(k,1)); cy = round(xy(k,2));
    cc = cx + ox; rrr = cy + oy;
    ok = cc >= 1 & cc <= W & rrr >= 1 & rrr <= H;
    v = raw(sub2ind([H W], rrr(ok), cc(ok)));
    if isempty(v), m(k,:) = [0 0 0]; else, m(k,:) = [mean(v) max(v) sum(v)]; end
end
end
