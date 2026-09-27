function dcWriteIJ(path, labels, pages, px)
%DCWRITEIJ  A 16-bit TIFF carrying ImageJ slice labels, assembled byte by byte.
%
% Shared by the test suite, not by the toolkit. imwrite cannot write the unknown tags (50838/50839)
% that carry ImageJ's per-page labels, and a test that fed the label reader a buffer it had also
% built would prove nothing about reading a FILE — so the fixtures write real files and hand them to
% imfinfo like any other.
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

function xy = walk(start, steps, S)
% A random walk that stays in the frame. A particle that wanders off the edge stops being detected,
% which fragments its track — and a test fixture should exercise the pipeline, not the boundary.
xy = cumsum([start; steps], 1);
m = 8;                                  % keep clear of the edge by a couple of PSF widths
for d = 1:2
    v = xy(:,d);
    over = v > S-m;  v(over) = 2*(S-m) - v(over);      % reflect
    under = v < m;   v(under) = 2*m - v(under);
    xy(:,d) = min(max(v, m), S-m);
end
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

