function dc_ijtiff_write(path, labels)
%DC_IJTIFF_WRITE  Write a tiny TIFF carrying ImageJ slice labels — a fixture, not a pipeline part.
%
%   dc_ijtiff_write(path, {'c:2/4 t:1/2000 - cell #1', ...})
%
% One 1x1 page per label, with tags 50839 (IJMetadata) and 50838 (IJMetadataByteCounts) laid out the
% way ImageJ lays them out: the magic 'IJIJ', a (type, count) table, then one block per item whose
% LENGTHS live in the companion tag with the header's own length first.
%
% This exists because the reader's fiddly part is that block walk, and a test that fed it a buffer it
% had also built would prove nothing about reading a FILE. imwrite cannot write unknown tags, so the
% file is assembled byte by byte here and handed to imfinfo like any other.

n = numel(labels);
ijm = uint8([]); cnt = [];
hdr = uint8([uint8('IJIJ'), uint8('labl'), be32(n)]);
cnt(1) = numel(hdr);                              % first count is the header itself
for k = 1:n
    b = utf16beBytes(labels{k});
    ijm = [ijm, b]; cnt(end+1) = numel(b); %#ok<AGROW>
end
ijm = [hdr, ijm];

% ---- data area, then the IFDs, so every offset is known before an IFD is written -----------------
out = uint8([uint8('II'), lo16(42), lo32(8)]);     % little-endian, first IFD patched below
pxOff = zeros(1,n);
for k = 1:n
    pxOff(k) = numel(out); out = [out, uint8(mod(k*40, 250))]; %#ok<AGROW>
end
if mod(numel(out),2), out = [out, uint8(0)]; end
ijmOff = numel(out); out = [out, ijm];
if mod(numel(out),2), out = [out, uint8(0)]; end
cntOff = numel(out); out = [out, lo32arr(cnt)];

ifdOff = zeros(1,n);
for k = 1:n
    ifdOff(k) = numel(out);
    E = [entry(256,3,1,1); entry(257,3,1,1); entry(258,3,1,8); entry(259,3,1,1); ...
         entry(262,3,1,1); entry(273,4,1,pxOff(k)); entry(277,3,1,1); ...
         entry(278,3,1,1); entry(279,4,1,1)];
    if k == 1
        E = [E; entry(50838,4,numel(cnt),cntOff); entry(50839,1,numel(ijm),ijmOff)]; %#ok<AGROW>
        E = sortrows(E, 1);                       % IFD entries must ascend by tag
    end
    body = uint8([]);
    for r = 1:size(E,1)
        body = [body, lo16(E(r,1)), lo16(E(r,2)), lo32(E(r,3)), lo32(E(r,4))]; %#ok<AGROW>
    end
    out = [out, lo16(size(E,1)), body, lo32(0)];   % next-IFD pointer patched below
end
for k = 1:n
    nxt = 0; if k < n, nxt = ifdOff(k+1); end
    p = ifdOff(k) + 2 + 12*tern(k==1, 11, 9);   % 0-based offset of this IFD's next-IFD pointer
    out(p+1 : p+4) = lo32(nxt);
end
out(5:8) = lo32(ifdOff(1));

fid = fopen(path,'w'); assert(fid > 0, 'cannot write %s', path);
fwrite(fid, out, 'uint8'); fclose(fid);
end

% =================================================================================================
function e = entry(tag, type, count, val), e = [tag type count val]; end
function b = lo16(v), b = uint8([mod(v,256), floor(v/256)]); end
function b = lo32(v)
v = double(v);
b = uint8([mod(v,256), mod(floor(v/256),256), mod(floor(v/65536),256), mod(floor(v/16777216),256)]);
end
function b = lo32arr(v), b = uint8([]); for i=1:numel(v), b = [b, lo32(v(i))]; end, end
function b = be32(v)
v = double(v);
b = uint8([floor(v/16777216), mod(floor(v/65536),256), mod(floor(v/256),256), mod(v,256)]);
end
function b = utf16beBytes(s)
u = double(s); b = uint8(zeros(1, 2*numel(u)));
b(1:2:end) = floor(u/256); b(2:2:end) = mod(u,256);
end
function v = tern(c,a,b), if c, v=a; else, v=b; end, end
