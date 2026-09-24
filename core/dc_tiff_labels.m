function L = dc_tiff_labels(tiffPath)
%DC_TIFF_LABELS  What the acquisition itself says each page of a stack is.
%
%   L = dc_tiff_labels(path)
%
% ImageJ keeps a label per page in TIFF tag 50839, and an acquisition that wrote a hyperstack leaves
% its own indexing in them:
%
%       c:2/4 t:1/2000 - HVK-3C-Plate1-296-011 #1
%       c:4/4 t:1/2000 - HVK-3C-Plate1-296-011 #1
%       c:2/4 t:2/2000 - ...
%
% CHANNEL AND TIMEPOINT, PER PAGE, AS INTEGERS, written by the microscope rather than inferred from
% page parity. That is worth more than it looks. Deducing "odd pages are colour A" from a stride
% assumes the acquisition never dropped a frame, never started on the other channel, and never
% changed order — assumptions that are invisible when they break, because the pages still alternate
% plausibly. Here the file is asked instead.
%
% There are NO per-page timestamps in these files: every page carries the same ImageJ header with a
% single nominal `finterval`. The timepoint index is what exists, and it is exact, which is why the
% rest of this toolkit indexes by it.
%
% OUTPUT L
%   .ok        true when labels were found and parsed
%   .source    'labels' | 'none'
%   .nPages    pages in the file
%   .ch        nPages x 1, the channel index of each page (NaN where unlabelled)
%   .tp        nPages x 1, the timepoint index of each page (NaN where unlabelled)
%   .nCh .nTp  what the labels say the acquisition HAD (c:2/4 -> nCh 4)
%   .channels  the distinct channel indices present, in page order of first appearance
%   .text      one line describing what was found
%   .labels    the raw label strings, for anything this does not parse
%
% Absent or unreadable labels are not an error: .ok is false and the caller falls back to whatever
% it knew before. A plain single-channel stack usually has none, and that is fine.

L = struct('ok',false, 'source','none', 'nPages',0, 'ch',[], 'tp',[], ...
           'nCh',NaN, 'nTp',NaN, 'channels',[], 'text','no slice labels in this stack', 'labels',{{}});
try
    I = imfinfo(tiffPath);
catch ME
    L.text = ['could not read the stack: ' ME.message]; return
end
L.nPages = numel(I);
if ~isfield(I(1),'UnknownTags') || isempty(I(1).UnknownTags), return; end

u = I(1).UnknownTags;
ids = [u.ID];
kMeta = find(ids == 50839, 1);      % IJMetadata
kCnt  = find(ids == 50838, 1);      % IJMetadataByteCounts
if isempty(kMeta) || isempty(kCnt), return; end

labels = decodeLabels(uint8(u(kMeta).Value(:))', double(u(kCnt).Value(:))');
if isempty(labels), return; end
L.labels = labels;

ch = nan(L.nPages,1); tp = nan(L.nPages,1); nCh = NaN; nTp = NaN;
for k = 1:min(numel(labels), L.nPages)
    t = regexp(labels{k}, 'c:(\d+)/(\d+)\s+t:(\d+)/(\d+)', 'tokens', 'once');
    if isempty(t), continue; end
    ch(k) = str2double(t{1}); nCh = str2double(t{2});
    tp(k) = str2double(t{3}); nTp = str2double(t{4});
end
if ~any(isfinite(ch))
    L.text = sprintf('%d slice labels, none carrying c:/t: indices', numel(labels));
    return
end

L.ok = true; L.source = 'labels';
L.ch = ch; L.tp = tp; L.nCh = nCh; L.nTp = nTp;
L.channels = unique(ch(isfinite(ch)), 'stable')';
parts = arrayfun(@(c) sprintf('c:%d (%d pages)', c, nnz(ch == c)), L.channels, 'uni', 0);
L.text = sprintf('%d pages: %s, timepoints %d-%d of %d', L.nPages, strjoin(parts, ', '), ...
    min(tp(isfinite(tp))), max(tp(isfinite(tp))), nTp);
end

% =================================================================================================
function labels = decodeLabels(b, cnt)
% ImageJ's metadata blob: 'IJIJ', then (type, count) pairs, then one block per counted item. The
% block LENGTHS are in the companion tag, first entry being the header itself — which is why both
% tags are needed and why walking the blob by eye does not work.
labels = {};
if numel(b) < 8 || ~strcmp(char(b(1:4)), 'IJIJ') || isempty(cnt), return; end
hdrLen = cnt(1);
if hdrLen < 12 || hdrLen > numel(b), return; end
nTypes = floor((hdrLen - 4) / 8);
types = cell(1,nTypes); counts = zeros(1,nTypes);
for i = 1:nTypes
    o = 4 + (i-1)*8;
    types{i} = char(b(o+1:o+4));
    counts(i) = be32(b(o+5:o+8));
end
pos = hdrLen; blk = 2;                       % blocks start after the header; cnt(1) was the header
for i = 1:nTypes
    for j = 1:counts(i)
        if blk > numel(cnt), return; end
        len = cnt(blk); blk = blk + 1;
        if pos + len > numel(b), return; end
        chunk = b(pos+1 : pos+len); pos = pos + len;
        if strcmp(types{i}, 'labl')
            labels{end+1} = utf16be(chunk); %#ok<AGROW>
        end
    end
end
end

function v = be32(q)
v = double(q(1))*2^24 + double(q(2))*2^16 + double(q(3))*2^8 + double(q(4));
end

function s = utf16be(q)
if mod(numel(q),2) == 1, q = q(1:end-1); end
hi = double(q(1:2:end)); lo = double(q(2:2:end));
s = char(hi*256 + lo);
end
