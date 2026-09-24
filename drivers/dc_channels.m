function C = dc_channels(varargin)
%DC_CHANNELS  The colours of an acquisition: which pages each one is, and which timepoint each of
%its frames belongs to.
%
%   C = dc_channels('fromStack', path)          read it from the stack's own slice labels
%   C = dc_channels('interleaved', nPages)      no labels: assume alternating pages
%   C = dc_channels('manual', specs)            or state it outright
%   C = dc_channels('load'/'save', ...)
%
% EVERYTHING HERE IS AN INTEGER. A colour is described by which PAGES of the stack it occupies and
% which TIMEPOINT of the acquisition each of those pages is:
%
%   key        'c2' — from the acquisition's own channel number where there is one
%   label      what to call it
%   pages      nFrames x 1, the 1-based page of the stack for each of this colour's frames
%   tp         nFrames x 1, the acquisition timepoint each of those frames belongs to
%   dt_s       seconds per frame OF THIS COLOUR — carried, but used only to report a result
%   nFrames
%
% WHY PAGES AND TIMEPOINTS RATHER THAN SECONDS. Frame numbers are exact; seconds are a frame number
% multiplied by a nominal interval, and that interval is nominal — these files carry one
% `finterval` for the whole stack and no per-page timestamps at all. Multiplying an imprecise dt by
% a large frame index accumulates: the 1.1% between 0.054 and 0.0534237 is half a second across a
% 1,791-frame track. Integers also say exactly what a tracking gap is — frame k+2 following frame k
% is a gap of one, with nothing to compare against a tolerance.
%
% THE TIMEPOINT IS THE SHARED CLOCK. Two colours of one acquisition have their own frame numbering
% but the same timepoints: c:2 t:5 and c:4 t:5 are one moment of the experiment imaged twice. That
% is the acquisition's own notion, read from the file, and it is what makes a merged view exact
% rather than approximate.
%
% dt_s IS STILL PER COLOUR, because a diffusion coefficient has to come out in um^2/s eventually.
% It is applied once, at the end, to a frame count — never accumulated along an axis.

switch lower(char(varargin{1}))
    case 'fromstack',   C = fromStack(varargin{2:end});
    case 'interleaved', C = interleaved(varargin{2:end});
    case 'manual',      C = manual(varargin{2});
    case 'load',        C = loadC(varargin{2});
    case 'save',        C = saveC(varargin{2:end});
    otherwise, error('dc_channels:mode','unknown mode ''%s''', varargin{1});
end
end

% =================================================================================================
function C = fromStack(path, dtPage)
% The acquisition's own account of itself. dtPage (seconds between successive PAGES) is optional and
% only used to fill each colour's dt_s; the indices do not depend on it.
if nargin < 2, dtPage = NaN; end
L = dc_tiff_labels(path);
assert(L.ok, 'dc_channels:noLabels', ...
    ['%s carries no c:/t: slice labels, so which pages are which colour cannot be read from it. ' ...
     'Use dc_channels(''interleaved'', nPages) to assume alternating pages, or ''manual'' to state it.'], path);
C = emptyC();
for c = L.channels
    m = find(L.ch == c);                       % the pages of this colour, in order
    e = entry(sprintf('c%d', c), sprintf('channel %d', c), m(:), L.tp(m), NaN);
    if isfinite(dtPage)
        % Frames of this colour are as far apart as its pages are — read, not assumed: an
        % acquisition that dropped a page makes this uneven, and the median is the honest summary.
        e.dt_s = median(diff(m)) * dtPage;
    end
    C(end+1) = e; %#ok<AGROW>
end
assert(numel(C) >= 1, 'dc_channels:noChannels', 'the labels named no channels');
end

function C = interleaved(nPages, dtPage)
% No labels: the fallback. Alternating pages, which is what an interleaved acquisition usually is —
% but it is an ASSUMPTION, and fromStack should be preferred wherever the file will say.
if nargin < 2, dtPage = NaN; end
C = emptyC();
for i = 0:1
    pg = ((i+1) : 2 : nPages)';
    tp = (1:numel(pg))';
    e = entry(sprintf('c%d', i+1), sprintf('colour %d (%s pages)', i+1, tern(i==0,'odd','even')), pg, tp, 2*dtPage);
    C(end+1) = e; %#ok<AGROW>
end
end

function C = manual(specs)
C = emptyC();
for i = 1:numel(specs)
    s = specs(i);
    pg = s.pages(:); tp = (1:numel(pg))';
    if isfield(s,'tp') && ~isempty(s.tp), tp = s.tp(:); end
    C(end+1) = entry(s.key, gf(s,'label',s.key), pg, tp, gf(s,'dt_s',NaN)); %#ok<AGROW>
end
end

function e = entry(key, label, pages, tp, dt_s)
pages = double(pages(:)); tp = double(tp(:));
assert(numel(pages) == numel(tp), 'dc_channels:ragged', 'a colour needs one timepoint per frame');
e = struct('key',char(key), 'label',char(label), 'pages',pages, 'tp',tp, ...
           'dt_s',dt_s, 'nFrames',numel(pages));
end

function C = emptyC()
C = struct('key',{},'label',{},'pages',{},'tp',{},'dt_s',{},'nFrames',{});
end

function C = loadC(projectDir)
p = fullfile(char(projectDir),'dc_channels.json');
assert(isfile(p), 'dc_channels:noFile','no %s', p);
raw = jsondecode(fileread(p));
if isstruct(raw) && isfield(raw,'channels'), raw = raw.channels; end
if iscell(raw), raw = [raw{:}]; end
C = emptyC();
for i = 1:numel(raw)
    C(end+1) = entry(raw(i).key, gf(raw(i),'label',raw(i).key), raw(i).pages, raw(i).tp, gf(raw(i),'dt_s',NaN)); %#ok<AGROW>
end
end

function p = saveC(projectDir, C)
p = fullfile(char(projectDir),'dc_channels.json');
s = struct('channels', {arrayfun(@(e) e, C(:)', 'uni', 0)});
fid = fopen(p,'w'); assert(fid > 0, 'dc_channels:noWrite','cannot write %s', p);
fprintf(fid, '%s\n', jsonencode(s, 'PrettyPrint', true)); fclose(fid);
end

function v = gf(s,f,d), if isstruct(s)&&isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
function v = tern(c,a,b), if c, v=a; else, v=b; end, end
