function C = dc_channels(varargin)
%DC_CHANNELS  The colours of a dual-colour experiment: where each one's frames are, and when.
%
%   C = dc_channels('interleaved', dtPage)      two colours alternating in one stack
%   C = dc_channels('paired', dtA, dtB)         two colours in two stacks, their own rates
%   C = dc_channels('load', projectDir) / ('save', projectDir, C)
%
% In this toolkit a COLOUR IS NOT AN AFTERTHOUGHT, so it is described completely and in one place:
%
%   key      'a' / 'b' — the token in every file this colour produces
%   label    what to call it on screen
%   stride   page step within its stack. 2 = one colour of an interleaved pair
%   offset   0-based page offset of its first frame
%   dt_s     THIS COLOUR'S frame interval, in seconds
%   t0_s     when its first frame happened, on the experiment's clock
%   file     '' = the cell's own stack (interleaved); otherwise a suffix naming its own stack
%
% WHY t0 IS A FIELD AND NOT AN ASSUMPTION. Interleaved colours do not start together: the colour on
% the even pages begins one page interval after the odd one. Half a frame sounds like nothing until
% it is compared against the thing being measured — a molecule moving at 0.1 µm²/s covers ~60 nm in
% 10 ms, which is the same size as the colocalization distances this toolkit exists to report. The
% single-colour pipeline this came from writes T = FRAME x dt and so silently asserts t0 = 0 for
% both colours; here the origin is carried, written and testable.
%
% dt IS PER COLOUR for the same reason. Deriving it as page interval x stride assumes the only
% reason a colour has fewer frames is that it shares pages with another. A strobed colour — every
% other exposure, in its own stack — has stride 1 and a dt unrelated to the page interval. Getting
% it wrong scales every diffusion coefficient, dwell time and rate by that factor, silently.

mode = lower(char(varargin{1}));
switch mode
    case 'interleaved'
        dtPage = varargin{2};
        C = [entry('a','colour A (odd pages)',  2, 0, 2*dtPage, 0), ...
             entry('b','colour B (even pages)', 2, 1, 2*dtPage, dtPage)];
    case 'paired'
        dtA = varargin{2}; dtB = varargin{3};
        C = [entry('a','colour A', 1, 0, dtA, 0), ...
             entry('b','colour B', 1, 0, dtB, 0)];
        C(1).file = '_a'; C(2).file = '_b';
    case 'load'
        C = loadC(varargin{2});
    case 'save'
        C = saveC(varargin{2}, varargin{3});
    otherwise
        error('dc_channels:mode','unknown mode ''%s''', mode);
end
end

% =================================================================================================
function e = entry(key, label, stride, offset, dt_s, t0_s)
e = struct('key',key, 'label',label, 'stride',stride, 'offset',offset, ...
           'dt_s',dt_s, 't0_s',t0_s, 'file','');
end

function C = loadC(projectDir)
C = dc_channels('interleaved', NaN);          % the shape, with nothing decided
p = cfgPath(projectDir);
if isempty(p) || ~isfile(p), C = C([]); return; end
try
    raw = jsondecode(fileread(p));
catch ME
    error('dc_channels:badJson','%s could not be read: %s', p, ME.message);
end
if isstruct(raw) && isfield(raw,'channels'), raw = raw.channels; end
if iscell(raw), raw = [raw{:}]; end
C = C([]);
for i = 1:numel(raw)
    r = raw(i);
    e = entry(lower(strtrim(char(gf(r,'key','')))), char(gf(r,'label','')), ...
        max(1,round(double(gf(r,'stride',1)))), max(0,round(double(gf(r,'offset',0)))), ...
        double(gf(r,'dt_s',NaN)), double(gf(r,'t0_s',0)));
    e.file = char(gf(r,'file',''));
    if isempty(e.label), e.label = e.key; end
    assert(~isempty(regexp(e.key,'^[a-z][a-z0-9_]*$','once')), ...
        'dc_channels:badKey','channel key ''%s'' is not [a-z][a-z0-9_]*', e.key);
    C(end+1) = e; %#ok<AGROW>
end
assert(numel(unique({C.key})) == numel(C), 'dc_channels:dupKey','two channels share a key');
end

function p = saveC(projectDir, C)
p = cfgPath(projectDir);
assert(~isempty(p), 'dc_channels:noProject','no project folder to save into');
s = struct('channels', {arrayfun(@(e) e, C(:)', 'uni', 0)});
fid = fopen(p,'w'); assert(fid > 0, 'dc_channels:noWrite','cannot write %s', p);
fprintf(fid, '%s\n', jsonencode(s, 'PrettyPrint', true)); fclose(fid);
end

function p = cfgPath(projectDir)
p = '';
if isempty(projectDir), return; end
p = fullfile(char(projectDir), 'dc_channels.json');
end

function v = gf(s,f,d), if isstruct(s)&&isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
