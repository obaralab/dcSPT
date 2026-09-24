function c = dc_tiff_calib(tiffPath)
%DC_TIFF_CALIB  Pixel size, frame interval and display range from a TIFF, ImageJ-aware.
%
%   c = dc_tiff_calib(tiffPath)
%
% OUTPUT struct (every field NaN / '' when it could not be determined — never a guess):
%   .pixUm      microns per pixel
%   .dt_s       frame interval, seconds
%   .dispLo     display range low  (raw intensity units, as Fiji last saved it)
%   .dispHi     display range high
%   .nFrames    frame count from the metadata (NOT the page count)
%   .width      image width  in pixels (NaN when unreadable)
%   .height     image height in pixels — the FOV is (width-1)*pixUm, so these are what makes a FOV
%               explicable rather than one bare number on a toolbar
%   .unit       the raw 'unit=' string, for reporting
%   .src        struct of provenance per field: 'imagej' | 'resunit' | 'missing'
%
% WHY THIS EXISTS
% spt_pixel_size read only XResolution + ResolutionUnit, and returned nothing when ResolutionUnit was
% None. That is precisely how ImageJ/Fiji writes a calibrated file, so the one case the button most
% needed to handle was the one it rejected. Fiji stores the SCALE in XResolution and the UNIT in a
% plain-text block in ImageDescription:
%
%     ImageJ=1.54f
%     images=5703
%     frames=5703
%     unit=micron
%     finterval=0.010519135743379593
%     min=131.0
%     max=1904.0
%
% with ResolutionUnit = None, because the unit is not one of TIFF's two (inch/cm). XResolution is then
% PIXELS PER UNIT, so a micron file with XResolution 6.25 is 1/6.25 = 0.16 um/px.
%
% The same block carries the frame interval and the display range Fiji was showing when the file was
% saved, which is why this returns all three rather than only the pixel size.
%
% A WRONG number here is much worse than no number: it silently rescales every downstream distance,
% area and diffusion coefficient. So everything is range-checked, and anything implausible comes back
% NaN so the caller keeps whatever the user typed.

PIX_LO = 0.005; PIX_HI = 5;        % um/px — below is beyond any light microscope, above is not SPT
DT_LO  = 1e-6;  DT_HI  = 3600;     % s

c = struct('pixUm',NaN,'dt_s',NaN,'dispLo',NaN,'dispHi',NaN,'nFrames',NaN, ...
           'width',NaN,'height',NaN,'unit','', ...
           'src',struct('pixUm','missing','dt_s','missing','disp','missing'));
if nargin < 1 || isempty(tiffPath) || ~isfile(tiffPath), return; end

try, info = imfinfo(tiffPath); catch, return; end
if isempty(info), return; end
i1 = info(1);

% Dimensions come from the SAME imfinfo the calibration is read from. spt_project_calib used to call
% imfinfo a second time just for the width, which on a 5 700-page stack is not a free call.
if isfield(i1,'Width')  && ~isempty(i1.Width)  && i1.Width  > 0, c.width  = double(i1.Width);  end
if isfield(i1,'Height') && ~isempty(i1.Height) && i1.Height > 0, c.height = double(i1.Height); end

% ---- the ImageJ metadata block ------------------------------------------------------------------
ij = struct();
if isfield(i1,'ImageDescription') && ~isempty(i1.ImageDescription)
    ij = parse_ij(char(i1.ImageDescription));
end
isIJ = isfield(ij,'imagej');

if isfield(ij,'unit'), c.unit = ij.unit; end
if isfield(ij,'frames') && ij.frames > 0,      c.nFrames = ij.frames;
elseif isfield(ij,'images') && ij.images > 0,  c.nFrames = ij.images; end

% ---- pixel size ---------------------------------------------------------------------------------
% ImageJ first: it is the only source that knows the unit is microns.
xr = NaN;
if isfield(i1,'XResolution') && ~isempty(i1.XResolution), xr = double(i1.XResolution(1)); end
if isIJ && isfinite(xr) && xr > 0
    f = unit_to_um(c.unit);                       % um per ImageJ unit; NaN if not a length
    if isfinite(f)
        cand = f / xr;                            % XResolution is PIXELS PER UNIT
        if inrange(cand, PIX_LO, PIX_HI), c.pixUm = cand; c.src.pixUm = 'imagej'; end
    end
end
% Otherwise fall back to the TIFF resolution unit, which is what the old reader did.
if ~isfinite(c.pixUm) && isfinite(xr) && xr > 0
    umPerUnit = resunit_to_um(i1);
    if isfinite(umPerUnit)
        cand = umPerUnit / xr;
        if inrange(cand, PIX_LO, PIX_HI), c.pixUm = cand; c.src.pixUm = 'resunit'; end
    end
end

% ---- frame interval -----------------------------------------------------------------------------
% finterval is authoritative. fps is the reciprocal and only used when finterval is absent or zero —
% a stack saved without a time calibration carries finterval=0, which must NOT become dt=0.
if isfield(ij,'finterval') && inrange(ij.finterval, DT_LO, DT_HI)
    c.dt_s = ij.finterval; c.src.dt_s = 'imagej';
elseif isfield(ij,'fps') && ij.fps > 0 && inrange(1/ij.fps, DT_LO, DT_HI)
    c.dt_s = 1/ij.fps;     c.src.dt_s = 'imagej';
end

% ---- display range ------------------------------------------------------------------------------
% What Fiji was showing when the file was written. Not necessarily Fiji's Auto result — it is
% whatever the range happened to be — so callers should offer it, not force it.
if isfield(ij,'min') && isfield(ij,'max') && isfinite(ij.min) && isfinite(ij.max) && ij.max > ij.min
    c.dispLo = ij.min; c.dispHi = ij.max; c.src.disp = 'imagej';
end
end

% =================================================================================================
function s = parse_ij(txt)
% ImageJ's block is 'key=value' one per line. Numeric values are returned as doubles, the rest as
% char. Keys are lowercased; anything unparseable is skipped rather than erroring.
s = struct();
lines = regexp(txt, '\r\n|\n|\r', 'split');
for k = 1:numel(lines)
    t = strtrim(lines{k}); if isempty(t), continue; end
    e = find(t == '=', 1); if isempty(e), continue; end
    key = lower(strtrim(t(1:e-1)));
    val = strtrim(t(e+1:end));
    if isempty(key), continue; end
    key = regexprep(key, '[^a-z0-9_]', '');        % 'ImageJ' -> imagej; drop stray punctuation
    % A key must be a legal MATLAB identifier before it can be a field. Vendor blocks carry things
    % like '2016=...' and dynamic-field assignment ERRORS on those rather than ignoring them — and
    % this parser runs against whatever TIFF the user points at, on every cell selection.
    if isempty(key) || ~isvarname(key), continue; end
    v = str2double(val);
    if ~isnan(v), s.(key) = v; else, s.(key) = val; end
end
% 'ImageJ=1.54f' parses to NaN as a number and stays char — either way the field exists, which is
% all the caller uses it for.
if isfield(s,'imagej') && isnumeric(s.imagej) && isnan(s.imagej), s.imagej = '?'; end
end

function f = unit_to_um(u)
% Microns per ImageJ length unit. NaN for 'pixel'/'' and anything unrecognised — an uncalibrated
% stack must not be read as if it were calibrated.
f = NaN;
if isempty(u), return; end
u = lower(strtrim(char(u)));
% The micro sign reaches us as UTF-8 (0xC2 0xB5) or Greek mu (0xCE 0xBC) as often as it does as a
% single char — ImageJ rewrites a typed 'um' to U+00B5 and writes the block in the platform charset.
% Collapsing any non-ASCII run to 'u' turns every spelling into 'um'.
u = regexprep(u, '[^\x00-\x7F]+', 'u');
switch u
    case {'micron','microns','um','micrometer','micrometers','micrometre','micrometres'}, f = 1;
    case {'nm','nanometer','nanometers','nanometre','nanometres'},                        f = 1e-3;
    case {'mm','millimeter','millimeters','millimetre','millimetres'},                    f = 1e3;
    case {'cm','centimeter','centimeters','centimetre','centimetres'},                    f = 1e4;
    case {'m','meter','meters','metre','metres'},                                         f = 1e6;
    case {'inch','inches','in','"'},                                                      f = 25400;
    otherwise, f = NaN;                            % 'pixel', 'a.u.', unset, anything odd
end
end

function um = resunit_to_um(i1)
% Microns per TIFF ResolutionUnit. 'None' is genuinely undeterminable WITHOUT an ImageJ block, so it
% returns NaN here — the ImageJ path above is what rescues that case.
um = NaN;
if ~isfield(i1,'ResolutionUnit') || isempty(i1.ResolutionUnit), return; end
ru = i1.ResolutionUnit;
if ischar(ru) || isstring(ru)
    switch lower(strtrim(char(ru)))
        case {'inch','inches'},                  um = 25400;
        case {'centimeter','centimetre','cm'},   um = 10000;
        otherwise,                               um = NaN;
    end
else
    switch double(ru)
        case 2, um = 25400;      % inch
        case 3, um = 10000;      % cm
        otherwise, um = NaN;     % 1 = none
    end
end
end

function tf = inrange(v, lo, hi)
tf = isscalar(v) && isnumeric(v) && isfinite(v) && v >= lo && v <= hi;
end
