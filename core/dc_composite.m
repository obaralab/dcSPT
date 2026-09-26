function [rgb, info] = dc_composite(imA, imB, opts)
%DC_COMPOSITE  Two raw frames as one colour image, the way a merged channel is meant to be read.
%
%   rgb = dc_composite(imA, imB)
%   [rgb, info] = dc_composite(imA, imB, struct('mode','composite','loA',[],'hiA',[]))
%
% MAGENTA AND GREEN, NOT RED AND GREEN. The overlap of magenta and green is white and both channels
% stay visible to a red-green colour-blind reader, which red/green does not: to a deuteranope a
% red/green merge is two shades of the same muddy colour and the whole point of the picture is gone.
% Magenta also reads as "channel A" at a glance because it is the brighter of the two on a dark
% ground.
%
% EACH CHANNEL IS SCALED ON ITS OWN. A dim channel beside a bright one would be invisible under a
% shared display range, and the merge would look like the bright channel alone. Percentile limits per
% channel, so "how much of each is here" is readable even when the two differ by an order of
% magnitude in photons. That makes the composite a picture of WHERE things are, not of how bright
% they are — read intensity off the traces in the pair panel, never off this.
%
% OPTIONS
%   .mode   'composite' (default) | 'a' | 'b'   — one channel alone still goes through the same
%           scaling, so stepping between views does not also change the contrast
%   .loA .hiA .loB .hiB   display limits; empty means the percentiles below
%   .pct    [50 99.9] the percentiles used when a limit is not given. The LOW one is the median, not
%           the 1st percentile: in a sparse single-molecule frame almost every pixel is background,
%           so a 1st-percentile floor maps that background to mid-grey and the picture is a wall of
%           speckle with the molecules barely above it. Putting the floor at the median sends half
%           the pixels to black and leaves the range for what is actually bright.
%   .gamma  1  applied after scaling, for pulling dim spots up without touching the limits
%
% OUTPUT
%   rgb   H x W x 3 double in [0 1]
%   info  the limits actually used, so a caller can show them or reuse them across frames

if nargin < 3 || ~isstruct(opts), opts = struct(); end
mode = lower(char(getf(opts,'mode','composite')));
pct  = getf(opts,'pct',[50 99.9]);
gam  = getf(opts,'gamma',1);

A = double(imA); B = double(imB);
assert(isequal(size(A), size(B)), 'dc_composite:size', ...
    'the two frames must be the same size (%s and %s)', mat2str(size(A)), mat2str(size(B)));

[loA, hiA] = limits(A, getf(opts,'loA',[]), getf(opts,'hiA',[]), pct);
[loB, hiB] = limits(B, getf(opts,'loB',[]), getf(opts,'hiB',[]), pct);
a = scale(A, loA, hiA, gam);
b = scale(B, loB, hiB, gam);

switch mode
    case 'a',  rgb = cat(3, a, zeros(size(a)), a);      % magenta alone
    case 'b',  rgb = cat(3, zeros(size(b)), b, zeros(size(b)));   % green alone
    otherwise, rgb = cat(3, a, b, a);                   % magenta + green, overlap -> white
end
info = struct('loA',loA,'hiA',hiA,'loB',loB,'hiB',hiB,'mode',mode, ...
    'text', sprintf('magenta %s, green %s (%.4g–%.4g / %.4g–%.4g)', ...
                    'A', 'B', loA, hiA, loB, hiB));
end

function [lo, hi] = limits(X, lo, hi, pct)
if isempty(lo), lo = prctile(X(:), pct(1)); end
if isempty(hi), hi = prctile(X(:), pct(2)); end
if ~(hi > lo), hi = lo + max(1, abs(lo)*0.01); end
end

function y = scale(X, lo, hi, gam)
y = (X - lo) / (hi - lo);
y = min(max(y, 0), 1);
if gam ~= 1, y = y .^ gam; end
end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
