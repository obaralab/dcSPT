function varargout = dc_masks(mode, varargin)
%DC_MASKS  The organelle masks: ONE set, shared by both colours, indexed by TIME.
%
%   D    = dc_masks('set', D, struct('er',pathE,'mito',pathM,'dt_s',dtPage,'t0_s',0))
%   [pg, t] = dc_masks('page', D, t_s)            which page of the mask stack that instant is
%   m    = dc_masks('at', D, 'mito', t_s)         that page, as a logical mask (cached)
%   C    = dc_masks('check', D)                   does the shared-mask assumption still hold?
%
% WHY ONE SET SERVES BOTH COLOURS. The ER and mitochondria move on a far slower timescale than the
% particles: seconds to minutes against 27 ms frames. So the mask that is right for one colour at an
% instant is right for the other, even though the two colours were not imaged at exactly that instant
% — the organelle did not go anywhere in the 13 ms between them. That is a statement about the
% biology, so it is recorded here and CHECKED rather than assumed ('check' below).
%
% WHY BY TIME AND NOT BY FRAME. A mask page belongs to a moment of the acquisition; each colour
% samples the acquisition differently, so `page = frame + 1` gives two different answers for one
% instant and neither is the mask's own numbering. Time is the only index both colours share. (The
% single-colour pipeline still has this bug on de-interleaved runs, where every viewer is off by the
% stride.)
%
% WHY THE PAGE RATE IS DECLARED AND NOT DERIVED. Deriving it from page counts needs an exact integer
% ratio between the stacks, which real acquisitions do not provide — the first attempt here carried
% a 1.1% rounding — and when the ratio is not whole the fallback is silently 1, which leaves every
% frame past the end of the organelle stack with no mask at all.

switch lower(char(mode))
    case 'set',   varargout{1} = setM(varargin{:});
    case 'page',  [varargout{1}, varargout{2}] = pageOf(varargin{:});
    case 'at',    varargout{1} = maskAt(varargin{:});
    case 'check', varargout{1} = check(varargin{:});
    otherwise,    error('dc_masks:mode','unknown mode ''%s''', mode);
end
end

% =================================================================================================
function D = setM(D, M)
assert(isfield(M,'dt_s') && isfinite(M.dt_s) && M.dt_s > 0, 'dc_masks:noRate', ...
    ['the mask stack needs its own page interval in seconds. It cannot be derived from page counts: ' ...
     'that needs an exact integer ratio, and a real acquisition rarely gives one.']);
D.masks = struct('er', gf(M,'er',''), 'mito', gf(M,'mito',''), ...
                 'dt_s', M.dt_s, 't0_s', gf(M,'t0_s',0));
end

function [pg, tPage] = pageOf(D, t_s)
% 0-based page index in the mask stack, and the instant that page actually is.
m = D.masks;
pg = max(0, round((t_s - m.t0_s) / m.dt_s));
tPage = m.t0_s + pg * m.dt_s;
end

function mk = maskAt(D, key, t_s)
p = '';
if isfield(D.masks, key), p = D.masks.(key); end
assert(~isempty(p) && isfile(p), 'dc_masks:noStack', 'no %s stack for this cell', key);
[pg, ~] = pageOf(D, t_s);
info = imfinfo(p);
idx = min(max(pg + 1, 1), numel(info));          % clamp: a mask outlasting the movie is normal
im = imread(p, idx);
if size(im,3) == 3, im = rgb2gray(im); end
v = unique(im(:)); nz = v(v > 0); fg = 1; if ~isempty(nz), fg = double(min(nz)); end
mk = (im == fg);
end

function C = check(D)
% Is sharing one mask between the two colours still safe? It is safe exactly while the organelle
% does not move appreciably in the time between the colours. That is a ratio, so it can be measured.
C = struct('ok',true, 'lag_s',0, 'maskDt_s',NaN, 'ratio',NaN, 'text','');
if ~isfield(D,'masks') || ~isfield(D.masks,'dt_s') || ~isfinite(D.masks.dt_s)
    C.ok = false; C.text = 'no mask page interval declared'; return;
end
C.maskDt_s = D.masks.dt_s;
if numel(D.channels) >= 2
    C.lag_s = abs(D.channels(2).t0_s - D.channels(1).t0_s);
end
% The worst case is not the lag between the colours but the worst gap between any colour's frame and
% the mask page it will be given.
C.ratio = C.maskDt_s / max(max(arrayfun(@(c) c.dt_s, D.channels)), eps);
C.ok = C.maskDt_s >= 10 * max([C.lag_s, arrayfun(@(c) c.dt_s, D.channels)]);
if C.ok
    C.text = sprintf(['one mask serves both colours: it changes every %.3g s, %.0fx slower than the ' ...
        'frames and %.0fx the %.1f ms between the colours, so the organelle has not moved between them'], ...
        C.maskDt_s, C.ratio, C.maskDt_s/max(C.lag_s,eps), 1000*C.lag_s);
else
    C.text = sprintf(['the mask changes every %.3g s, which is NOT much slower than the frames ' ...
        '(%.3g s) or the %.1f ms between the colours. Sharing one mask assumes the organelle is ' ...
        'still between them, and at this rate that assumption is no longer safe.'], ...
        C.maskDt_s, max(arrayfun(@(c) c.dt_s, D.channels)), 1000*C.lag_s);
end
end

function v = gf(s,f,d), if isstruct(s)&&isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
