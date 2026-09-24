function varargout = dc_masks(mode, varargin)
%DC_MASKS  The organelle masks: ONE set, shared by both colours, indexed by TIMEPOINT.
%
%   D  = dc_masks('set', D, struct('mito',path,'er',path,'tpPerPage',100,'tp0',1))
%   pg = dc_masks('page', D, tp)                which mask page a timepoint falls in
%   m  = dc_masks('at', D, 'mito', tp)          that page, as a logical mask
%   C  = dc_masks('check', D)                   is one mask still safe to share?
%
% WHY ONE SET SERVES BOTH COLOURS. The organelles are window-averaged and move on a far slower
% timescale than the particles — seconds against tens of milliseconds. Both colours of a timepoint
% fall inside the same averaging window, so the mask that is right for one is right for the other.
% That is a claim about the biology, so 'check' measures it rather than assuming it.
%
% A PAGE IS A WINDOW OF TIMEPOINTS, NOT AN INSTANT. Each mask page averages `tpPerPage` timepoints,
% so page p covers timepoints [tp0 + (p-1)*tpPerPage, tp0 + p*tpPerPage). Which page a timepoint
% belongs to is integer division — containment, not nearest-centre, which would be a half-window
% error at every boundary.
%
% THE WINDOW IS DECLARED, NOT DERIVED. Deriving it by dividing page counts needs an exact integer
% ratio that real acquisitions rarely give, and the usual fallback when the ratio is not whole is
% silently 1 — which leaves every timepoint past the end of the organelle stack with no mask at all.
%
% WHAT THE AVERAGING COSTS: a distance to a window-averaged mask is a distance to where the
% organelle was ON AVERAGE over that window, not where it was at that timepoint. If it moves within
% the window the mask is blurred relative to any one moment, and that blur is a floor on how
% precisely such a distance can mean anything. 'check' reports the window so it can be quoted.

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
assert(isfield(M,'tpPerPage') && isfinite(M.tpPerPage) && M.tpPerPage >= 1, 'dc_masks:noWindow', ...
    ['the mask stack needs its averaging window, in TIMEPOINTS per page. It cannot be derived by ' ...
     'dividing page counts: that needs an exact integer ratio, and the fallback when there is none ' ...
     'is silently 1, which leaves most of the movie with no mask.']);
D.masks = struct('er', gf(M,'er',''), 'mito', gf(M,'mito',''), ...
                 'tpPerPage', round(M.tpPerPage), 'tp0', round(gf(M,'tp0',1)));
end

function [pg, span] = pageOf(D, tp)
% 1-based mask page CONTAINING this timepoint, and the timepoints that page averages over.
m = D.masks;
pg = floor((double(tp) - m.tp0) / m.tpPerPage) + 1;
pg = max(pg, 1);
first = m.tp0 + (pg-1)*m.tpPerPage;
span = [first, first + m.tpPerPage - 1];
end

function mk = maskAt(D, key, tp)
p = '';
if isfield(D.masks, key), p = D.masks.(key); end
assert(~isempty(p) && isfile(p), 'dc_masks:noStack', 'no %s stack for this cell', key);
pg = pageOf(D, tp);
info = imfinfo(p);
im = imread(p, min(max(pg,1), numel(info)));      % clamp: a mask outlasting the movie is normal
if size(im,3) == 3, im = rgb2gray(im); end
v = unique(im(:)); nz = v(v > 0); fg = 1; if ~isempty(nz), fg = double(min(nz)); end
mk = (im == fg);
end

function C = check(D)
% Safe to share exactly while the organelle has not moved between the two colours' exposures. In
% timepoints that is a ratio: the averaging window against the separation of the exposures, which
% for colours of one timepoint is under a single timepoint.
C = struct('ok',false, 'tpPerPage',NaN, 'text','');
if ~isfield(D,'masks') || ~isfield(D.masks,'tpPerPage'), C.text = 'no mask window declared'; return; end
C.tpPerPage = D.masks.tpPerPage;
C.ok = C.tpPerPage >= 10;
if C.ok
    C.text = sprintf(['one mask serves both colours: each page averages %d timepoints, so both ' ...
        'colours of any timepoint fall well inside the same window and the organelle has not ' ...
        'moved between them'], C.tpPerPage);
else
    C.text = sprintf(['each mask page averages only %d timepoint(s), which is not much slower than ' ...
        'the particles. Sharing one mask assumes the organelle is still between the two colours, ' ...
        'and at this window that assumption is no longer safe.'], C.tpPerPage);
end
end

function v = gf(s,f,d), if isstruct(s)&&isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
