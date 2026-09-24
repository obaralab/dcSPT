function [a, b] = dc_stem(varargin)
%DC_STEM  The one rule for putting a tracked channel into a file name, and taking it out.
%
%   stem        = dc_stem(base, key)   % compose:  'cellA' + 'ch2' -> 'cellA__ch2'
%   [base, key] = dc_stem(stem)        % split:    'cellA__ch2'    -> 'cellA', 'ch2'
%
% Tracking one stack twice — once per colour — wrote both results over the same
% <base>_tracks.xml, with no warning and no way to tell them apart afterwards. A channel token in
% the name is what keeps them separate, and this is the only place that knows its shape, because a
% name composed by one rule and split by another is worse than no rule at all.
%
% THE SEPARATOR IS A DOUBLE UNDERSCORE, and that is not decoration. Cell names in this pipeline are
% full of single underscores — HVK-3C-NoLigand-Baseline_Plate2-002_ch24_spt — so a single-underscore
% token would make '..._spt' look like a channel called 'spt' on every dataset that exists. A double
% underscore does not occur in any name the matcher produces.
%
% AN EMPTY OR DEFAULT KEY COMPOSES TO THE BARE BASE. A single-colour project keeps writing exactly
% the names it wrote before this existed, which is what makes the change safe to land: nothing has
% to be renamed, and nothing downstream sees a token until a second colour is declared.

if nargin == 2
    base = char(varargin{1}); key = char(varargin{2});
    if isempty(key) || strcmp(key, dc_stem_default())
        a = base;
    else
        a = [base '__' key];
    end
    return
end

stem = char(varargin{1});
t = regexp(stem, '^(.*)__([a-z][a-z0-9_]*)$', 'tokens', 'once');
if isempty(t)
    a = stem; b = '';
else
    a = t{1}; b = t{2};
end
end

function k = dc_stem_default()
k = '';   % the single-colour key: no token in the name
end
