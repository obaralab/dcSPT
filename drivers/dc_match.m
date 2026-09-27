function [cells, info] = dc_match(sptDir, opts)
%DC_MATCH  Pair each cell's TWO colour stacks in a folder, by the token that distinguishes them.
%
%   cells = dc_match(folder)
%   cells = dc_match(folder, struct('tokens',{{'Ch1','Ch2'}}))
%   cells = dc_match(folder, struct('pattern','(?<base>.*)_(?<ch>c[24])_spt'))
%   [cells, info] = dc_match(...)
%
% A cell here is a PAIR of files. Everything downstream is built on the pair being the unit, so the
% first thing a project needs is a rule for deciding which two files belong together — and that rule
% has to be visible and editable, because every lab names its files differently and a wrong pairing
% is silent: two colours of different cells produce tracks, separations and a co-motion number that
% all look completely normal.
%
% TWO WAYS TO SAY IT.
%
%   tokens   {'Ch1','Ch2'}  the simple case. The token is removed from the name and what is left is
%                           the cell key; two files with the same key and different tokens are a
%                           pair. Matching is case-insensitive.
%
%   pattern  a regexp with named groups <base> and <ch>. For names the token rule cannot express —
%            'HVK-3C-Plate1-296-011_ch24_spt.tif' against '..._ch3_spt.tif', or a channel that only
%            appears in the middle. The pattern wins when both are given.
%
% ONE FILE HOLDING BOTH COLOURS is also a pair, and is detected without either rule: a stack whose
% ImageJ slice labels name two channels is one cell on its own, and dc_channels reads the split from
% the labels. That is the case the token rule cannot see, because there is no second file to match.
%
% WHAT IT REFUSES TO GUESS. A key with one file, or with three, is reported in info.unpaired rather
% than paired with whatever is nearest. A folder where nothing pairs comes back empty with a message
% naming the tokens it looked for and the keys it found, because "0 cells" with no explanation sends
% people looking at the data when the answer is the rule.
%
% OUTPUT cells : 1xN struct
%   key        the cell identifier, the name with the token taken out
%   stacks     1x2 cellstr, the two files in token order
%   chTokens   1x2 cellstr, which token each one carried
%   interleaved  true when the pair is one file holding both colours
%   ok         true
% info : .unpaired (key + the files found) .tokens .nFiles .text

if nargin < 2 || ~isstruct(opts), opts = struct(); end
tokens  = getf(opts,'tokens', {'Ch1','Ch2'});
pattern = getf(opts,'pattern', '');
exts    = getf(opts,'ext', {'*.tif','*.tiff'});

d = [];
for e = exts(:)', d = [d; dir(fullfile(sptDir, e{1}))]; end %#ok<AGROW>
d = d(~[d.isdir]);
cells = emptyCells(); info = struct('unpaired',{{}}, 'tokens',{tokens}, 'nFiles',numel(d), 'text','');
if isempty(d)
    info.text = sprintf('no TIFF stacks in %s', sptDir);
    return
end

names = string({d.name});
paths = arrayfun(@(q) string(fullfile(q.folder, q.name)), d);
stem  = regexprep(names, '\.tiff?$', '', 'ignorecase');

% ---- a single file carrying both colours is already a pair ---------------------------------------
isDual = false(size(names));
for k = 1:numel(names)
    try
        L = dc_tiff_labels(char(paths(k)));
        isDual(k) = L.ok && numel(L.channels) > 1;
    catch
    end
end
for k = find(isDual)
    cells(end+1) = struct('key', char(stem(k)), 'stacks', {{char(paths(k)), char(paths(k))}}, ...
        'chTokens', {{'',''}}, 'interleaved', true, 'ok', true); %#ok<AGROW>
end

% ---- the rest are matched by token or pattern ------------------------------------------------------
rest = find(~isDual);
key = strings(size(rest)); tok = strings(size(rest)); got = false(size(rest));
for q = 1:numel(rest)
    k = rest(q);
    if ~isempty(pattern)
        m = regexp(char(stem(k)), pattern, 'names', 'once');
        if ~isempty(m) && isfield(m,'base') && isfield(m,'ch')
            key(q) = string(m.base); tok(q) = string(m.ch); got(q) = true;
        end
    else
        for ti = 1:numel(tokens)
            t = tokens{ti};
            if contains(stem(k), t, 'IgnoreCase', true)
                key(q) = erase(lower(stem(k)), lower(t));
                tok(q) = string(t); got(q) = true; break
            end
        end
    end
end

uk = unique(key(got));
for u = uk(:)'
    m = find(got & key == u);
    if numel(m) ~= 2
        info.unpaired{end+1} = struct('key',char(u), 'files',{cellstr(names(rest(m)))}); %#ok<AGROW>
        continue
    end
    % put them in the order the tokens were given, so colour 1 is colour 1 in every cell
    tks = tok(m);
    if ~isempty(pattern)
        [~, ord] = sort(tks);
    else
        ord = zeros(1,2);
        for ti = 1:numel(tokens)
            j = find(strcmpi(tks, tokens{ti}), 1);
            if ~isempty(j), ord(ti) = j; end
        end
        if any(ord == 0), [~, ord] = sort(tks); end
    end
    mm = m(ord);
    cells(end+1) = struct('key', char(u), ...
        'stacks', {cellstr(paths(rest(mm))')}, 'chTokens', {cellstr(tks(ord)')}, ...
        'interleaved', false, 'ok', true); %#ok<AGROW>
end

% files that matched no rule at all
lost = rest(~got);
if ~isempty(lost)
    info.unpaired{end+1} = struct('key','(no token)', 'files',{cellstr(names(lost))});
end

if isempty(cells)
    info.text = sprintf(['%d file(s) in %s, none of them a pair. Looked for the tokens %s. The keys ' ...
        'found were: %s. If the two colours are distinguished some other way, give dc_match a ' ...
        'pattern with <base> and <ch> groups.'], numel(d), sptDir, strjoin(tokens,', '), ...
        strjoin(cellstr(unique(stem))', ', '));
else
    info.text = sprintf('%d cell(s) from %d file(s)%s', numel(cells), numel(d), ...
        tern(isempty(info.unpaired), '', sprintf('; %d key(s) unpaired', numel(info.unpaired))));
end
end

% =================================================================================================
function C = emptyCells()
C = struct('key',{},'stacks',{},'chTokens',{},'interleaved',{},'ok',{});
end
function v = getf(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
function y = tern(c,a,b), if c, y=a; else, y=b; end, end
