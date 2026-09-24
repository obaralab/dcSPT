function A = dc_align(varargin)
%DC_ALIGN  How two colours' frames correspond — asked and answered in acquisition TIMEPOINTS.
%
%   A = dc_align(D)        the two colours of a dataset
%   A = dc_align(a, b)     or two channel entries directly
%
% "Do the two colours have the same frame numbers" is really "do their frames fall on the same
% timepoints of the acquisition". A timepoint is an integer the microscope recorded, so this is set
% comparison rather than a tolerance:
%
%   matched     both colours have a frame at every timepoint the other does — frame k of one
%               corresponds to frame k of the other, exactly
%   subsampled  one colour's timepoints are every k-th of the other's
%   partial     they overlap, but neither contains the other
%   disjoint    no timepoint is shared; nothing can be compared frame to frame
%
% Note what 'matched' does NOT mean. Two colours of an interleaved acquisition share every timepoint
% — the acquisition calls that one moment imaged twice — but the two exposures are consecutive
% PAGES, not simultaneous. The page separation is reported so that difference is visible rather than
% implied by the word.
%
% OUTPUT A
%   .relation .nFrames [nA nB] .nShared .ratio .pageGap
%   .text        one sentence for the panel

if nargin == 1
    D = varargin{1};
    assert(numel(D.channels) >= 2, 'dc_align:oneChannel', 'need two colours to align');
    a = D.channels(1); b = D.channels(2);
else
    a = varargin{1}; b = varargin{2};
end
tpA = a.tp(:); tpB = b.tp(:);
shared = intersect(tpA, tpB);
A = struct('relation','', 'nFrames',[a.nFrames b.nFrames], 'nShared',numel(shared), ...
           'ratio', a.nFrames / max(b.nFrames,1), 'pageGap',NaN, 'text','');

if isempty(shared)
    A.relation = 'disjoint';
elseif isequal(unique(tpA), unique(tpB))
    A.relation = 'matched';
else
    big = tpA; small = tpB;
    if numel(tpB) > numel(tpA), big = tpB; small = tpA; end
    if all(ismember(small, big))
        k = numel(big) / max(numel(small),1);
        if abs(k - round(k)) < 1e-9, A.relation = 'subsampled'; else, A.relation = 'partial'; end
    else
        A.relation = 'partial';
    end
end

if ~isempty(shared)
    [~, ia] = ismember(shared, tpA); [~, ib] = ismember(shared, tpB);
    A.pageGap = median(abs(b.pages(ib) - a.pages(ia)));
end
A.text = sentence(A, a, b);
end

% =================================================================================================
function s = sentence(A, a, b)
switch A.relation
    case 'matched'
        s = sprintf(['%s and %s have a frame at every one of the same %d timepoints, so their ' ...
            'frame numbers correspond exactly. The two exposures of a timepoint are %g page(s) ' ...
            'apart — consecutive, not simultaneous.'], a.key, b.key, A.nShared, A.pageGap);
    case 'subsampled'
        fast = a; slow = b;
        if b.nFrames > a.nFrames, fast = b; slow = a; end
        k = round(fast.nFrames / max(slow.nFrames,1));
        s = sprintf(['%s has a frame at only %d of %s''s %d timepoints — 1 in %d. Frame numbers do ' ...
            'not correspond: frame k of %s falls at frame %dk of %s.'], slow.key, slow.nFrames, ...
            fast.key, fast.nFrames, k, slow.key, k, fast.key);
    case 'partial'
        s = sprintf(['%s and %s share %d timepoints but neither covers the other (%d and %d ' ...
            'frames). Frame numbers do not correspond; only the shared timepoints can be compared.'], ...
            a.key, b.key, A.nShared, A.nFrames(1), A.nFrames(2));
    otherwise
        s = sprintf(['%s and %s share no timepoint at all (%d and %d frames). Nothing can be ' ...
            'compared frame to frame.'], a.key, b.key, A.nFrames(1), A.nFrames(2));
end
end
