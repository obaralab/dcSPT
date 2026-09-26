function [S, info] = dc_steps(D, opts)
%DC_STEPS  Every track's step vectors, keyed by the TIMEPOINT the step starts at.
%
%   S = dc_steps(D)
%   S = dc_steps(D, struct('span', 1, 'minLen', 3))
%
% A step is the displacement between two consecutive localizations of one track. Two steps from
% different tracks are SIMULTANEOUS when they start at the same acquisition timepoint and cover the
% same number of timepoints — which is why both are recorded and why the default keeps only span 1.
%
% WHY SPAN MATTERS. A track with a gap has a localization at tp k and the next at tp k+2; the
% displacement between them covers twice as long as its neighbours' steps and is on average sqrt(2)
% longer for the same diffusion. Correlating it against a one-timepoint step compares two different
% quantities, and the longer step also overlaps two of the other track's steps, so it is not even a
% pairing. Mixing spans inflates the apparent step size of whichever track gapped more, which biases a
% dot product (um^2) while leaving a cosine alone. The default drops them; 'span', [] keeps all and
% records the span so a caller can decide.
%
% GAP CLOSING AND CO-MOTION PULL AGAINST EACH OTHER, and `info` is how you see it. A tracker allowed
% to close gaps produces steps that span 2, 3 or more timepoints, and every one of those is dropped
% here. Raise maxGap to stop tracks fragmenting and you can lose most of the step population the
% pairing needs — on a 26-track fixture, maxGap 3 left 285 pairable step pairs where maxGap 1 leaves
% far more. `info.kept` and `info.bySpan` report the trade so it is a decision rather than a surprise.
%
% OUTPUT S — one row per step, a table:
%   .trackId   global track id (unique across the cell; see dc_dataset)
%   .ch        the colour, categorical
%   .tp0 .tp1  the timepoints the step runs between
%   .span      tp1 - tp0, in timepoints
%   .x .y      the step's START position, um
%   .ux .uy    the step vector, um
%   .len       |u|, um
%
% Steps are what a co-motion measurement is made of, so this is deliberately a flat table rather than
% a per-track cell: the pairing in dc_comotion is a join on .tp0, and a table makes that one line.

if nargin < 2 || ~isstruct(opts), opts = struct(); end
span   = getf(opts, 'span', 1);            % [] keeps every span
minLen = getf(opts, 'minLen', 2);          % a track needs this many localizations to have a step

P = D.spots;
P = P(isfinite(P.trackId), :);
assert(~isempty(P), 'dc_steps:noTracks', 'no tracked spots in this dataset');
P = sortrows(P, {'trackId', 'tp'});

tid = P.trackId;
same = [false; diff(tid) == 0];            % this row continues the previous row's track

tp0 = P.tp(1:end-1);   tp1 = P.tp(2:end);
x0  = P.x(1:end-1);    y0  = P.y(1:end-1);
ux  = P.x(2:end) - x0; uy  = P.y(2:end) - y0;
keep = same(2:end);                        % only pairs of rows inside one track

S = table(tid(1:end-1), P.ch(1:end-1), tp0, tp1, tp1 - tp0, x0, y0, ux, uy, hypot(ux,uy), ...
    'VariableNames', {'trackId','ch','tp0','tp1','span','x','y','ux','uy','len'});
S = S(keep, :);

nAll = height(S);
spans = S.span;
if ~isempty(span)
    S = S(S.span == span, :);
end
info = struct('nAll',nAll, 'nKept',height(S), 'span',span, ...
              'bySpan', [unique(spans), accumarray(findgroups(spans), 1)], 'text','');
if ~isempty(span) && nAll > 0
    info.text = sprintf(['%d of %d steps span exactly %d timepoint(s) and are usable; %d were gap-' ...
        'closed over more and cannot be paired with a single-timepoint step'], ...
        info.nKept, nAll, span, nAll - info.nKept);
else
    info.text = sprintf('%d steps, every span kept', nAll);
end
if minLen > 2
    n = groupcounts(S, 'trackId');
    ok = n.trackId(n.GroupCount >= minLen - 1);
    S = S(ismember(S.trackId, ok), :);
end
S = sortrows(S, {'tp0','trackId'});
end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f), v = s.(f); else, v = d; end
end
