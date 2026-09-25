function [S, Dr] = dc_drift(S, opts)
%DC_DRIFT  The motion the whole field shares, estimated per timepoint and taken out of the steps.
%
%   [S, Dr] = dc_drift(S)                 estimate and subtract
%   Dr      = dc_drift(S, struct('apply', false))    estimate only
%
% WHY THIS COMES FIRST. A stage that moves, a sample that settles, or a focus that wanders adds the
% SAME vector to every step in that frame. Two molecules on opposite sides of the cell then have
% correlated steps for a reason that has nothing to do with them. Drift is the dominant false positive
% in any co-motion measurement, and unlike localization noise it does not average away with more
% steps — it biases the mean, so collecting more data makes a drift artefact MORE significant, not
% less.
%
% HOW. At each timepoint the drift is the MEDIAN step vector over every track with a step there.
% Median, not mean: a handful of genuinely co-moving molecules, or one bright aggregate, would drag a
% mean. A timepoint with fewer than minN steps gets no estimate and is left alone rather than being
% handed a number made from two tracks.
%
% WHAT IT CANNOT DO. This removes only what is common to the WHOLE field. Motion shared by a local
% group — the thing this analysis is looking for — survives it, which is the point. But so does any
% large-scale flow that varies across the field, so a strong gradient would still show up as
% correlation at long range. dc_comotion reports the far-field level for exactly that reason: read it
% as the floor, and if subtracting drift does not bring it near zero, the residual is structured and
% the far-field null is doing the work instead.
%
% OUTPUT
%   S   the same step table with .ux .uy corrected (and .len recomputed), plus .dx .dy, what was taken
%       out, so the correction is visible rather than silently folded in
%   Dr  .tp .dx .dy .n per timepoint, .medianSpeed, .totalUm, .text

if nargin < 2 || ~isstruct(opts), opts = struct(); end
minN  = getf(opts, 'minN', 5);
apply = getf(opts, 'apply', true);

[g, tp] = findgroups(S.tp0);
n  = splitapply(@numel, S.ux, g);
dx = splitapply(@median, S.ux, g);
dy = splitapply(@median, S.uy, g);
weak = n < minN;
dx(weak) = 0; dy(weak) = 0;               % not enough tracks to say anything: correct nothing

Dr = struct('tp', tp, 'dx', dx, 'dy', dy, 'n', n, ...
            'medianSpeed', median(hypot(dx,dy)), ...
            'totalUm', hypot(sum(dx), sum(dy)), ...
            'nEstimated', nnz(~weak), 'nTimepoints', numel(tp), 'text', '');
Dr.text = sprintf(['drift: median %.1f nm per step over %d of %d timepoints, %.3f um net across the ' ...
    'movie (timepoints with fewer than %d steps were left uncorrected)'], ...
    1000*Dr.medianSpeed, Dr.nEstimated, Dr.nTimepoints, Dr.totalUm, minN);

if nargout == 1 && ~apply, S = Dr; return; end
if apply
    S.dx = dx(g); S.dy = dy(g);
    S.ux = S.ux - S.dx;  S.uy = S.uy - S.dy;
    S.len = hypot(S.ux, S.uy);
end
end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f), v = s.(f); else, v = d; end
end
