function A = dc_align(D, opts)
%DC_ALIGN  How the two colours' clocks relate — which is what "do they have the same frame numbers"
%really asks.
%
%   A = dc_align(D)              D from dc_dataset
%   A = dc_align(chA, chB)       or two registry entries directly
%
% FOUR ANSWERS, and only the first makes frame k of one colour the same instant as frame k of the
% other. The merged view is driven by TIME precisely because the other three exist and look exactly
% like the first if you go by index.
%
%   simultaneous  same dt, same t0            frames correspond 1:1
%   interleaved   same dt, t0 differs by <dt  SAME COUNT, never simultaneous — a constant lag
%   subsampled    one dt is ~k x the other    only 1 in k frames has a counterpart
%   unequal       anything else               nearest partner within a stated tolerance
%
% 'interleaved' is the one worth labouring. The two colours have the same number of frames, so every
% index lines up and a panel stepping by frame looks perfectly sensible — while showing two things
% that were never imaged at the same moment. On the data here that lag is 13.4 ms, during which a
% molecule at 0.1 um^2/s moves ~70 nm: the size of the colocalization distances being measured.
%
% OUTPUT A
%   .relation .sameCount .nFrames [nA nB] .dt_s [..] .t0_s [..] .ratio .lag_s
%   .maxGap_s   worst |time gap| when every frame of the finer colour takes its nearest partner
%   .matched    fraction of the finer colour's frames with a partner within .tol_s
%   .tol_s      what counts as simultaneous here — half the coarser colour's frame interval
%   .text       one sentence, for the panel to print

if nargin >= 2 && isstruct(opts) && isfield(opts,'dt_s')
    chA = D; chB = opts; nA = fr(chA); nB = fr(chB);
else
    assert(numel(D.channels) >= 2, 'dc_align:oneChannel', 'need two colours to align');
    chA = D.channels(1); chB = D.channels(2);
    nA = frOf(D, chA.key); nB = frOf(D, chB.key);
end

dtA = chA.dt_s; dtB = chB.dt_s;
lag = chB.t0_s - chA.t0_s;
ratio = dtB / dtA;
tol = max(dtA, dtB) / 2;

sameDt = abs(dtB - dtA) <= 1e-9 * max(dtA, dtB);
kRatio = round(max(ratio, 1/ratio));
nearInt = abs(max(ratio, 1/ratio) - kRatio) <= 0.02 * kRatio;   % 2%: a rounded dt should still count

if sameDt && abs(lag) <= 1e-9
    rel = 'simultaneous';
elseif sameDt
    rel = 'interleaved';
elseif nearInt && kRatio > 1
    rel = 'subsampled';
else
    rel = 'unequal';
end

% Every frame of the FINER colour, taking its nearest partner: the worst gap and how many find one.
if dtA <= dtB, tFine = chA.t0_s + (0:nA-1)'*dtA; tCoarse = chB.t0_s + (0:nB-1)'*dtB;
else,          tFine = chB.t0_s + (0:nB-1)'*dtB; tCoarse = chA.t0_s + (0:nA-1)'*dtA; end
M = dc_time_match(tFine, tCoarse, struct('tol_s', tol));
maxGap = 0; if any(M.matched), maxGap = max(abs(M.dt_s(M.matched))); end

A = struct('relation',rel, 'sameCount', nA == nB, 'nFrames',[nA nB], ...
           'dt_s',[dtA dtB], 't0_s',[chA.t0_s chB.t0_s], 'ratio',ratio, 'lag_s',lag, ...
           'maxGap_s',maxGap, 'matched', mean(M.matched), 'tol_s',tol, 'text','');
A.text = sentence(A, chA, chB);
end

% =================================================================================================
function s = sentence(A, chA, chB)
a = chA.key; b = chB.key;
switch A.relation
    case 'simultaneous'
        s = sprintf(['%s and %s share a clock: %d frames each, %.4g s apart, imaged at the same ' ...
            'instants. Frame numbers correspond.'], a, b, A.nFrames(1), A.dt_s(1));
    case 'interleaved'
        s = sprintf(['%s and %s have the SAME number of frames (%d) and the same %.4g s interval, ' ...
            'but %s is offset by %.1f ms and so is never imaged at the same instant. Frame %d of ' ...
            'one is not frame %d of the other.'], a, b, A.nFrames(1), A.dt_s(1), b, 1000*A.lag_s, 7, 7);
    case 'subsampled'
        k = round(max(A.ratio, 1/A.ratio));
        slow = b; if A.dt_s(1) > A.dt_s(2), slow = a; end
        s = sprintf(['%s is imaged %gx more slowly (%d frames against %d, %.4g s against %.4g s), ' ...
            'so only 1 in %g frames of the faster colour has a partner imaged at the same time; ' ...
            'the rest are up to %.1f ms away.'], slow, k, min(A.nFrames), max(A.nFrames), ...
            max(A.dt_s), min(A.dt_s), k, 1000*A.maxGap_s);
    otherwise
        s = sprintf(['%s and %s run at unrelated rates (%.4g s and %.4g s, %d and %d frames). ' ...
            'Frame numbers mean nothing across them: %.0f%% of the faster colour''s frames have a ' ...
            'partner within %.1f ms, the worst gap being %.1f ms.'], a, b, A.dt_s(1), A.dt_s(2), ...
            A.nFrames(1), A.nFrames(2), 100*A.matched, 1000*A.tol_s, 1000*A.maxGap_s);
end
if ~A.sameCount && ~strcmp(A.relation,'subsampled') && ~strcmp(A.relation,'unequal')
    s = [s sprintf(' (Frame counts differ: %d against %d.)', A.nFrames(1), A.nFrames(2))];
end
end

function n = fr(ch), n = 0; if isfield(ch,'nFrames'), n = ch.nFrames; end, end

function n = frOf(D, key)
n = 0;
if isempty(D.tracks), return; end
k = find(strcmp({D.tracks.key}, key), 1);
if ~isempty(k), n = D.tracks(k).nFrames; end
end
