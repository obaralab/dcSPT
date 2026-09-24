function [r, keep, wMaj] = dc_ridge(dog, px, py, ridgeMax, scale, sizeMax)
%DC_RIDGE  Principal-curvature ratio at candidate peaks — rejects RIDGE (filament) responses.
%
%   r                = dc_ridge(dog, px, py)
%   [r, keep, wMaj]  = dc_ridge(dog, px, py, ridgeMax, scale, sizeMax)
%
% THE PROBLEM THIS SOLVES. A DoG detector answers "is this a bright blob at the spot scale", and an
% extended structure — a mitochondrion bleeding through into the SPT channel, an ER tubule — is not
% one blob but a RIDGE. Along the ridge the intensity barely changes; across it, it falls off like a
% spot. Local-maximum detection therefore breaks the filament into a CHAIN of "spots" strung along
% its crest, each of them a perfectly ordinary-looking detection: right size, right brightness, and
% no threshold separates them from real single molecules, because they are equally bright.
%
% What separates them is CURVATURE, not brightness. At a point source the DoG surface curves down
% steeply in every direction, so the two eigenvalues of its Hessian are both large and similar. On a
% ridge crest it curves down steeply across the ridge and hardly at all along it, so one eigenvalue
% dwarfs the other. The ratio of the two is the discriminator, and it does not care how bright the
% structure is.
%
% This is SIFT's edge-response elimination (Lowe 2004 §4.1) applied to a 2-D DoG frame: rather than
% compute the eigenvalues, use
%
%       r = trace(H)^2 / det(H) = (R+1)^2 / R ,  R = |lambda_max| / |lambda_min|
%
% which is monotonic in R and needs no eigendecomposition. r is 4 at a perfectly round peak and
% grows without bound as the response becomes ridge-like.
%
% INPUT
%   dog      : the DoG detection image (from dc_dog)
%   px, py   : integer column/row indices of the candidate peaks (must be >=2 px from the border,
%              which dc_detect already guarantees — the Hessian is a 3x3 central difference)
%   ridgeMax : the largest curvature RATIO R to accept. [] or <=0 disables the shape test.
%   scale    : dc_dog's `scale` output — this spot's size relative to the reference. Required for
%              the SIZE test, which is expressed relative to the expected peak width.
%   sizeMax  : the largest peak width to accept, as a MULTIPLE of the width a real point source of
%              the chosen diameter produces. [] or <=0 disables the size test. 1.6 is a good value.
% OUTPUT
%   r    : Nx1 curvature ratio (R+1)^2/R at each peak; Inf where the Hessian is not a clean maximum
%   keep : Nx1 logical, false for peaks rejected as ridge-like OR too wide
%   wMaj : Nx1 major-axis peak width in px, implied by the Hessian (Inf at a saddle)
%
% TWO CRITERIA, AND THEY ARE NOT THE SAME TEST.
%   SHAPE (ridgeMax) asks "is this round or elongated" and is SIZE-BLIND — deliberately so, since
%   r = tr^2/det is a ratio of eigenvalues. It catches the crest of a filament.
%   SIZE (sizeMax) asks "is this as big as a spot of the diameter I set". It catches everything the
%   shape test cannot: a filament END, a CROSSING, a focal blob of bleedthrough — all of them round
%   enough to pass a curvature ratio, and all of them too WIDE to be a single molecule.
%
% THE STEP SIZE IS NOT A KNOB, AND MUST NOT BECOME ONE. An earlier version scaled the central
% difference with the spot, on the reasoning that a large spot needs a wider stencil. Measured, that
% changes nothing at all: enlarging the step h multiplies Dxx, Dyy and Dxy by ~h^2, so tr scales by
% h^2 and det by h^4, and r = tr^2/det is invariant — 4.001 vs 4.001 at scale 1, 4.000 vs 4.000 at
% scales 2 and 3 on the same peaks. It is a 1-px difference here because that is both the cheapest
% and the least biased choice, and the size dependence lives in the SIZE test, where it belongs.
%
% THE WIDTH. Near its peak the DoG surface is c + (1/2)*lambda*d^2 along each principal axis, so the
% distance over which it falls away is sqrt(c/|lambda|). The SMALLER |lambda| gives the WIDER axis,
% which is the one to test. For a real point source that width is proportional to `scale`: measured
% at 1.53, 3.06 and 4.58 px for scales 1, 2 and 3, i.e. W_REF*scale with W_REF = 1.53 px. The same
% filament measured 6.16, 9.07 and 11.03 px, and never dipped below 3.12, 5.29 and 8.92 — so a cut
% anywhere from about 1.3 to 1.9 times the expected width separates them at every scale.
%
% A SADDLE IS ALSO REJECTED. det(H) <= 0 means the surface curves UP along one axis — a point on a
% crest between two brighter blobs, never a spot — and r is meaningless there, so it is reported as
% Inf and rejected. That is deliberate: those saddle points are a large part of what a filament
% contributes.
%
% This test is SHAPE-ONLY and knows nothing about the mito channel. It cannot be used to ask "is
% this detection on a mitochondrion" — and must not be, since a real single molecule sitting on a
% mitochondrion is exactly the thing these experiments exist to measure. It rejects things that are
% shaped like a filament, and a molecule on a filament is still shaped like a point.

W_REF = 1.53;                                % px, the DoG peak width of a point source at scale 1
if nargin < 4, ridgeMax = []; end
if nargin < 5 || isempty(scale) || scale <= 0, scale = 1; end
if nargin < 6, sizeMax = []; end
px = px(:); py = py(:);
n = numel(px);
r = inf(n,1);
wMaj = inf(n,1);
keep = true(n,1);
if n == 0, return; end

for i = 1:n
    x = px(i); y = py(i);
    c   = dog(y,   x);
    Dxx = dog(y,   x+1) - 2*c + dog(y,   x-1);
    Dyy = dog(y+1, x)   - 2*c + dog(y-1, x);
    Dxy = (dog(y+1,x+1) - dog(y+1,x-1) - dog(y-1,x+1) + dog(y-1,x-1)) / 4;
    tr  = Dxx + Dyy;
    dt  = Dxx*Dyy - Dxy*Dxy;
    if dt > 0 && tr < 0                      % a genuine local maximum of the DoG surface
        r(i) = tr*tr / dt;
        lMin = (abs(tr) - sqrt(max(tr*tr - 4*dt, 0))) / 2;    % smaller |eigenvalue| = wider axis
        if lMin > 0 && c > 0, wMaj(i) = sqrt(c / lMin); end
    end                                      % else both stay Inf: saddle, or a minimum
end

if ~isempty(ridgeMax) && ridgeMax > 0
    % R is |lambda_max|/|lambda_min| and so is NEVER below 1. A value under 1 is not "very strict":
    % (R+1)^2/R is symmetric about R=1, so 0.1 gives exactly the cut 10 does — a very PERMISSIVE
    % one, the opposite of what anyone typing 0.1 intends. Clamp to 1 and say so rather than
    % silently honouring a number that means its own reciprocal.
    if ridgeMax < 1
        warning('dc_ridge:ratioBelowOne', ...
            ['ridgeMax = %g is below 1, but a curvature RATIO cannot be: (R+1)^2/R is symmetric ' ...
             'about 1, so %g would behave exactly like %g — far more permissive than intended. ' ...
             'Using 1 (the strictest setting: only round peaks survive).'], ...
            ridgeMax, ridgeMax, 1/ridgeMax);
        ridgeMax = 1;
    end
    rMax = (ridgeMax + 1)^2 / ridgeMax;      % the r that corresponds to that eigenvalue ratio
    keep = keep & (r <= rMax);
end
if ~isempty(sizeMax) && sizeMax > 0
    keep = keep & (wMaj <= sizeMax * W_REF * scale);
end
end
