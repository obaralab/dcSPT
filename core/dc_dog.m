function [dog, hp, scale] = dc_dog(frame, diamUm, pxUm)
%DC_DOG  Difference-of-Gaussians spot filter, scaled by spot diameter. Port of ERAware _dog().
%
%   [dog, hp] = dc_dog(frame, diamUm, pxUm)
%
% Background-subtract with a wide Gaussian, then take a DoG at the spot scale. The three sigmas
% scale with the spot's pixel radius (diamUm/pxUm/2); at the reference (0.5 µm @ 0.10785 µm/px,
% radius ≈ 2.32 px) they equal ERAware's exact values (background σ=6, DoG σ=1.0 and 2.2 px), so
% detection matches ERAware at the default and scales sensibly for other diameters/cameras.
%
%   dog   : the detection image (local maxima above threshold are spots)
%   hp    : the high-pass (background-subtracted) image, used for the sub-pixel centroid
%   scale : this spot's size relative to the reference (1.0 at 0.5 µm / 0.10785 µm/px). Every
%           structure in `dog` is this many times bigger than at the reference, so anything that
%           MEASURES the DoG surface — dc_ridge's curvature step — has to scale by it too. Handed
%           out rather than recomputed by the caller: REF_RPX belongs in one place.
% Gaussian blur matches scipy (symmetric padding, kernel radius ceil(4σ)) for detection parity.
if nargin<2 || isempty(diamUm), diamUm = 0.5; end
if nargin<3 || isempty(pxUm),   pxUm   = 0.10785; end
f = double(frame);
REF_RPX = (0.5/0.10785)/2;                      % reference spot radius in px (≈2.318) at ERAware sigmas
rpx   = (max(diamUm,1e-3)/max(pxUm,1e-6)) / 2;  % this spot's radius in px
scale = rpx / REF_RPX;                          % 1.0 at the reference
sBg = 6.0*scale; s1 = 1.0*scale; s2 = 2.2*scale;
hp  = f - gblur(f, sBg);
dog = gblur(hp, s1) - gblur(hp, s2);
end

function y = gblur(x, s)
% A separable Gaussian: two 1-D convolutions over a symmetrically padded frame. This is exactly what
% imgaussfilt does in its SPATIAL path — same kernel (radius ceil(4s), normalised), same padding, and
% the result is bit-identical, which detection parity requires.
%
% It is here because imgaussfilt chooses its domain by kernel size, and the background blur is wide:
% at 0.4 um spots on a 0.0968 um/px camera, sigma is 5.3 px and the kernel 45 taps, for which it
% switches to the FREQUENCY domain. That costs 4.1 ms a frame against 1.3 ms spatial and 0.8 ms
% here — and detection runs three of these on every frame, so on a 5,000-frame movie the choice of
% domain alone was about four minutes a cell.
r = ceil(4*s);
v = (-r:r)';
g = exp(-(v.^2) / (2*s^2));
g = g / sum(g);
% Symmetric padding by indexing rather than padarray: [r..1, 1..n, n..n-r+1] IS symmetric padding
% (the edge pixel is mirrored, not skipped), and padarray's argument parsing costs more than the
% copy it performs — 0.6 s of an 8.5 s cell, for three calls on every frame.
[h, w] = size(x);
iy = [min(r,h-1):-1:1, 1:h, h:-1:max(1,h-r+1)];
ix = [min(r,w-1):-1:1, 1:w, w:-1:max(1,w-r+1)];
y = conv2(g, g', x(iy, ix), 'valid');
end
