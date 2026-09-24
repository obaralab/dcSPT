function t = dc_times(ch, frames)
%DC_TIMES  The acquisition clock of one colour: when each of its FRAMES happened.
%
%   t = dc_times(ch)           every frame of a processed colour (needs .nFrames)
%   t = dc_times(ch, frames)   those 0-based frame numbers
%
% A colour's clock is a property of the ACQUISITION, not of what was seen: frame 12 happened at
% t0 + 12*dt whether or not a molecule was detected in it. That distinction is the whole reason this
% is a function rather than a column on the detections. Matching two colours' DETECTION times says
% "where was the other colour's molecule when mine was seen"; matching their FRAME times says "was
% the other colour even imaged then". They are different questions and both get asked — the first
% by the partner-distance stage, the second whenever a gap has to be told apart from a blind frame.
if nargin < 2 || isempty(frames)
    n = 0; if isfield(ch,'nFrames'), n = ch.nFrames; end
    frames = (0:n-1)';
end
t = ch.t0_s + double(frames(:)) * ch.dt_s;
end
