function info = dc_comotion_view(ax, st, mode, opts)
%DC_COMOTION_VIEW  A pair's motion with the space between them taken out, three ways.
%
%   info = dc_comotion_view(ax, st, 'centred')
%   info = dc_comotion_view(ax, st, 'rose',      struct('unit',true))
%   info = dc_comotion_view(ax, st, 'decompose')
%
% st : the rows of R.steps for ONE pair (from dc_comotion), in timepoint order.
%
% WHY NOT JUST THE QUIVER ON THE IMAGE. Arrows drawn where the molecules are make you compare two
% directions across the gap between them, by eye, at whatever angle the pair happens to sit. Two
% molecules 500 nm apart moving identically look like two unrelated arrow fields. Taking the offset
% out is the whole trick: once both are drawn from the same origin, "together" is a picture of
% overlap and needs no measurement.
%
% THE THREE VIEWS, and what each is for.
%
%   'centred'    Both trajectories translated so they start at the same point. Co-motion is then two
%                paths of the SAME SHAPE, which is the most direct thing to look at and the one to
%                show someone who has not seen the statistic. It keeps the magnitudes, so a pair
%                where one is dragging the other is visible as two similar shapes of different size.
%
%   'rose'       Every step vector drawn from one origin, both colours over each other. This is the
%                cosine made visible: the statistic is the mean of cos between paired arrows, and
%                the plot is those arrows. With 'unit' the lengths are thrown away and only the
%                directions remain, which is exactly what the cosine measures and stops one long
%                excursion dominating the picture the way it dominates a dot product. The thick
%                arrow on each colour is its mean direction; the angle between those two thick
%                arrows is the pair's co-motion at a glance.
%
%   'decompose'  Each pair of steps split into what the two SHARE and what they do not:
%
%                    common c = (u+v)/2        the motion both made
%                    relative d = (v-u)/2      the motion only one made
%
%                and u = c-d, v = c+d exactly, so nothing is lost. Drawn as two clouds. A pair that
%                moves together has a broad common cloud and a tight relative one; two independent
%                molecules have clouds of the same size. This is the view that separates "they move
%                together" from "they both move a lot".
%
%                It also gives a number that needs no null: shared = (<|c|^2> - <|d|^2>) /
%                (<|c|^2> + <|d|^2>), which is +1 for identical motion, 0 for independent and -1 for
%                exactly opposed. It equals the correlation of the step vectors, so it is the same
%                quantity the cosine estimates, computed without dividing by lengths.
%
% OUTPUT info: .shared .commonRms .relRms (um) .n .text, plus .meanAngleDeg between the two mean
%              directions and .meanResolved saying whether that angle means anything — for a pair
%              with no net drift the mean of its steps is noise and the angle between two such means
%              is uniformly random, so it is reported as NaN rather than as a measurement.
%
% .shared IS THE CORRELATION ITSELF, not a proxy for it. Measured against injected values on
% synthetic pairs it reads 0.016, 0.263, 0.608 and 0.889 for true 0, 0.3, 0.6 and 0.9. The mean
% cosine that dc_comotion reports estimates (pi/4)*rho instead, so the two differ by about 21% by
% construction; this one is the number to quote as a correlation.

if nargin < 4 || ~isstruct(opts), opts = struct(); end
cA = getf(opts,'colourA',[0.85 0.25 0.85]);
cB = getf(opts,'colourB',[0.20 0.70 0.30]);
unit = getf(opts,'unit',false);

info = struct('shared',NaN,'commonRms',NaN,'relRms',NaN,'meanAngleDeg',NaN, ...
              'meanResolved',false,'n',0,'text','');
if isempty(st), cla(ax); title(ax,'no shared steps'); return; end
st = sortrows(st,'tp');
u = [st.uax st.uay];  v = [st.ubx st.uby];
info.n = size(u,1);

c = (u + v)/2;  d = (v - u)/2;
pc = mean(sum(c.^2,2)); pd = mean(sum(d.^2,2));
info.shared    = (pc - pd) / max(pc + pd, eps);
info.commonRms = sqrt(pc);
info.relRms    = sqrt(pd);
% The MEAN direction of a set of steps only means something when the steps have a net drift. For a
% pair going nowhere on average the mean is a small vector made of noise, and the angle between two
% such vectors is uniformly random — a number that looks like a measurement and is not. Report it
% only when each mean clears its own standard error.
mu = mean(u,1); mv = mean(v,1);
seU = std(u,0,1)/sqrt(max(info.n,1));  seV = std(v,0,1)/sqrt(max(info.n,1));
info.meanResolved = hypot(mu(1),mu(2)) > 2*hypot(seU(1),seU(2)) && ...
                    hypot(mv(1),mv(2)) > 2*hypot(seV(1),seV(2));
if info.meanResolved
    info.meanAngleDeg = rad2deg(atan2(mu(1)*mv(2)-mu(2)*mv(1), mu*mv'));
else
    info.meanAngleDeg = NaN;
end

cla(ax); hold(ax,'on');
switch lower(char(mode))
    case 'centred'
        % Both paths from a common origin. Positions, not steps: the cumulative sum of the steps is
        % the trajectory each molecule made, with where it happened to be subtracted off.
        pa = [0 0; cumsum(u,1)];
        pb = [0 0; cumsum(v,1)];
        plot(ax, pa(:,1)*1000, pa(:,2)*1000, '-', 'Color',cA, 'LineWidth',1.3);
        plot(ax, pb(:,1)*1000, pb(:,2)*1000, '-', 'Color',cB, 'LineWidth',1.3);
        plot(ax, 0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize',4);
        plot(ax, pa(end,1)*1000, pa(end,2)*1000, 'o', 'Color',cA, 'MarkerFaceColor',cA, 'MarkerSize',5);
        plot(ax, pb(end,1)*1000, pb(end,2)*1000, 'o', 'Color',cB, 'MarkerFaceColor',cB, 'MarkerSize',5);
        xlabel(ax,'nm from start'); ylabel(ax,'nm from start');
        title(ax, sprintf('centred paths — same shape means together (shared %+.2f)', info.shared), ...
            'FontSize',9);

    case 'rose'
        a = u; b = v;
        if unit
            a = a ./ max(hypot(a(:,1),a(:,2)), eps);
            b = b ./ max(hypot(b(:,1),b(:,2)), eps);
            sc = 1; lbl = 'unit steps';
        else
            sc = 1000; lbl = 'nm';
        end
        quiver(ax, zeros(info.n,1), zeros(info.n,1), a(:,1)*sc, a(:,2)*sc, 0, ...
            'Color',[cA 0.35], 'LineWidth',0.7, 'MaxHeadSize',0.05);
        quiver(ax, zeros(info.n,1), zeros(info.n,1), b(:,1)*sc, b(:,2)*sc, 0, ...
            'Color',[cB 0.35], 'LineWidth',0.7, 'MaxHeadSize',0.05);
        % The mean direction of each, thick — but only when there IS one. Drawing a thick arrow
        % through a mean that is indistinguishable from zero invites reading a direction off noise.
        ma = mean(a,1)*sc; mb = mean(b,1)*sc;
        if info.meanResolved
            quiver(ax, 0,0, ma(1), ma(2), 0, 'Color',cA, 'LineWidth',2.4, 'MaxHeadSize',0.9);
            quiver(ax, 0,0, mb(1), mb(2), 0, 'Color',cB, 'LineWidth',2.4, 'MaxHeadSize',0.9);
            ttl = sprintf('steps from one origin — mean directions %.0f° apart', abs(info.meanAngleDeg));
        else
            ttl = 'steps from one origin — no net drift, so no mean direction to compare';
        end
        axis(ax,'equal');
        xlabel(ax, lbl); ylabel(ax, lbl);
        title(ax, ttl, 'FontSize',9);

    case 'decompose'
        % Two clouds: what they shared, and what they did not. Same axes, so their relative SIZE is
        % the answer — a tight relative cloud inside a broad common one is co-motion.
        plot(ax, c(:,1)*1000, c(:,2)*1000, '.', 'Color',[0.15 0.45 0.75], 'MarkerSize',9);
        plot(ax, d(:,1)*1000, d(:,2)*1000, '.', 'Color',[0.85 0.45 0.10], 'MarkerSize',9);
        th = linspace(0,2*pi,64);
        plot(ax, info.commonRms*1000*cos(th), info.commonRms*1000*sin(th), '-', ...
            'Color',[0.15 0.45 0.75], 'LineWidth',1.6);
        plot(ax, info.relRms*1000*cos(th), info.relRms*1000*sin(th), '-', ...
            'Color',[0.85 0.45 0.10], 'LineWidth',1.6);
        axis(ax,'equal');
        legend(ax, {'common (u+v)/2','relative (v-u)/2'}, 'Location','best','Box','off','FontSize',7);
        xlabel(ax,'nm'); ylabel(ax,'nm');
        title(ax, sprintf('shared %+.2f — common %.0f nm rms, relative %.0f nm', ...
            info.shared, 1000*info.commonRms, 1000*info.relRms), 'FontSize',9);

    otherwise
        error('dc_comotion_view:mode','unknown mode ''%s''', mode);
end
grid(ax,'on'); hold(ax,'off');

info.text = sprintf('%d steps; shared %+.3f (common %.0f nm rms against relative %.0f nm)%s', ...
    info.n, info.shared, 1000*info.commonRms, 1000*info.relRms, ...
    tern(info.meanResolved, sprintf('; mean directions %.0f° apart', abs(info.meanAngleDeg)), ...
         '; no net drift, so no mean direction'));
end

function y = tern(c,a,b), if c, y=a; else, y=b; end, end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
