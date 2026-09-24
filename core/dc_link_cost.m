function C = dc_link_cost(P, Q, R)
%DC_LINK_COST  Frame-to-frame link cost (np x nq): Euclidean distance, Inf beyond R.
%
%   C = dc_link_cost(P, Q, R)      P, Q: Nx2+ [x y ...] in px; R: max link distance in px
%
% Deliberately just distance. SPTinMatlab's linker can also bias a link by how much of it lies off a
% segmented ER, which is a statement about that biology; this toolkit tracks two particle channels
% and has no such mask. If a cost model is ever needed here it should be about the CHANNEL — a link
% that would cross into the other colour's territory — not about an organelle.
dx = P(:,1) - Q(:,1).';
dy = P(:,2) - Q(:,2).';
C  = hypot(dx, dy);
C(C > R) = Inf;
end
