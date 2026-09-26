function H = dc_app(opts)
%DC_APP  dcSPT — two colours, tracked independently, then compared.
%
%   dc_app
%   H = dc_app(struct('folder', '/path/to/cell', 'visible', 'off'))
%
% Four tabs, in the order the work happens, following the same shape as SPTinMatlab's tools: a top
% bar carrying the calibration that everything downstream reads, then
%
%   1 Cells        which stacks, what the acquisition says each page is, how the two colours line up
%   2 Detect       per-colour detection on a real frame, then tracking
%   3 Co-motion    the correlation-vs-separation curve, the null, and how many steps it takes
%   4 Pair         one cross-colour pair: each colour's bleaching trace and their co-motion
%
% WHY THE CALIBRATION IS IN THE TOP BAR AND NOT A TAB. Pixel size and frame interval are read from
% each stack's own metadata and are then used by every stage. Putting them where they are always
% visible means a wrong one is noticed before a run rather than after; the source of each is shown
% beside it, so a value that was typed is distinguishable from one the file supplied.
%
% WHAT THIS TOOL WILL NOT DO. It will not pair two colours that do not share timepoints, and it will
% not draw a co-motion panel for a same-colour pair. Both refusals are in the drivers rather than
% here, so a script gets them too.
%
% H (for tests and scripting): .fig and the handles the callbacks use, plus .api — a struct of
% function handles (.loadFolder .detectPreview .runTracking .runComotion .showPair) so the whole tool
% can be driven headlessly. Every button calls one of those and nothing else, which is what makes a
% test of the api a test of the app.

if nargin < 1 || ~isstruct(opts), opts = struct(); end
vis = getf(opts,'visible','on');

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root,'core'), fullfile(root,'drivers'));

% ---- state ---------------------------------------------------------------------------------------
St = struct('folder','', 'stacks',{{}}, 'C',[], 'D',[], 'Draw',[], 'R',[], 'N',[], 'S',[], ...
            'rejected',[], 'nDropped',0, ...
            'pxUm',0.0967821, 'dtS',NaN, 'pxSrc','default', 'dtSrc','default', 'align',[]);
% Thresholding follows SPTinMatlab: pool the DoG candidate qualities over sampled frames, then cut
% at a top percentile of THAT distribution. The histogram is the control; the number is read off it.
prm = struct('diamUm',0.4, 'thr',[10 10], 'thrMode','pct', 'topPct',[2 2], ...
             'linkUm',0.6, 'gapUm',0.9, 'maxGap',1, 'maxFrames',[], 'minLen',5);
pool = {[] []};            % the pooled candidate qualities, per colour
% rMaxUm is deliberately generous. It caps which pairs are kept at all, and the far-field null is
% the shell between rFarUm and rMaxUm — set it tight and there is no null to compare the near field
% against. Cost is quadratic in the tracks alive at once, so lower it if a run gets slow.
co  = struct('rMaxUm',5, 'nMin',20, 'rNearUm',0.5, 'rFarUm',2, 'drift',true);

% ---- window --------------------------------------------------------------------------------------
fig = uifigure('Name','dcSPT — two colours, tracked independently, then compared', ...
    'Position',[70 70 1180 760], 'Visible', vis);
gl  = uigridlayout(fig,[2 1],'RowHeight',{32,'1x'},'Padding',[8 8 8 8],'RowSpacing',6);
top = uigridlayout(gl,[1 9],'ColumnWidth',{'1x',94,76,116,84,58,70,64,60}, ...
    'Padding',[0 0 0 0],'ColumnSpacing',6);
uilabel(top,'Text','dcSPT — dual colour','FontWeight','bold','FontColor',[0.25 0.25 0.3]);
uilabel(top,'Text','Pixel (µm/px)','HorizontalAlignment','right');
eCalPx = uieditfield(top,'numeric','Value',St.pxUm,'ValueDisplayFormat','%.5g','Limits',[1e-4 10], ...
    'Tooltip','µm per pixel. Read from the stack''s own metadata when it carries it; edit to override.', ...
    'ValueChangedFcn',@(s,e) onCal('px'));
uilabel(top,'Text','Frame interval (s)','HorizontalAlignment','right');
eCalDt = uieditfield(top,'numeric','Value',0.012,'ValueDisplayFormat','%.5g','Limits',[1e-6 3600], ...
    'Tooltip',['Seconds per frame OF ONE COLOUR. Nothing in detection, tracking or co-motion indexes ' ...
               'by this — it is applied once, at the end, to turn a frame count into seconds.'], ...
    'ValueChangedFcn',@(s,e) onCal('dt'));
lblCal = uilabel(top,'Text','(no stack)','FontSize',11,'FontColor',[0.45 0.45 0.5]);
uibutton(top,'Text','Auto','Tooltip','Re-read pixel size and frame interval from the stacks'' own metadata.', ...
    'ButtonPushedFcn',@(s,e) onAuto());
uibutton(top,'Text','Folder…','Tooltip','Pick a folder holding this cell''s two colour stacks.', ...
    'ButtonPushedFcn',@(s,e) onPick());
uibutton(top,'Text','❓ Help','ButtonPushedFcn',@(s,e) onHelp());

tg = uitabgroup(gl); tg.Layout.Row = 2;
t1 = uitab(tg,'Title','1 · Cells');      buildCells(t1);
t2 = uitab(tg,'Title','2 · Detect');     buildDetect(t2);
t3 = uitab(tg,'Title','3 · Track');      buildTrack(t3);
t4 = uitab(tg,'Title','4 · Co-motion');  buildComotion(t4);
t5 = uitab(tg,'Title','5 · Pair');       buildPair(t5);

% ---- handles out ---------------------------------------------------------------------------------
H = struct('fig',fig, 'tabs',[t1 t2 t3 t4 t5], 'tg',tg);
% state/params/coParams are NESTED functions, not anonymous ones. An anonymous handle captures its
% variables BY VALUE where it is created, so @() St would hand back the empty struct this line sees
% and never the state the callbacks have since filled in.
H.api = struct('loadFolder',@loadFolder, 'detectPreview',@detectPreview, ...
               'runTracking',@runTracking, 'runComotion',@runComotion, ...
               'showPair',@showPair, 'state',@getState, 'params',@getParams, 'coParams',@getCo, ...
               'setParam',@setParam, 'setCo',@setCo, 'export',@onExport, ...
               'poolQuality',@poolQuality, 'showFrame',@showTrackFrame, 'curate',@curateTrack, ...
               'pool',@getPool);

if isfield(opts,'folder') && ~isempty(opts.folder), loadFolder(opts.folder); end

% ================================ tab 1: cells ====================================================
    function buildCells(parent)
        g = uigridlayout(parent,[3 1],'RowHeight',{92,'1x',120},'Padding',[10 10 10 10],'RowSpacing',8);
        hp = uigridlayout(g,[2 1],'RowHeight',{24,'1x'},'Padding',[0 0 0 0],'RowSpacing',4);
        uilabel(hp,'Text','The stacks, and what the acquisition says each page is','FontWeight','bold');
        lblFolder = uilabel(hp,'Text','No folder chosen. Use Folder… above, or drop a folder holding two colour stacks.', ...
            'FontColor',[0.45 0.45 0.5],'WordWrap','on');

        tblCh = uitable(g,'ColumnName',{'colour','stack','pages','frames','timepoints','dt (s)','source'}, ...
            'ColumnWidth',{62,'auto',70,70,110,84,90},'RowName',{});

        ap = uigridlayout(g,[2 1],'RowHeight',{22,'1x'},'Padding',[0 0 0 0],'RowSpacing',4);
        uilabel(ap,'Text','How the two colours line up','FontWeight','bold');
        lblAlign = uilabel(ap,'Text','—','WordWrap','on','FontColor',[0.2 0.2 0.25]);

        cellsCtl = struct('folder',lblFolder,'tbl',tblCh,'align',lblAlign);
        setappdata(fig,'cells',cellsCtl);
    end

% ================================ tab 2: detect & track ===========================================
    function buildDetect(parent)
        % Laid out like SPTinMatlab's Detect tab, doubled: the controls, then per colour a preview
        % over its pooled-quality histogram, then the log.
        g = uigridlayout(parent,[3 1],'RowHeight',{92,'1x',82},'Padding',[8 8 8 8],'RowSpacing',6);

        cp = uigridlayout(g,[3 1],'RowHeight',{26,26,24},'Padding',[0 0 0 0],'RowSpacing',5);
        r1 = uigridlayout(cp,[1 14],'ColumnWidth',{96,58,74,86,70,54,70,54,8,64,58,78,74,'1x'}, ...
            'Padding',[0 0 0 0],'ColumnSpacing',5);
        uilabel(r1,'Text','Spot diameter (µm)','HorizontalAlignment','right');
        spnDiam = uieditfield(r1,'numeric','Value',prm.diamUm,'Limits',[0.05 5], ...
            'ValueChangedFcn',@(s,e) onDiam(s.Value));
        uilabel(r1,'Text','Threshold','HorizontalAlignment','right');
        ddMode = uidropdown(r1,'Items',{'Top %','Quality ≥'},'ItemsData',{'pct','qual'}, ...
            'Value',prm.thrMode, ...
            'Tooltip',['Top %% cuts at a percentile of the POOLED candidate qualities — the ' ...
                       'histogram below. Quality ≥ sets the absolute DoG value directly, in the ' ...
                       'same units as that histogram''s x axis.'], ...
            'ValueChangedFcn',@(s,e) onThrMode(s.Value));
        uilabel(r1,'Text','c1','HorizontalAlignment','right');
        spnP1 = uieditfield(r1,'numeric','Value',prm.topPct(1),'Limits',[0.001 100], ...
            'ValueChangedFcn',@(s,e) onPct(1,s.Value));
        uilabel(r1,'Text','c2','HorizontalAlignment','right');
        spnP2 = uieditfield(r1,'numeric','Value',prm.topPct(2),'Limits',[0.001 100], ...
            'ValueChangedFcn',@(s,e) onPct(2,s.Value));
        uilabel(r1,'Text','');
        uilabel(r1,'Text','Frame','HorizontalAlignment','right');
        spnFrame = uieditfield(r1,'numeric','Value',1,'Limits',[1 1e6],'RoundFractionalValues',true, ...
            'ValueChangedFcn',@(s,e) detectPreview());
        uibutton(r1,'Text','Pool quality','Tooltip', ...
            ['Sample frames across the movie and pool their DoG candidate qualities, per colour. ' ...
             'That distribution is what a Top %% threshold is a percentile OF, so it has to be ' ...
             'built before the number means anything.'], ...
            'ButtonPushedFcn',@(s,e) poolQuality());
        uibutton(r1,'Text','Preview','ButtonPushedFcn',@(s,e) detectPreview());

        r2 = uigridlayout(cp,[1 6],'ColumnWidth',{140,'1x',110,150,110,150}, ...
            'Padding',[0 0 0 0],'ColumnSpacing',6);
        sldFrame = uislider(r2,'Limits',[1 100],'Value',1,'MajorTicks',[], ...
            'ValueChangingFcn',@(s,e) onSlide(e.Value));
        sldFrame.Layout.Column = [1 2];
        uilabel(r2,'Text','c1 threshold','HorizontalAlignment','right');
        lblT1 = uilabel(r2,'Text','—','FontName','Menlo','FontSize',11);
        uilabel(r2,'Text','c2 threshold','HorizontalAlignment','right');
        lblT2 = uilabel(r2,'Text','—','FontName','Menlo','FontSize',11);

        r3 = uigridlayout(cp,[1 1],'Padding',[0 0 0 0]);
        lblDet = uilabel(r3,'Text','Pick a folder, then Pool quality.','FontColor',[0.35 0.35 0.4]);

        mid = uigridlayout(g,[1 2],'ColumnWidth',{'1x','1x'},'Padding',[0 0 0 0],'ColumnSpacing',8);
        col = gobjects(2,2);
        for i = 1:2
            cg = uigridlayout(mid,[2 1],'RowHeight',{'1x',150},'Padding',[0 0 0 0],'RowSpacing',4);
            pn = uipanel(cg,'BorderType','none');          % see the note on the Track tab's axes
            col(i,1) = uiaxes(pn); col(i,1).Units = 'normalized';
            col(i,1).Position = [0.09 0.08 0.88 0.84];
            col(i,2) = uiaxes(cg);
        end

        lg = uigridlayout(g,[1 1],'Padding',[0 0 0 0]);
        txtLog = uitextarea(lg,'Editable','off','Value',{'ready.'},'FontName','Menlo','FontSize',11);

        setappdata(fig,'detect', struct('ax',col(:,1)','hist',col(:,2)','lbl',lblDet,'log',txtLog, ...
            'frame',spnFrame,'sld',sldFrame,'diam',spnDiam,'mode',ddMode, ...
            'pct',[spnP1 spnP2],'thrLbl',[lblT1 lblT2]));
    end

    function onDiam(v)
        prm.diamUm = v;
        pool = {[] []};                 % the candidate set depends on the diameter: the old pool is stale
        d = getappdata(fig,'detect');
        say('spot diameter %.3g µm — the pooled quality is stale, pool again before trusting Top %%', v);
        cla(d.hist(1)); cla(d.hist(2));
        detectPreview();
    end
    function onThrMode(m)
        prm.thrMode = m;
        d = getappdata(fig,'detect');
        isQ = strcmp(m,'qual');
        for i = 1:2
            d.pct(i).Limits = tern(isQ, [0 1e5], [0.001 100]);
            d.pct(i).Value  = tern(isQ, prm.thr(i), prm.topPct(i));
        end
        detectPreview();
    end
    function onPct(i, v)
        if strcmp(prm.thrMode,'qual'), prm.thr(i) = v; else, prm.topPct(i) = v; end
        detectPreview();
    end
    function onSlide(v)
        d = getappdata(fig,'detect');
        d.frame.Value = round(v);
        detectPreview();
    end

    function q = getPool(), q = pool; end

    function poolQuality()
        assert(~isempty(St.C), 'dc_app:noCells', 'pick a folder first');
        d = getappdata(fig,'detect');
        for i = 1:2
            say('pooling %s candidate qualities…', St.C(i).key);
            pool{i} = dc_pool_quality(St.stacks{i}, prm.diamUm, St.pxUm, 40, ...
                struct('pages', St.C(i).pages));
        end
        say('pooled %d and %d candidates', numel(pool{1}), numel(pool{2}));
        detectPreview();
    end

    function t = curThr(i)
        % The absolute DoG value in force for this colour, whichever way it was set.
        if strcmp(prm.thrMode,'qual'), t = prm.thr(i); return; end
        if isempty(pool{i}), t = prm.thr(i); return; end          % not pooled yet: last known
        t = prctile(pool{i}, 100 - prm.topPct(i));
        prm.thr(i) = t;
    end

    function out = detectPreview()
        if isempty(St.C), return; end
        d = getappdata(fig,'detect');
        k = max(1, round(d.frame.Value));
        out = cell(1,2);
        for i = 1:2
            pg = St.C(i).pages(min(k, numel(St.C(i).pages)));
            raw = double(imread(St.stacks{i}, pg));
            thr = curThr(i);
            xy = dc_detect(raw, prm.diamUm, St.pxUm, thr, struct());
            out{i} = xy;

            % Keep the view the user set. cla + imagesc resets the limits, so zooming in and then
            % stepping a frame would throw the zoom away — which is the whole reason to zoom.
            [xl, yl] = keepView(d.ax(i), size(raw));
            cla(d.ax(i));
            lo = prctile(raw(:),1); hi = prctile(raw(:),99.9);
            imagesc(d.ax(i), raw, [lo hi]); colormap(d.ax(i), gray);
            set(d.ax(i), 'DataAspectRatio',[1 1 1], 'YDir','reverse', ...
                         'XLim',[0.5 size(raw,2)+0.5], 'YLim',[0.5 size(raw,1)+0.5]);
            hold(d.ax(i),'on');
            dc_draw(d.ax(i), 'spots', xy, (prm.diamUm/St.pxUm)/2, tern(i==1,[1 0.3 1],[0.3 1 0.4]));
            hold(d.ax(i),'off');
            if ~isempty(xl), xlim(d.ax(i), xl); ylim(d.ax(i), yl); end
            title(d.ax(i), sprintf('%s  page %d  thr %.4g  %d spots', St.C(i).key, pg, thr, size(xy,1)));
            d.thrLbl(i).Text = sprintf('%.4g', thr);
            drawHist(d.hist(i), pool{i}, thr, St.C(i).key);
        end
        if isgraphics(d.sld) && numel(St.C(1).pages) > 1
            d.sld.Limits = [1 numel(St.C(1).pages)];
            d.sld.Value = min(max(k,1), numel(St.C(1).pages));
        end
        gate = tern(strcmp(prm.thrMode,'qual'), ...
            sprintf('quality ≥ %.4g / %.4g', prm.thr(1), prm.thr(2)), ...
            sprintf('top %.3g%% / %.3g%%', prm.topPct(1), prm.topPct(2)));
        d.lbl.Text = sprintf('frame %d · %s · %d spots in %s, %d in %s · pooled n = %d / %d', ...
            k, gate, size(out{1},1), St.C(1).key, size(out{2},1), St.C(2).key, ...
            numel(pool{1}), numel(pool{2}));
    end

    function drawHist(ax, q, thr, key)
        cla(ax);
        if isempty(q)
            title(ax, sprintf('%s — press Pool quality', key), 'FontSize', 9);
            return
        end
        histogram(ax, q, 60, 'FaceColor',[0.45 0.55 0.75], 'EdgeColor','none');
        set(ax,'YScale','log');
        xline(ax, thr, 'r-', 'LineWidth', 1.4);
        xlabel(ax,'DoG quality'); ylabel(ax,'count');
        title(ax, sprintf('%s pooled candidates n=%d — line is the threshold', key, numel(q)), 'FontSize', 9);
    end

    function [xl, yl] = keepView(ax, sz)
        % Return the current limits only if they are a genuine zoom, so a first draw still autoscales.
        xl = []; yl = [];
        if ~isgraphics(ax) || isempty(ax.Children), return; end
        x = ax.XLim; y = ax.YLim;
        if x(1) > 0.6 || y(1) > 0.6 || x(2) < sz(2)-0.4 || y(2) < sz(1)-0.4
            xl = x; yl = y;
        end
    end

% ================================ tab 3: track ====================================================
    function buildTrack(parent)
        % The tracks over the raw data, as a two-colour composite. Tracking is judged by watching it,
        % not by a count — a linker that swaps identities produces exactly the right number of tracks.
        g = uigridlayout(parent,[3 1],'RowHeight',{84,'1x',92},'Padding',[8 8 8 8],'RowSpacing',6);

        cp = uigridlayout(g,[3 1],'RowHeight',{26,26,22},'Padding',[0 0 0 0],'RowSpacing',5);
        r1 = uigridlayout(cp,[1 14],'ColumnWidth',{72,56,92,56,100,56,74,56,8,96,104,86,78,'1x'}, ...
            'Padding',[0 0 0 0],'ColumnSpacing',5);
        uilabel(r1,'Text','Link (µm)','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',prm.linkUm,'Limits',[0.01 20], ...
            'ValueChangedFcn',@(s,e) setParam('linkUm',s.Value));
        uilabel(r1,'Text','Gap close (µm)','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',prm.gapUm,'Limits',[0 20], ...
            'ValueChangedFcn',@(s,e) setParam('gapUm',s.Value));
        uilabel(r1,'Text','Max gap (frames)','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',prm.maxGap,'Limits',[0 20],'RoundFractionalValues',true, ...
            'ValueChangedFcn',@(s,e) setParam('maxGap',s.Value));
        uilabel(r1,'Text','Min length','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',prm.minLen,'Limits',[2 1e4],'RoundFractionalValues',true, ...
            'Tooltip','Tracks shorter than this are hidden and excluded from the analysis (curation, not deletion).', ...
            'ValueChangedFcn',@(s,e) onMinLen(s.Value));
        uilabel(r1,'Text','');
        uilabel(r1,'Text','Frames (0 = all)','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',0,'Limits',[0 1e6],'RoundFractionalValues',true, ...
            'ValueChangedFcn',@(s,e) setParam('maxFrames', tern(s.Value>0, s.Value, [])));
        uibutton(r1,'Text','Track both','FontWeight','bold','ButtonPushedFcn',@(s,e) runTracking());
        uibutton(r1,'Text','▶ Play','ButtonPushedFcn',@(s,e) onPlay());

        r2 = uigridlayout(cp,[1 8],'ColumnWidth',{'1x',120,104,74,96,74,96,86}, ...
            'Padding',[0 0 0 0],'ColumnSpacing',6);
        sldT = uislider(r2,'Limits',[1 100],'Value',1,'MajorTicks',[], ...
            'ValueChangingFcn',@(s,e) showTrackFrame(round(e.Value)));
        ddView = uidropdown(r2,'Items',{'composite','c1 only','c2 only'}, ...
            'ItemsData',{'composite','a','b'},'Value','composite', ...
            'ValueChangedFcn',@(s,e) showTrackFrame([]));
        uilabel(r2,'Text','Tail (frames)','HorizontalAlignment','right');
        spnTail = uieditfield(r2,'numeric','Value',40,'Limits',[1 1e5],'RoundFractionalValues',true, ...
            'Tooltip','How much of each track''s history to draw behind it.', ...
            'ValueChangedFcn',@(s,e) showTrackFrame([]));
        uilabel(r2,'Text','Gamma','HorizontalAlignment','right');
        spnGam = uieditfield(r2,'numeric','Value',0.7,'Limits',[0.1 3], ...
            'Tooltip','Display gamma on the composite. Pulls dim spots up without changing the limits.', ...
            'ValueChangedFcn',@(s,e) showTrackFrame([]));
        chkPaths = uicheckbox(r2,'Text','Paths','Value',true, ...
            'ValueChangedFcn',@(s,e) showTrackFrame([]));
        uibutton(r2,'Text','Reject track','Tooltip', ...
            ['Exclude the track nearest the last click from the analysis. Curation here BLANKS a ' ...
             'track rather than deleting it, so every other track keeps its number.'], ...
            'ButtonPushedFcn',@(s,e) onRejectClick());

        r3 = uigridlayout(cp,[1 1],'Padding',[0 0 0 0]);
        lblTrk = uilabel(r3,'Text','Track both colours to see them here.','FontColor',[0.35 0.35 0.4]);

        % The axes lives in a uipanel, not directly in the grid. A uiaxes that is a grid child
        % resizes ITSELF to the data aspect ratio and will happily grow past its cell, ending up
        % underneath whatever is in the next row; inside a panel it letterboxes within a fixed box.
        pnlT = uipanel(g,'BorderType','none');
        axT = uiaxes(pnlT); axT.Units = 'normalized'; axT.Position = [0.05 0.06 0.90 0.88];
        title(axT,'composite');
        axT.Toolbar.Visible = 'on';

        lg = uigridlayout(g,[1 1],'Padding',[0 0 0 0]);
        txtT = uitextarea(lg,'Editable','off','Value',{''},'FontName','Menlo','FontSize',11);

        setappdata(fig,'track', struct('ax',axT,'sld',sldT,'view',ddView,'tail',spnTail, ...
            'gamma',spnGam,'paths',chkPaths,'lbl',lblTrk,'log',txtT,'tp',1,'timer',[]));
    end

    function onMinLen(v)
        prm.minLen = v;
        applyCuration();
        showTrackFrame([]);
    end

    function applyCuration()
        % Curation BLANKS rather than deletes: a rejected track's spots keep their row and lose their
        % track id. Deleting would renumber every track after it, and every pair id already quoted in
        % the Co-motion tab would then point somewhere else.
        if isempty(St.D), return; end
        if ~isfield(St,'rejected') || isempty(St.rejected), St.rejected = []; end
        S = St.Draw.spots;                       % always re-derive from the untouched build
        keep = true(height(S),1);
        n = groupcounts(S(isfinite(S.trackId),:), 'trackId');
        short = n.trackId(n.GroupCount < prm.minLen);
        drop = unique([short(:); St.rejected(:)]);
        blank = ismember(S.trackId, drop);
        S.trackId(blank) = NaN; S.trackLocal(blank) = NaN;
        St.D.spots = S;
        St.nDropped = numel(drop);
    end

    function curateTrack(id)
        if ~isfield(St,'rejected'), St.rejected = []; end
        if ismember(id, St.rejected), St.rejected(St.rejected == id) = [];
        else, St.rejected(end+1) = id; end
        applyCuration();
        sayT('%d track(s) excluded (min length %d, %d rejected by hand)', ...
            St.nDropped, prm.minLen, numel(St.rejected));
        showTrackFrame([]);
    end

    function onRejectClick()
        t = getappdata(fig,'track');
        cp = t.ax.CurrentPoint;
        if isempty(cp) || isempty(St.D), return; end
        xy = cp(1,1:2) * St.pxUm;
        S = St.D.spots(St.D.spots.tp == t.tp & isfinite(St.D.spots.trackId), :);
        if isempty(S), sayT('no track in this frame to reject'); return; end
        [dmin, j] = min(hypot(S.x - xy(1), S.y - xy(2)));
        if dmin > 1.0, sayT('nearest track is %.2f µm away — click closer', dmin); return; end
        curateTrack(S.trackId(j));
    end

    function showTrackFrame(tp)
        t = getappdata(fig,'track');
        if isempty(St.D) || ~isgraphics(t.ax), return; end
        if isempty(tp), tp = t.tp; end
        tp = max(1, round(tp));
        nTp = max(St.C(1).tp);
        tp = min(tp, nTp);
        t.tp = tp; setappdata(fig,'track',t);

        pgA = pageFor(1, tp); pgB = pageFor(2, tp);
        if isempty(pgA) || isempty(pgB), return; end
        A = double(imread(St.stacks{1}, pgA));
        B = double(imread(St.stacks{2}, pgB));
        rgb = dc_composite(A, B, struct('mode', t.view.Value, 'gamma', t.gamma.Value));

        [xl, yl] = keepView(t.ax, size(A));
        cla(t.ax);
        image(t.ax, rgb);
        % 'axis image' resizes the AXES to the data's aspect, and inside a uigridlayout that lets it
        % grow past its cell and sit under the log box. Fix the data aspect and the limits instead,
        % and let the layout keep owning the box.
        set(t.ax, 'DataAspectRatio',[1 1 1], 'YDir','reverse', ...
                  'XLim',[0.5 size(rgb,2)+0.5], 'YLim',[0.5 size(rgb,1)+0.5]);
        hold(t.ax,'on');
        if t.paths.Value, tail = t.tail.Value; else, tail = 1; end
        hA = dc_draw(t.ax, 'tracks', St.D, tp, struct('ch',chKey(1),'colour',[1 0.45 1], ...
            'pxUm',St.pxUm, 'rPx',(prm.diamUm/St.pxUm)/2, 'tail',tail));
        hB = dc_draw(t.ax, 'tracks', St.D, tp, struct('ch',chKey(2),'colour',[0.45 1 0.5], ...
            'pxUm',St.pxUm, 'rPx',(prm.diamUm/St.pxUm)/2, 'tail',tail));
        hold(t.ax,'off');
        if ~isempty(xl), xlim(t.ax, xl); ylim(t.ax, yl); end
        if isgraphics(t.sld), t.sld.Limits = [1 max(nTp,2)]; t.sld.Value = tp; end
        title(t.ax, sprintf('timepoint %d of %d — magenta %s, green %s', tp, nTp, chKey(1), chKey(2)));
        nDrop = 0; if isfield(St,'nDropped'), nDrop = St.nDropped; end
        t.lbl.Text = sprintf(['timepoint %d · %d %s tracks and %d %s tracks visible · %d excluded ' ...
            '(shorter than %d, or rejected)'], tp, hA.n, chKey(1), hB.n, chKey(2), nDrop, prm.minLen);
    end

    function k = chKey(i)
        k = ''; if ~isempty(St.C) && numel(St.C) >= i, k = St.C(i).key; end
    end
    function pg = pageFor(i, tp)
        pg = []; if isempty(St.C), return; end
        j = find(St.C(i).tp == tp, 1);
        if ~isempty(j), pg = St.C(i).pages(j); end
    end

    function onPlay()
        t = getappdata(fig,'track');
        if ~isempty(t.timer) && isvalid(t.timer)
            stop(t.timer); delete(t.timer); t.timer = []; setappdata(fig,'track',t);
            sayT('stopped'); return
        end
        t.timer = timer('ExecutionMode','fixedSpacing','Period',0.08,'TimerFcn',@(~,~) tick());
        setappdata(fig,'track',t); start(t.timer);
    end
    function tick()
        t = getappdata(fig,'track');
        if ~isgraphics(fig) || isempty(St.D)
            if ~isempty(t.timer) && isvalid(t.timer), stop(t.timer); delete(t.timer); end
            return
        end
        nTp = max(St.C(1).tp);
        showTrackFrame(mod(t.tp, nTp) + 1);
    end
    function sayT(varargin)
        t = getappdata(fig,'track');
        if isempty(t) || ~isgraphics(t.log), return; end
        t.log.Value = [t.log.Value; {sprintf(varargin{:})}];
        drawnow limitrate;
    end

% ================================ tab 4: co-motion ================================================
    function buildComotion(parent)
        g = uigridlayout(parent,[3 1],'RowHeight',{62,'1x',96},'Padding',[10 10 10 10],'RowSpacing',8);

        r = uigridlayout(g,[2 1],'RowHeight',{26,26},'Padding',[0 0 0 0],'RowSpacing',5);
        r1 = uigridlayout(r,[1 12],'ColumnWidth',{84,60,94,60,84,60,80,60,104,'1x',126,86}, ...
            'Padding',[0 0 0 0],'ColumnSpacing',6);
        uilabel(r1,'Text','Max r (µm)','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',co.rMaxUm,'Limits',[0.1 20], ...
            'ValueChangedFcn',@(s,e) setCo('rMaxUm',s.Value));
        uilabel(r1,'Text','Near ≤ (µm)','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',co.rNearUm,'Limits',[0.02 20], ...
            'ValueChangedFcn',@(s,e) setCo('rNearUm',s.Value));
        uilabel(r1,'Text','Far > (µm)','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',co.rFarUm,'Limits',[0.05 20], ...
            'Tooltip','Pairs beyond this build the null: same movie, same drift, too far to interact.', ...
            'ValueChangedFcn',@(s,e) setCo('rFarUm',s.Value));
        uilabel(r1,'Text','Min steps','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',co.nMin,'Limits',[3 1e5],'RoundFractionalValues',true, ...
            'ValueChangedFcn',@(s,e) setCo('nMin',s.Value));
        chkDrift = uicheckbox(r1,'Text','Remove drift','Value',co.drift, ...
            'Tooltip',['Subtract the per-timepoint median step. Drift correlates EVERY pair at every ' ...
                       'separation and biases the mean, so more data makes it look more significant.'], ...
            'ValueChangedFcn',@(s,e) setCo('drift',s.Value));
        uilabel(r1,'Text','');
        uibutton(r1,'Text','Run co-motion','FontWeight','bold','ButtonPushedFcn',@(s,e) runComotion());
        uibutton(r1,'Text','Export…','ButtonPushedFcn',@(s,e) onExport());

        r2 = uigridlayout(r,[1 1],'Padding',[0 0 0 0]);
        lblCo = uilabel(r2,'Text','Cross-colour pairs only.','FontColor',[0.35 0.35 0.4],'WordWrap','on');

        mid = uigridlayout(g,[1 2],'ColumnWidth',{'1.15x','1x'},'Padding',[0 0 0 0],'ColumnSpacing',8);
        axC = uiaxes(mid); title(axC,'co-motion vs separation');
        rp  = uigridlayout(mid,[2 1],'RowHeight',{'1x',150},'Padding',[0 0 0 0],'RowSpacing',6);
        tblP = uitable(rp,'ColumnName',{'A','B','n','r (nm)','cos','z'}, ...
            'ColumnWidth',{52,52,56,66,66,56},'RowName',{}, ...
            'CellSelectionCallback',@(s,e) onPickPair(e));
        axN = uiaxes(rp); title(axN,'steps needed');

        lg = uigridlayout(g,[1 1],'Padding',[0 0 0 0]);
        txtCo = uitextarea(lg,'Editable','off','Value',{''},'FontName','Menlo','FontSize',11);

        setappdata(fig,'co', struct('axC',axC,'axN',axN,'tbl',tblP,'lbl',lblCo,'log',txtCo,'drift',chkDrift));
    end

% ================================ tab 4: pair =====================================================
    function buildPair(parent)
        g = uigridlayout(parent,[2 1],'RowHeight',{30,'1x'},'Padding',[8 8 8 8],'RowSpacing',6);
        hr = uigridlayout(g,[1 4],'ColumnWidth',{300,90,90,'1x'},'Padding',[0 0 0 0],'ColumnSpacing',6);
        lblPair = uilabel(hr,'Text','Select a pair in the Co-motion tab.','FontColor',[0.35 0.35 0.4]);
        uibutton(hr,'Text','◀ prev','ButtonPushedFcn',@(s,e) stepPair(-1));
        uibutton(hr,'Text','next ▶','ButtonPushedFcn',@(s,e) stepPair(+1));
        host = uipanel(g,'BorderType','none');
        setappdata(fig,'pair', struct('host',host,'lbl',lblPair,'k',0));
    end

% ================================ actions =========================================================
    function loadFolder(folder)
        St.folder = char(folder);
        d = [dir(fullfile(St.folder,'*.tif')); dir(fullfile(St.folder,'*.tiff'))];
        assert(numel(d) >= 2, 'dc_app:needTwo', ...
            '%s holds %d TIFF(s); two colour stacks are needed.', St.folder, numel(d));
        [~, ord] = sort({d.name}); d = d(ord);
        St.stacks = arrayfun(@(q) fullfile(q.folder,q.name), d(1:2), 'uni', 0);

        specs = struct('key',{},'label',{},'pages',{},'tp',{},'dt_s',{});
        rows = {};
        dtAll = [];
        for i = 1:2
            p = St.stacks{i};
            L = dc_tiff_labels(p);
            info = imfinfo(p); nPg = numel(info);
            key = sprintf('c%d', i); src = 'page order';
            if L.ok && ~isempty(L.channels)
                key = sprintf('c%d', L.channels(1)); src = 'slice labels';
                tp = L.tp(:); pages = (1:nPg)';
                if numel(L.channels) > 1
                    % one file holding both colours: take this file's own first channel
                    m = L.ch == L.channels(1); pages = find(m); tp = L.tp(m);
                end
            else
                pages = (1:nPg)'; tp = (1:nPg)';
            end
            dtp = NaN;
            try, cb = dc_tiff_calib(p);
                 if isfield(cb,'pixUm') && isfinite(cb.pixUm), St.pxUm = cb.pixUm; St.pxSrc = 'metadata'; end
                 if isfield(cb,'dt_s') && isfinite(cb.dt_s), dtp = cb.dt_s; end
            catch, end
            dtAll(end+1) = dtp; %#ok<AGROW>
            specs(i) = struct('key',key,'label',key,'pages',pages,'tp',tp,'dt_s',dtp);
            [~,nm,ex] = fileparts(p);
            % Counts as integers and dt to its real precision: a uitable renders a double as
            % "60.0000", which reads as a measurement rather than a count.
            rows(end+1,:) = {key, [nm ex], sprintf('%d', nPg), sprintf('%d', numel(pages)), ...
                sprintf('%d–%d', min(tp), max(tp)), sprintf('%.6g', dtp), src}; %#ok<AGROW>
        end
        St.C = dc_channels('manual', specs);
        if any(isfinite(dtAll)), St.dtS = median(dtAll(isfinite(dtAll))); St.dtSrc = 'metadata'; end

        c = getappdata(fig,'cells');
        c.folder.Text = St.folder;
        c.tbl.Data = rows;
        A = dc_align(St.C(1), St.C(2));
        St.align = A;
        c.align.Text = A.text;
        eCalPx.Value = St.pxUm;
        if isfinite(St.dtS), eCalDt.Value = St.dtS; end
        lblCal.Text = sprintf('%s / %s', St.pxSrc, St.dtSrc);
        say('loaded %s: %s and %s', St.folder, specs(1).key, specs(2).key);
        say('%s', A.text);
    end

    function runTracking()
        assert(~isempty(St.C), 'dc_app:noCells', 'pick a folder first');
        d = getappdata(fig,'detect');
        say('tracking both colours…');
        cel = struct('stack', St.stacks{1}, 'base', 'cell', 'stacks', struct());
        for i = 1:2, cel.stacks.(St.C(i).key) = St.stacks{i}; end
        P = struct('diamUm',prm.diamUm, 'thrAbs',curThr(1), 'pxUm',St.pxUm, ...
                   'linkUm',prm.linkUm, 'gapUm',prm.gapUm, 'maxGap',prm.maxGap);
        if ~isempty(prm.maxFrames), P.maxFrames = prm.maxFrames; end
        % A cell, not a struct array: assigning a filled struct into repmat(struct(),1,2) fails with
        % "dissimilar structures" the moment the two have different fields.
        R = cell(1,2);
        for i = 1:2
            Pi = P; Pi.thrAbs = curThr(i);        % each colour's own threshold, however it was set
            R{i} = dc_process_cell(cel, St.C(i), Pi);
            say('  %s: %d detections, %d tracks over %d frames', R{i}.key, R{i}.nDets, R{i}.nTracks, R{i}.nFrames);
        end
        D = dc_dataset('new', 'cell', St.C);
        for i = 1:2, D = dc_dataset('addChannel', D, R{i}); end
        St.Draw = D;          % the build as tracked, never edited
        St.D    = D;          % the curated view every analysis reads
        St.rejected = [];
        applyCuration();
        say('dataset: %d spots, %d tracks (%d excluded: shorter than %d)', ...
            height(D.spots), numel(unique(D.spots.trackId(isfinite(D.spots.trackId)))), ...
            St.nDropped, prm.minLen);
        showTrackFrame(1);
        tg.SelectedTab = t3;
    end

    function runComotion()
        assert(~isempty(St.D), 'dc_app:noTracks', 'track first (tab 2)');
        c = getappdata(fig,'co');
        [S, si] = dc_steps(St.D);
        sayCo('%s', si.text);
        if si.nAll > 0 && si.nKept < 0.5*si.nAll
            sayCo(['  NOTE: over half the steps were gap-closed across more than one timepoint, so ' ...
                   'they cannot be paired. Lower Max gap on the Track tab and re-track if the pair ' ...
                   'count is thin.']);
        end
        if co.drift, [S, Dr] = dc_drift(S); sayCo('%s', Dr.text); end
        St.S = S;
        R = dc_comotion(S, struct('rMaxUm',co.rMaxUm, 'nMin',co.nMin, 'classes',"cross"));
        St.R = R;
        % The null can legitimately refuse: too few distant pairs to build one from. That is a
        % statement about the movie, not a crash, so it belongs in the log next to the curve it
        % could not produce.
        N = [];
        try
            N = dc_comotion_null(R, struct('rNearUm',co.rNearUm, 'rFarUm',co.rFarUm));
        catch ME
            sayCo('no null: %s', ME.message);
        end
        St.N = N;

        b = R.bins(R.bins.class=="cross", :);
        cla(c.axC); hold(c.axC,'on');
        errorbar(c.axC, b.rMid*1000, b.meanCos, 2*b.seCos, 'o-', 'Color',[0.80 0.20 0.20], ...
            'MarkerFaceColor',[0.80 0.20 0.20], 'MarkerSize',4, 'LineWidth',1.2);
        yline(c.axC, 0, '-', 'Color',[0.75 0.75 0.75]);
        ttl = 'cross-colour co-motion';
        if ~isempty(N)
            yline(c.axC, N.far.meanCos, '--', 'Color',[0.35 0.4 0.45]);
            ttl = sprintf('%s — dashed line is the far-field null (%+.4f)', ttl, N.far.meanCos);
        else
            ttl = sprintf('%s — no far-field null (too few distant pairs)', ttl);
        end
        hold(c.axC,'off'); grid(c.axC,'on');
        xlabel(c.axC,'separation (nm)'); ylabel(c.axC,'mean cos\theta');
        title(c.axC, ttl);

        cla(c.axN);
        if ~isempty(N)
            hold(c.axN,'on');
            plot(c.axN, N.curve.n, N.curve.seFar, 'o-','Color',[0.15 0.45 0.70],'MarkerSize',3);
            plot(c.axN, N.curve.n, N.k.ideal./sqrt(N.curve.n), '--','Color',[0.6 0.6 0.6]);
            if isfinite(N.excess) && N.excess > 0, yline(c.axN, N.excess, ':', 'Color',[0.8 0.3 0.1]); end
            set(c.axN,'XScale','log','YScale','log'); grid(c.axN,'on'); hold(c.axN,'off');
            xlabel(c.axN,'steps in a pair'); ylabel(c.axN,'null SE');
            title(c.axN, sprintf('SE = %.3f/\\surdn (ideal %.3f)', N.k.fitted, N.k.ideal));
        else
            title(c.axN, 'steps needed — needs a far-field null');
        end

        P = R.pairs;
        c.tbl.Data = [num2cell(P.trackA), num2cell(P.trackB), num2cell(P.n), ...
                      num2cell(round(P.rMedian*1000)), num2cell(round(P.meanCos,3)), ...
                      num2cell(round(P.z,2))];
        c.lbl.Text = sprintf('%d cross-colour pairs with %d+ shared steps; %d step pairs in total', ...
            height(P), co.nMin, height(R.steps));
        if ~isempty(N), sayCo('%s', N.text); end
        if ~isempty(P), showPair(1); end
    end

    function showPair(k)
        if isempty(St.R) || isempty(St.R.pairs) || k < 1 || k > height(St.R.pairs), return; end
        pp = getappdata(fig,'pair'); pp.k = k; setappdata(fig,'pair',pp);
        dtv = []; if isfinite(St.dtS), dtv = St.dtS; end
        Hp = dc_pair_panel(St.D, St.R, k, struct('parent',pp.host,'null',St.N,'dtS',dtv));
        pp.lbl.Text = sprintf('pair %d of %d — tracks %g and %g, %d shared steps, cos %+.3f', ...
            k, height(St.R.pairs), Hp.trackA, Hp.trackB, Hp.n, Hp.meanCos);
        tg.SelectedTab = t4;
    end

    function stepPair(d)
        pp = getappdata(fig,'pair');
        if isempty(St.R) || isempty(St.R.pairs), return; end
        showPair(min(max(pp.k + d, 1), height(St.R.pairs)));
    end

    function onPickPair(e)
        if isempty(e.Indices), return; end
        showPair(e.Indices(1));
    end

    function onExport()
        if isempty(St.R), say('nothing to export yet'); return; end
        f = uigetdir(tern(isempty(St.folder), pwd, St.folder), 'Where to write the co-motion tables');
        if isequal(f,0), return; end
        writetable(St.R.pairs, fullfile(f,'dc_comotion_pairs.csv'));
        writetable(St.R.bins,  fullfile(f,'dc_comotion_bins.csv'));
        if ~isempty(St.N), writetable(St.N.curve, fullfile(f,'dc_comotion_null.csv')); end
        say('wrote dc_comotion_pairs.csv, _bins.csv and _null.csv to %s', f);
    end

% ================================ small helpers ===================================================
    function v = getState(), v = St;  end
    function v = getParams(), v = prm; end
    function v = getCo(),     v = co;  end
    function setParam(f,v), prm.(f) = v; end
    function setThr(i,v), prm.thr(i) = v; end
    function setCo(f,v), co.(f) = v; end
    function onCal(which)
        if strcmp(which,'px'), St.pxUm = eCalPx.Value; St.pxSrc = 'typed';
        else, St.dtS = eCalDt.Value; St.dtSrc = 'typed'; end
        lblCal.Text = sprintf('%s / %s', St.pxSrc, St.dtSrc);
    end
    function onAuto()
        if isempty(St.stacks), return; end
        loadFolder(St.folder);
    end
    function onPick()
        f = uigetdir(pwd, 'Folder holding this cell''s two colour stacks');
        if isequal(f,0), return; end
        loadFolder(f);
    end
    function onHelp()
        p = fullfile(root,'README.md');
        if isfile(p), web(p,'-browser'); end
    end
    function say(varargin)
        d = getappdata(fig,'detect');
        if isempty(d) || ~isgraphics(d.log), return; end
        d.log.Value = [d.log.Value; {sprintf(varargin{:})}];
        scroll(d.log,'bottom'); drawnow limitrate;
    end
    function sayCo(varargin)
        c = getappdata(fig,'co');
        if isempty(c) || ~isgraphics(c.log), return; end
        c.log.Value = [c.log.Value; {sprintf(varargin{:})}];
        drawnow limitrate;
    end
end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
function y = tern(c,a,b), if c, y=a; else, y=b; end, end
function i = argmin(v), [~,i] = min(v); end
function i = argmax(v), [~,i] = max(v); end
