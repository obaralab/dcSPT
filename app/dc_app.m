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
St = struct('folder','', 'stacks',{{}}, 'C',[], 'D',[], 'R',[], 'N',[], 'S',[], ...
            'pxUm',0.0967821, 'dtS',NaN, 'pxSrc','default', 'dtSrc','default', 'align',[]);
prm = struct('diamUm',0.4, 'thr',[10 10], 'noiseK',5, 'linkUm',0.6, 'gapUm',0.9, ...
             'maxGap',1, 'maxFrames',[]);
co  = struct('rMaxUm',2, 'nMin',20, 'rNearUm',0.5, 'rFarUm',1.5, 'drift',true);

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
t3 = uitab(tg,'Title','3 · Co-motion');  buildComotion(t3);
t4 = uitab(tg,'Title','4 · Pair');       buildPair(t4);

% ---- handles out ---------------------------------------------------------------------------------
H = struct('fig',fig, 'tabs',[t1 t2 t3 t4], 'tg',tg);
% state/params/coParams are NESTED functions, not anonymous ones. An anonymous handle captures its
% variables BY VALUE where it is created, so @() St would hand back the empty struct this line sees
% and never the state the callbacks have since filled in.
H.api = struct('loadFolder',@loadFolder, 'detectPreview',@detectPreview, ...
               'runTracking',@runTracking, 'runComotion',@runComotion, ...
               'showPair',@showPair, 'state',@getState, 'params',@getParams, 'coParams',@getCo, ...
               'setParam',@setParam, 'setCo',@setCo, 'calibrate',@onCalibrate, 'export',@onExport);

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
        g = uigridlayout(parent,[3 1],'RowHeight',{104,'1x',96},'Padding',[10 10 10 10],'RowSpacing',8);

        cp = uigridlayout(g,[3 1],'RowHeight',{28,28,28},'Padding',[0 0 0 0],'RowSpacing',6);
        r1 = uigridlayout(cp,[1 14],'ColumnWidth',{88,64,88,58,88,58,62,50,58,58,8,66,72,'1x'}, ...
            'Padding',[0 0 0 0],'ColumnSpacing',5);
        uilabel(r1,'Text','Spot diameter (µm)','HorizontalAlignment','right');
        spnDiam = uieditfield(r1,'numeric','Value',prm.diamUm,'Limits',[0.05 5], ...
            'ValueChangedFcn',@(s,e) setParam('diamUm',s.Value));
        uilabel(r1,'Text','Threshold c1','HorizontalAlignment','right');
        spnT1 = uieditfield(r1,'numeric','Value',prm.thr(1),'Limits',[0 1e5], ...
            'Tooltip','DoG response a local maximum must reach. Calibrate sets it from the response distribution.', ...
            'ValueChangedFcn',@(s,e) setThr(1,s.Value));
        uilabel(r1,'Text','Threshold c2','HorizontalAlignment','right');
        spnT2 = uieditfield(r1,'numeric','Value',prm.thr(2),'Limits',[0 1e5], ...
            'ValueChangedFcn',@(s,e) setThr(2,s.Value));
        uilabel(r1,'Text','Noise k','HorizontalAlignment','right');
        spnK = uieditfield(r1,'numeric','Value',prm.noiseK,'Limits',[1 30], ...
            'Tooltip', ['How many robust noise sigmas above the DoG background a spot must reach. ' ...
                        'Calibrate sets each colour''s threshold from its OWN noise, so a dim ' ...
                        'channel and a bright one are not held to the same absolute number.'], ...
            'ValueChangedFcn',@(s,e) setParam('noiseK',s.Value));
        uilabel(r1,'Text','Frame','HorizontalAlignment','right');
        spnFrame = uieditfield(r1,'numeric','Value',1,'Limits',[1 1e6],'RoundFractionalValues',true, ...
            'ValueChangedFcn',@(s,e) detectPreview());
        uilabel(r1,'Text','');
        uibutton(r1,'Text','Calibrate','Tooltip', ...
            ['Set each colour''s threshold from its own DoG response distribution, so a dim channel ' ...
             'and a bright one are not held to one number.'], ...
            'ButtonPushedFcn',@(s,e) onCalibrate());
        uibutton(r1,'Text','Preview','ButtonPushedFcn',@(s,e) detectPreview());

        r2 = uigridlayout(cp,[1 10],'ColumnWidth',{88,64,96,64,96,64,96,64,'1x',110}, ...
            'Padding',[0 0 0 0],'ColumnSpacing',6);
        uilabel(r2,'Text','Link (µm)','HorizontalAlignment','right');
        uieditfield(r2,'numeric','Value',prm.linkUm,'Limits',[0.01 20], ...
            'ValueChangedFcn',@(s,e) setParam('linkUm',s.Value));
        uilabel(r2,'Text','Gap close (µm)','HorizontalAlignment','right');
        uieditfield(r2,'numeric','Value',prm.gapUm,'Limits',[0 20], ...
            'ValueChangedFcn',@(s,e) setParam('gapUm',s.Value));
        uilabel(r2,'Text','Max gap (frames)','HorizontalAlignment','right');
        uieditfield(r2,'numeric','Value',prm.maxGap,'Limits',[0 20],'RoundFractionalValues',true, ...
            'ValueChangedFcn',@(s,e) setParam('maxGap',s.Value));
        uilabel(r2,'Text','Frames (blank = all)','HorizontalAlignment','right');
        uieditfield(r2,'numeric','Value',0,'Limits',[0 1e6],'RoundFractionalValues',true, ...
            'Tooltip','Track only the first N frames. 0 means all — use a few hundred to judge settings first.', ...
            'ValueChangedFcn',@(s,e) setParam('maxFrames', tern(s.Value>0, s.Value, [])));
        uilabel(r2,'Text','');
        btnRun = uibutton(r2,'Text','Track both colours','FontWeight','bold', ...
            'ButtonPushedFcn',@(s,e) runTracking());

        r3 = uigridlayout(cp,[1 1],'Padding',[0 0 0 0]);
        lblDet = uilabel(r3,'Text','Pick a folder, then Calibrate and Preview.','FontColor',[0.35 0.35 0.4]);

        axp = uigridlayout(g,[1 2],'ColumnWidth',{'1x','1x'},'Padding',[0 0 0 0],'ColumnSpacing',8);
        ax1 = uiaxes(axp); title(ax1,'c1'); ax2 = uiaxes(axp); title(ax2,'c2');

        lg = uigridlayout(g,[1 1],'Padding',[0 0 0 0]);
        txtLog = uitextarea(lg,'Editable','off','Value',{'ready.'},'FontName','Menlo','FontSize',11);

        setappdata(fig,'detect', struct('ax',[ax1 ax2],'lbl',lblDet,'log',txtLog,'frame',spnFrame, ...
            'thr',[spnT1 spnT2],'diam',spnDiam,'run',btnRun));
    end

% ================================ tab 3: co-motion ================================================
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
            rows(end+1,:) = {key, [nm ex], nPg, numel(pages), ...
                sprintf('%d–%d', min(tp), max(tp)), dtp, src}; %#ok<AGROW>
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

    function out = detectPreview()
        assert(~isempty(St.C), 'dc_app:noCells', 'pick a folder first');
        d = getappdata(fig,'detect');
        k = max(1, round(d.frame.Value));
        out = cell(1,2);
        for i = 1:2
            pg = St.C(i).pages(min(k, numel(St.C(i).pages)));
            raw = double(imread(St.stacks{i}, pg));
            xy = dc_detect(raw, prm.diamUm, St.pxUm, prm.thr(i), struct());
            out{i} = xy;
            cla(d.ax(i));
            lo = prctile(raw(:),1); hi = prctile(raw(:),99.9);
            imagesc(d.ax(i), raw, [lo hi]); colormap(d.ax(i), gray); axis(d.ax(i),'image');
            hold(d.ax(i),'on');
            if ~isempty(xy), plot(d.ax(i), xy(:,1), xy(:,2), 'o', 'Color',[1 0.35 0.1], 'MarkerSize',7); end
            hold(d.ax(i),'off');
            title(d.ax(i), sprintf('%s  page %d  thr %g  %d spots', St.C(i).key, pg, prm.thr(i), size(xy,1)));
        end
        d.lbl.Text = sprintf('frame %d: %d spots in %s, %d in %s', k, ...
            size(out{1},1), St.C(1).key, size(out{2},1), St.C(2).key);
    end

    function onCalibrate()
        assert(~isempty(St.C), 'dc_app:noCells', 'pick a folder first');
        % Each colour gets its own threshold from its OWN response distribution. A dim channel and a
        % bright one held to one number is how a weak colour ends up with no tracks and a strong one
        % with a field of noise.
        % Each colour's threshold from its OWN noise, not from a percentile of its maxima. A
        % percentile assumes how many spots there are: on a field with more particles than the
        % percentile allows for, it silently keeps the brightest few and drops the rest, and the
        % tracks then fragment because each particle is detected only intermittently. The robust
        % sigma of the DoG image says nothing about how many spots there are, only about the floor
        % they have to clear.
        yld = zeros(1,2);
        for i = 1:2
            ts = round(linspace(1, numel(St.C(i).pages), min(8, numel(St.C(i).pages))));
            mu = zeros(1,numel(ts)); sg = zeros(1,numel(ts));
            for j = 1:numel(ts)
                raw = double(imread(St.stacks{i}, St.C(i).pages(ts(j))));
                gg = dc_dog(raw, prm.diamUm, St.pxUm);
                mu(j) = median(gg(:));
                sg(j) = 1.4826 * median(abs(gg(:) - mu(j)));
            end
            prm.thr(i) = round(median(mu) + prm.noiseK*median(sg), 2);
            n = 0;
            for j = 1:numel(ts)
                raw = double(imread(St.stacks{i}, St.C(i).pages(ts(j))));
                n = n + size(dc_detect(raw, prm.diamUm, St.pxUm, prm.thr(i), struct()), 1);
            end
            yld(i) = n / numel(ts);
        end
        d = getappdata(fig,'detect');
        d.thr(1).Value = prm.thr(1); d.thr(2).Value = prm.thr(2);
        say('calibrated at k = %g sigma: %s -> %.2f (%.1f spots/frame), %s -> %.2f (%.1f spots/frame)', ...
            prm.noiseK, St.C(1).key, prm.thr(1), yld(1), St.C(2).key, prm.thr(2), yld(2));
        % A colour that yields far less than the other is the thing that limits a CROSS-colour
        % analysis, because every pair needs one track from each. Say so here rather than letting it
        % turn up as an empty pair table three tabs later.
        if min(yld) > 0 && max(yld)/min(yld) > 3
            say(['  NOTE: %s yields %.0fx fewer spots than %s. Cross-colour pairing needs a track ' ...
                 'from each, so the weaker colour sets the ceiling — lower k and check the preview ' ...
                 'for false positives before accepting it.'], ...
                 St.C(argmin(yld)).key, max(yld)/max(min(yld),eps), St.C(argmax(yld)).key);
        elseif min(yld) < 0.5
            say('  NOTE: %s yields almost nothing at k = %g. Lower k, or this colour cannot be paired.', ...
                St.C(argmin(yld)).key, prm.noiseK);
        end
        detectPreview();
    end

    function runTracking()
        assert(~isempty(St.C), 'dc_app:noCells', 'pick a folder first');
        d = getappdata(fig,'detect');
        say('tracking both colours…');
        cel = struct('stack', St.stacks{1}, 'base', 'cell', 'stacks', struct());
        for i = 1:2, cel.stacks.(St.C(i).key) = St.stacks{i}; end
        P = struct('diamUm',prm.diamUm, 'thrAbs',prm.thr(1), 'pxUm',St.pxUm, ...
                   'linkUm',prm.linkUm, 'gapUm',prm.gapUm, 'maxGap',prm.maxGap);
        if ~isempty(prm.maxFrames), P.maxFrames = prm.maxFrames; end
        % A cell, not a struct array: assigning a filled struct into repmat(struct(),1,2) fails with
        % "dissimilar structures" the moment the two have different fields.
        R = cell(1,2);
        for i = 1:2
            Pi = P; Pi.thrAbs = prm.thr(i);       % per-colour threshold
            R{i} = dc_process_cell(cel, St.C(i), Pi);
            say('  %s: %d detections, %d tracks over %d frames', R{i}.key, R{i}.nDets, R{i}.nTracks, R{i}.nFrames);
        end
        D = dc_dataset('new', 'cell', St.C);
        for i = 1:2, D = dc_dataset('addChannel', D, R{i}); end
        St.D = D;
        say('dataset: %d spots, %d tracks', height(D.spots), numel(unique(D.spots.trackId(isfinite(D.spots.trackId)))));
    end

    function runComotion()
        assert(~isempty(St.D), 'dc_app:noTracks', 'track first (tab 2)');
        c = getappdata(fig,'co');
        S = dc_steps(St.D);
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
