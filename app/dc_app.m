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
% One open TIFF handle per colour, and one set of display limits per cell. imread(path,k) reopens
% the file and walks its directory chain to page k every single call — measured at 52 ms a page on a
% 6320-page stack, so 104 ms per frame before anything is drawn. Through a held handle the same read
% is 0.6 ms. The display limits are cached for a second reason as well as speed: recomputing
% percentiles per frame rescales the contrast to whatever is in that frame, so a molecule that
% bleaches looks constant because the picture keeps adjusting to it.
St = struct('folder','', 'stacks',{{}}, 'C',[], 'D',[], 'Draw',[], 'R',[], 'N',[], 'S',[], ...
            'rd',{{}}, 'rdClose',{{}}, 'disp',[], 'prep',[], ...
            'rejected',[], 'nDropped',0, ...
            'pxUm',0.0967821, 'dtS',NaN, 'pxSrc','default', 'dtSrc','default', 'align',[]);
% MANY CELLS, ONE AT A TIME. Each row of `cells` is a folder holding one cell's two colour stacks,
% with whatever has been computed for it so far. Settings are shared, results are per cell, and the
% pooled quality is too — a percentile of another cell's candidate distribution is not this cell's
% threshold. `iCell` says which one the tabs are showing.
cells = struct('folder',{}, 'name',{}, 'stacks',{}, 'C',[], 'D',[], 'Draw',[], 'R',[], 'N',[], ...
               'pool',{}, 'thr',[], 'rejected',[], 'pxUm',[], 'dtS',[], 'align',[], 'status',{});
iCell = 0;
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
% A player left running after the window goes keeps firing into deleted handles, once per tick,
% forever. Stop both before anything is torn down.
fig.CloseRequestFcn = @(s,e) onClose();
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
               'pool',@getPool, 'addCell',@addCell, 'selectCell',@selectCell, ...
               'cells',@getCells, 'runAll',@onRunAll);

if isfield(opts,'folder') && ~isempty(opts.folder), loadFolder(opts.folder); end

% ================================ tab 1: cells ====================================================
    function buildCells(parent)
        % The cell list, in the shape SPTinMatlab's Experiment tab has: everything the project holds,
        % what stage each one has reached, and which is selected. Settings are shared across cells;
        % results are not, and neither is the pooled quality — a percentile of another cell's
        % candidate distribution is not this cell's threshold.
        g = uigridlayout(parent,[4 1],'RowHeight',{28,'1x',76,26},'Padding',[10 10 10 10],'RowSpacing',7);

        hr = uigridlayout(g,[1 6],'ColumnWidth',{'1x',108,120,104,96,104},'Padding',[0 0 0 0],'ColumnSpacing',6);
        lblProj = uilabel(hr,'Text','No cells yet. Add one folder per cell, or scan a parent folder.', ...
            'FontColor',[0.45 0.45 0.5]);
        uibutton(hr,'Text','Add cell…','ButtonPushedFcn',@(s,e) onAddCell());
        uibutton(hr,'Text','Scan folder…','Tooltip', ...
            'Add every subfolder that holds two TIFF stacks, as one cell each.', ...
            'ButtonPushedFcn',@(s,e) onScan());
        uibutton(hr,'Text','Remove','ButtonPushedFcn',@(s,e) onRemoveCell());
        uibutton(hr,'Text','Run all','FontWeight','bold','Tooltip', ...
            'Pool, detect, track and run co-motion on every cell with the current settings.', ...
            'ButtonPushedFcn',@(s,e) onRunAll());
        uibutton(hr,'Text','Export all…','ButtonPushedFcn',@(s,e) onExportAll());

        tblCells = uitable(g,'ColumnName',{'cell','colours','frames','tracks','pairs','status'}, ...
            'ColumnWidth',{'auto',80,72,72,64,140},'RowName',{}, ...
            'SelectionType','row', 'CellSelectionCallback',@(s,e) onPickCell(e));

        ap = uigridlayout(g,[2 1],'RowHeight',{20,'1x'},'Padding',[0 0 0 0],'RowSpacing',3);
        uilabel(ap,'Text','The selected cell: what each page is, and how the two colours line up', ...
            'FontWeight','bold');
        lblAlign = uilabel(ap,'Text','—','WordWrap','on','FontColor',[0.2 0.2 0.25]);

        tblCh = uitable(g,'ColumnName',{'colour','stack','pages','frames','timepoints','dt (s)','source'}, ...
            'ColumnWidth',{62,'auto',64,64,96,80,88},'RowName',{});

        setappdata(fig,'cells', struct('tbl',tblCells,'ch',tblCh,'align',lblAlign,'proj',lblProj));
    end

    function onAddCell()
        f = uigetdir(tern(isempty(St.folder), pwd, fileparts(St.folder)), ...
            'Folder holding ONE cell''s two colour stacks');
        if isequal(f,0), return; end
        addCell(f); refreshCells(); selectCell(numel(cells));
    end

    function onScan()
        root = uigetdir(pwd, 'Parent folder — every subfolder with two stacks becomes a cell');
        if isequal(root,0), return; end
        d = dir(root); d = d([d.isdir] & ~startsWith({d.name},'.'));
        n0 = numel(cells);
        for k = 1:numel(d)
            p = fullfile(root, d(k).name);
            if numel([dir(fullfile(p,'*.tif')); dir(fullfile(p,'*.tiff'))]) >= 2, addCell(p); end
        end
        if numel([dir(fullfile(root,'*.tif')); dir(fullfile(root,'*.tiff'))]) >= 2, addCell(root); end
        refreshCells();
        say('scan added %d cell(s)', numel(cells) - n0);
        if numel(cells) > n0, selectCell(n0 + 1); end
    end

    function addCell(folder)
        if any(strcmp({cells.folder}, folder)), return; end     % adding twice is always a slip
        [~, nm] = fileparts(folder);
        cells(end+1) = struct('folder',folder, 'name',nm, 'stacks',{{}}, 'C',[], 'D',[], 'Draw',[], ...
            'R',[], 'N',[], 'pool',{{[] []}}, 'thr',prm.thr, 'rejected',[], ...
            'pxUm',[], 'dtS',[], 'align',[], 'status','added');
    end

    function onRemoveCell()
        c = getappdata(fig,'cells');
        k = c.tbl.Selection;
        if isempty(k), return; end
        cells(k(1)) = [];
        if iCell > numel(cells), iCell = numel(cells); end
        refreshCells();
        if iCell >= 1, selectCell(iCell); end
    end

    function onPickCell(e)
        if isempty(e.Indices), return; end
        selectCell(e.Indices(1));
    end

    function selectCell(k)
        if k < 1 || k > numel(cells), return; end
        stash();                       % keep what the outgoing cell has computed
        iCell = k;
        c = cells(k);
        St.folder = c.folder; St.stacks = c.stacks; St.C = c.C;
        St.D = c.D; St.Draw = c.Draw; St.R = c.R; St.N = c.N;
        St.rejected = c.rejected; St.align = c.align;
        if ~isempty(c.pxUm), St.pxUm = c.pxUm; end
        if ~isempty(c.dtS),  St.dtS  = c.dtS;  end
        pool = c.pool; prm.thr = c.thr;
        if isempty(St.C), loadFolder(c.folder); else, openReaders(); end
        refreshCells(); showChannels();
        eCalPx.Value = St.pxUm;
        if isfinite(St.dtS), eCalDt.Value = St.dtS; end
        detectPreview();
        if ~isempty(St.D), showTrackFrame(1); end
    end

    function stash()
        if iCell < 1 || iCell > numel(cells), return; end
        cells(iCell).stacks = St.stacks; cells(iCell).C = St.C;
        cells(iCell).D = St.D; cells(iCell).Draw = St.Draw;
        cells(iCell).R = St.R; cells(iCell).N = St.N;
        cells(iCell).rejected = St.rejected; cells(iCell).align = St.align;
        cells(iCell).pxUm = St.pxUm; cells(iCell).dtS = St.dtS;
        cells(iCell).pool = pool; cells(iCell).thr = prm.thr;
        cells(iCell).status = stageOf(iCell);
    end

    function s2 = stageOf(k)
        c = cells(k); s2 = 'added';
        if ~isempty(c.C), s2 = 'read'; end
        if ~isempty(c.pool) && ~isempty(c.pool{1}), s2 = 'pooled'; end
        if ~isempty(c.D), s2 = 'tracked'; end
        if ~isempty(c.R), s2 = 'co-motion'; end
    end

    function refreshCells()
        c = getappdata(fig,'cells');
        if isempty(c) || ~isgraphics(c.tbl), return; end
        stash();
        rows = cell(numel(cells), 6);
        for k = 1:numel(cells)
            q = cells(k);
            nCol = numel(q.C);
            nFr = ''; if ~isempty(q.C), nFr = sprintf('%d', q.C(1).nFrames); end
            nTr = ''; if ~isempty(q.D)
                nTr = sprintf('%d', numel(unique(q.D.spots.trackId(isfinite(q.D.spots.trackId))))); end
            nPr = ''; if ~isempty(q.R), nPr = sprintf('%d', height(q.R.pairs)); end
            rows(k,:) = {q.name, sprintf('%d', nCol), nFr, nTr, nPr, stageOf(k)};
        end
        c.tbl.Data = rows;
        if iCell >= 1 && iCell <= numel(cells), c.tbl.Selection = iCell; end
        c.proj.Text = sprintf('%d cell(s); settings are shared, results and pooled quality are per cell', ...
            numel(cells));
    end

    function showChannels()
        c = getappdata(fig,'cells');
        if isempty(St.C), c.ch.Data = {}; c.align.Text = '—'; return; end
        rows = cell(numel(St.C), 7);
        for i = 1:numel(St.C)
            [~,nm,ex] = fileparts(St.stacks{i});
            rows(i,:) = {St.C(i).key, [nm ex], sprintf('%d', numel(St.C(i).pages)), ...
                sprintf('%d', St.C(i).nFrames), ...
                sprintf('%d–%d', min(St.C(i).tp), max(St.C(i).tp)), ...
                sprintf('%.6g', St.C(i).dt_s), 'slice labels'};
        end
        c.ch.Data = rows;
        if ~isempty(St.align), c.align.Text = St.align.text; end
    end

    function onRunAll()
        assert(~isempty(cells), 'dc_app:noCells', 'add at least one cell');
        for k = 1:numel(cells)
            selectCell(k);
            say('=== %s (%d of %d) ===', cells(k).name, k, numel(cells));
            try
                poolQuality(); runTracking(); runComotion();
            catch ME
                say('  %s failed: %s', cells(k).name, ME.message);
            end
            stash();
        end
        refreshCells();
        say('run all finished over %d cell(s)', numel(cells));
    end

    function onExportAll()
        if isempty(cells), return; end
        f = uigetdir(pwd, 'Where to write one folder per cell');
        if isequal(f,0), return; end
        n = 0;
        for k = 1:numel(cells)
            if isempty(cells(k).R), continue; end
            d = fullfile(f, cells(k).name); if ~isfolder(d), mkdir(d); end
            writetable(cells(k).R.pairs, fullfile(d,'dc_comotion_pairs.csv'));
            writetable(cells(k).R.bins,  fullfile(d,'dc_comotion_bins.csv'));
            if ~isempty(cells(k).N), writetable(cells(k).N.curve, fullfile(d,'dc_comotion_null.csv')); end
            n = n + 1;
        end
        say('exported %d cell(s) to %s', n, f);
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

    function onClose()
        closeReaders();
        for nm = {'track','pair'}
            d = getappdata(fig, nm{1});
            if isstruct(d) && isfield(d,'timer') && ~isempty(d.timer) && isvalid(d.timer)
                stop(d.timer); delete(d.timer);
            end
        end
        delete(fig);
    end

    function openReaders()
        closeReaders();
        St.rd = cell(1,2); St.rdClose = cell(1,2);
        for i = 1:2
            [St.rd{i}, St.rdClose{i}] = dc_tiff_pages(St.stacks{i});
        end
        % Display limits from a sample of frames, then held: fixed across the movie so a change in
        % brightness is a change in the data and not in the scaling.
        loA=[]; hiA=[]; loB=[]; hiB=[];
        ts = round(linspace(1, numel(St.C(1).pages), min(6, numel(St.C(1).pages))));
        for t = ts
            a = double(St.rd{1}(St.C(1).pages(t)));
            b = double(St.rd{2}(St.C(2).pages(min(t, numel(St.C(2).pages)))));
            loA(end+1) = prctile(a(:),50); hiA(end+1) = prctile(a(:),99.9); %#ok<AGROW>
            loB(end+1) = prctile(b(:),50); hiB(end+1) = prctile(b(:),99.9); %#ok<AGROW>
        end
        St.disp = struct('loA',median(loA),'hiA',median(hiA), ...
                         'loB',median(loB),'hiB',median(hiB));
    end

    function closeReaders()
        for i = 1:numel(St.rdClose)
            try, if ~isempty(St.rdClose{i}), St.rdClose{i}(); end, catch, end
        end
        St.rd = {}; St.rdClose = {};
    end

    function im = readPg(i, pg)
        if numel(St.rd) >= i && ~isempty(St.rd{i}), im = double(St.rd{i}(pg));
        else, im = double(imread(St.stacks{i}, pg)); end
    end

    function P = drawSrc()
        if isempty(St.prep) && ~isempty(St.D), St.prep = dc_draw([], 'prep', St.D); end
        P = St.prep; if isempty(P), P = St.D; end
    end

    function q = getPool(), q = pool; end
    function c = getCells(), stash(); c = cells; end

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
            raw = readPg(i, pg);
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
        btnPlayT = uibutton(r1,'Text','▶ Play','ButtonPushedFcn',@(s,e) onPlay());

        r2 = uigridlayout(cp,[1 8],'ColumnWidth',{'1x',120,104,74,96,74,96,86}, ...
            'Padding',[0 0 0 0],'ColumnSpacing',6);
        sldT = uislider(r2,'Limits',[1 100],'Value',1,'MajorTicks',[], ...
            'ValueChangingFcn',@(s,e) showTrackFrame(round(e.Value)));
        % No view dropdown: all three are on screen at once. Each colour alone in grey is how a
        % detection is judged; the merge is how a pair is. Making them alternatives meant flipping
        % back and forth to answer two questions about the same frame.
        ddView = uilabel(r2,'Text','c1 · c2 · merged','FontColor',[0.45 0.45 0.5]);
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

        % Each axes lives in a uipanel, not directly in the grid. A uiaxes that is a grid child
        % resizes ITSELF to the data aspect ratio and will happily grow past its cell, ending up
        % underneath whatever is in the next row; inside a panel it letterboxes within a fixed box.
        row3 = uigridlayout(g,[1 3],'ColumnWidth',{'1x','1x','1x'},'Padding',[0 0 0 0],'ColumnSpacing',6);
        axT = gobjects(1,3);
        for i = 1:3
            pn = uipanel(row3,'BorderType','none');
            axT(i) = uiaxes(pn); axT(i).Units = 'normalized'; axT(i).Position = [0.08 0.07 0.87 0.85];
            axT(i).Toolbar.Visible = 'on';
        end

        lg = uigridlayout(g,[1 1],'Padding',[0 0 0 0]);
        txtT = uitextarea(lg,'Editable','off','Value',{''},'FontName','Menlo','FontSize',11);

        setappdata(fig,'track', struct('ax',axT,'sld',sldT,'view',ddView,'tail',spnTail, ...
            'gamma',spnGam,'paths',chkPaths,'lbl',lblTrk,'log',txtT,'tp',1,'timer',[], ...
            'play',btnPlayT));
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
        St.prep = dc_draw([], 'prep', St.D);     % the draw-ready view, rebuilt only when D changes
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
        cp = t.ax(3).CurrentPoint;      % the merged panel is where a click is aimed
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
        if isempty(St.D) || ~all(isgraphics(t.ax)), return; end   % t.ax is the three panels
        if isempty(tp), tp = t.tp; end
        tp = max(1, round(tp));
        nTp = max(St.C(1).tp);
        tp = min(tp, nTp);
        t.tp = tp; setappdata(fig,'track',t);

        pgA = pageFor(1, tp); pgB = pageFor(2, tp);
        if isempty(pgA) || isempty(pgB), return; end
        A = readPg(1, pgA);  B = readPg(2, pgB);
        if t.paths.Value, tail = t.tail.Value; else, tail = 1; end
        px = St.pxUm; rPx = (prm.diamUm/px)/2;
        gam = t.gamma.Value;

        [xl, yl] = keepView(t.ax(1), size(A));
        if isempty(xl), xl = [0.5 size(A,2)+0.5]; yl = [0.5 size(A,1)+0.5]; end
        names = {sprintf('%s', chKey(1)), sprintf('%s', chKey(2)), 'merged'};
        nA = 0; nB = 0;
        % Update in place rather than cla + image. Clearing an axes and building a new image and new
        % line objects every frame is what makes a MATLAB player stutter: the pixels are the cheap
        % part. The image handle is kept and its CData replaced; only the overlays, which genuinely
        % change shape, are rebuilt — and there are four of them, not hundreds, because dc_draw puts
        % every path in one NaN-separated line.
        if ~isfield(t,'him') || numel(t.him) ~= 3 || ~all(isgraphics(t.him))
            t.him = gobjects(1,3); t.hov = {[] [] []};
            for i = 1:3
                cla(t.ax(i));
                t.him(i) = image(t.ax(i), zeros(size(A,1), size(A,2), 3));
                set(t.ax(i),'DataAspectRatio',[1 1 1],'YDir','reverse');
                title(t.ax(i), names{i}, 'FontSize', 9);
                hold(t.ax(i),'on');
            end
        end
        for i = 1:3
            switch i
                case 1, im = gray3f(A, gam, dispLim('A'));
                case 2, im = gray3f(B, gam, dispLim('B'));
                case 3, im = dc_composite(A, B, dispOpts(gam));
            end
            set(t.him(i), 'CData', im);
            if ~isempty(t.hov{i}), delete(t.hov{i}(isgraphics(t.hov{i}))); end
            h = gobjects(0);
            if i == 1 || i == 3
                hA = dc_draw(t.ax(i),'tracks',drawSrc(),tp,struct('ch',chKey(1),'colour',[1 0.45 1], ...
                    'pxUm',px,'rPx',rPx,'tail',tail));
                nA = hA.n; h = [h hA.paths hA.now];
            end
            if i == 2 || i == 3
                hB = dc_draw(t.ax(i),'tracks',drawSrc(),tp,struct('ch',chKey(2),'colour',[0.45 1 0.5], ...
                    'pxUm',px,'rPx',rPx,'tail',tail));
                nB = hB.n; h = [h hB.paths hB.now];
            end
            t.hov{i} = h(isgraphics(h));
            set(t.ax(i), 'XLim', xl, 'YLim', yl);
        end
        setappdata(fig,'track',t);
        if isgraphics(t.sld), t.sld.Limits = [1 max(nTp,2)]; t.sld.Value = tp; end
        nDrop = 0; if isfield(St,'nDropped'), nDrop = St.nDropped; end
        t.lbl.Text = sprintf(['timepoint %d of %d · %d %s tracks and %d %s tracks visible · %d ' ...
            'excluded (shorter than %d, or rejected)'], tp, nTp, nA, chKey(1), nB, chKey(2), ...
            nDrop, prm.minLen);
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
            stop(t.timer); delete(t.timer); t.timer = [];
            t.play.Text = '▶ Play'; setappdata(fig,'track',t); return
        end
        t.timer = timer('ExecutionMode','fixedSpacing','Period',0.080,'TimerFcn',@(~,~) tick());
        t.play.Text = '❚❚ Pause'; setappdata(fig,'track',t); start(t.timer);
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
        % Two questions, in the order you ask them: IS there co-motion at short range, and WHICH
        % pairs. The verdict goes at the top in words; the two plots under it are the evidence for
        % it; the list is what you click.
        %
        % The old third plot was the null's SE against n on log axes. That is a diagnostic about the
        % method, not an answer about the data, and it sat where the answer should be. The
        % near-versus-far histogram replaces it: two distributions of the same quantity, one for
        % pairs that were close and one for pairs that were not. If the near one is shifted right,
        % there is co-motion, and you can see it without reading an axis.
        g = uigridlayout(parent,[4 1],'RowHeight',{30,52,'1x',72},'Padding',[10 10 10 10],'RowSpacing',7);

        r1 = uigridlayout(g,[1 12],'ColumnWidth',{84,56,96,56,84,56,78,56,108,'1x',122,84}, ...
            'Padding',[0 0 0 0],'ColumnSpacing',5);
        uilabel(r1,'Text','Near ≤ (µm)','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',co.rNearUm,'Limits',[0.02 20], ...
            'Tooltip','The cutoff that defines "nearby". Only pairs inside it are listed — they are the measurement.', ...
            'ValueChangedFcn',@(s,e) setCo('rNearUm',s.Value));
        uilabel(r1,'Text','Null starts at (µm)','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',co.rFarUm,'Limits',[0.05 20], ...
            'Tooltip','Pairs beyond this build the null: same movie, same drift, too far to interact.', ...
            'ValueChangedFcn',@(s,e) setCo('rFarUm',s.Value));
        uilabel(r1,'Text','Max r (µm)','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',co.rMaxUm,'Limits',[0.1 40], ...
            'Tooltip','Outer cap. The null is the shell between "null starts at" and this, so leave room.', ...
            'ValueChangedFcn',@(s,e) setCo('rMaxUm',s.Value));
        uilabel(r1,'Text','Min steps','HorizontalAlignment','right');
        uieditfield(r1,'numeric','Value',co.nMin,'Limits',[3 1e5],'RoundFractionalValues',true, ...
            'ValueChangedFcn',@(s,e) setCo('nMin',s.Value));
        chkDrift = uicheckbox(r1,'Text','Remove drift','Value',co.drift, ...
            'Tooltip',['Subtract the per-timepoint median step. Drift correlates EVERY pair at every ' ...
                       'separation and biases the mean, so more data makes it look MORE significant.'], ...
            'ValueChangedFcn',@(s,e) setCo('drift',s.Value));
        uilabel(r1,'Text','');
        uibutton(r1,'Text','Run co-motion','FontWeight','bold','ButtonPushedFcn',@(s,e) runComotion());
        uibutton(r1,'Text','Export…','ButtonPushedFcn',@(s,e) onExport());

        vp = uigridlayout(g,[2 1],'RowHeight',{26,22},'Padding',[0 0 0 0],'RowSpacing',2);
        lblVerdict = uilabel(vp,'Text','Run co-motion to see whether nearby tracks move together.', ...
            'FontSize',15,'FontWeight','bold','FontColor',[0.25 0.25 0.3]);
        lblCo = uilabel(vp,'Text','','FontColor',[0.4 0.4 0.45],'WordWrap','on');

        mid = uigridlayout(g,[1 3],'ColumnWidth',{'1.1x','0.85x','1.05x'},'Padding',[0 0 0 0],'ColumnSpacing',8);
        axC = uiaxes(mid);
        axH = uiaxes(mid);
        tblP = uitable(mid,'ColumnName',{'A','B','n','r (nm)','cos','z'}, ...
            'ColumnWidth',{54,54,50,62,58,52},'RowName',{}, ...
            'SelectionType','row','CellSelectionCallback',@(s,e) onPickPair(e));

        lg = uigridlayout(g,[1 1],'Padding',[0 0 0 0]);
        txtCo = uitextarea(lg,'Editable','off','Value',{''},'FontName','Menlo','FontSize',11);

        setappdata(fig,'co', struct('axC',axC,'axH',axH,'tbl',tblP,'lbl',lblCo, ...
            'verdict',lblVerdict,'log',txtCo,'drift',chkDrift));
    end

% ================================ tab 4: pair =====================================================
    function buildPair(parent)
        % THREE panels, not one. A magenta/green merge is the right picture for "are these two in
        % the same place", and the wrong one for "is this detection real": on a sparse, noisy frame
        % the merge turns background speckle into coloured confetti and a single molecule has to
        % compete with it. Each colour on its own, in grey, is how you judge the spot; the merge is
        % how you judge the pair. Showing all three costs nothing and each answers its own question.
        g = uigridlayout(parent,[3 1],'RowHeight',{30,28,'1x'},'Padding',[8 8 8 8],'RowSpacing',5);

        hr = uigridlayout(g,[1 9],'ColumnWidth',{'1x',70,70,68,94,62,78,92,142},'Padding',[0 0 0 0],'ColumnSpacing',5);
        lblPair = uilabel(hr,'Text','Select a pair in the Co-motion tab.','FontColor',[0.35 0.35 0.4]);
        uibutton(hr,'Text','◀ prev','ButtonPushedFcn',@(s,e) stepPair(-1));
        uibutton(hr,'Text','next ▶','ButtonPushedFcn',@(s,e) stepPair(+1));
        btnPlayP = uibutton(hr,'Text','▶ Play','ButtonPushedFcn',@(s,e) onPairPlay());
        uilabel(hr,'Text','Window (steps)','HorizontalAlignment','right');
        spnWin = uieditfield(hr,'numeric','Value',40,'Limits',[4 1e4],'RoundFractionalValues',true, ...
            'ValueChangedFcn',@(s,e) drawPairImage());
        chkQuiv = uicheckbox(hr,'Text','Quiver','Value',true, ...
            'Tooltip','Each step as an arrow from where the molecule was. This IS the co-motion.', ...
            'ValueChangedFcn',@(s,e) drawPairImage());
        chkZoom = uicheckbox(hr,'Text','Zoom to pair','Value',true, ...
            'ValueChangedFcn',@(s,e) drawPairImage());
        ddView4 = uidropdown(hr,'Items',{'centred paths','step rose','unit rose','shared vs relative'}, ...
            'ItemsData',{'centred','rose','unitrose','decompose'},'Value','centred', ...
            'Tooltip',['What the fourth panel shows. All of them take the space between the two ' ...
                       'molecules out, so their motion can be compared directly instead of across ' ...
                       'the gap between them.'], ...
            'ValueChangedFcn',@(s,e) drawPairImage());

        sr = uigridlayout(g,[1 3],'ColumnWidth',{'1x',150,110},'Padding',[0 0 0 0],'ColumnSpacing',6);
        sldP = uislider(sr,'Limits',[1 100],'Value',1,'MajorTicks',[], ...
            'ValueChangingFcn',@(s,e) onPairSlide(e.Value));
        lblTp = uilabel(sr,'Text','—','FontName','Menlo','FontSize',11);
        spnFps = uieditfield(sr,'numeric','Value',12,'Limits',[1 60], ...
            'Tooltip','Playback rate, frames per second.');

        % Images left, traces right. The three image panels share the left half: each colour on its
        % own along the top, the merge beneath them spanning both — the merge is the one you read
        % the pair off, so it gets the wider box.
        % Images left, traces right. Top row: each colour alone, in grey. Bottom row: the merge
        % with the quiver in place, and beside it the same motion with the space between the two
        % taken out — which is where "do they move together" is actually legible.
        mn = uigridlayout(g,[1 2],'ColumnWidth',{'1.25x','1x'},'Padding',[0 0 0 0],'ColumnSpacing',8);
        imgs = uigridlayout(mn,[2 2],'RowHeight',{'1x','1.1x'},'ColumnWidth',{'1x','1x'}, ...
            'Padding',[0 0 0 0],'RowSpacing',4,'ColumnSpacing',4);
        axP = gobjects(1,4);
        for i = 1:4
            pn = uipanel(imgs,'BorderType','none');
            axP(i) = uiaxes(pn); axP(i).Units='normalized'; axP(i).Position=[0.13 0.11 0.82 0.78];
        end
        host = uipanel(mn,'BorderType','none');       % the traces get the whole right half

        setappdata(fig,'pair', struct('host',host,'ax',axP,'lbl',lblPair,'k',0,'tp',1, ...
            'win',spnWin,'quiv',chkQuiv,'zoom',chkZoom,'sld',sldP,'tpLbl',lblTp, ...
            'play',btnPlayP,'fps',spnFps,'timer',[],'view4',ddView4));
    end

    function onPairSlide(v)
        pp = getappdata(fig,'pair'); pp.tp = round(v); setappdata(fig,'pair',pp);
        drawPairImage();
    end

    function onPairPlay()
        pp = getappdata(fig,'pair');
        if ~isempty(pp.timer) && isvalid(pp.timer)
            stop(pp.timer); delete(pp.timer); pp.timer = [];
            pp.play.Text = '▶ Play'; setappdata(fig,'pair',pp); return
        end
        % Rounded to the millisecond: timer refuses finer precision and warns on every construction,
        % and 1/12 of a second is not a whole number of them.
        per = round(max(1/pp.fps.Value, 0.03), 3);
        pp.timer = timer('ExecutionMode','fixedSpacing','Period',per, 'TimerFcn',@(~,~) pairTick());
        pp.play.Text = '❚❚ Pause'; setappdata(fig,'pair',pp); start(pp.timer);
    end
    function pairTick()
        pp = getappdata(fig,'pair');
        if ~isgraphics(fig) || isempty(St.R) || pp.k < 1
            if ~isempty(pp.timer) && isvalid(pp.timer), stop(pp.timer); delete(pp.timer); end
            return
        end
        tps = pairTps(pp.k);
        if isempty(tps), return; end
        nxt = pp.tp + 1; if nxt > max(tps), nxt = min(tps); end
        pp.tp = nxt; setappdata(fig,'pair',pp);
        drawPairImage();
    end
    function tps = pairTps(k)
        tps = [];
        if isempty(St.R) || k < 1 || k > height(St.R.pairs), return; end
        m = St.R.steps.trackA == St.R.pairs.trackA(k) & St.R.steps.trackB == St.R.pairs.trackB(k);
        tps = sort(St.R.steps.tp(m));
    end

    function drawPairImage()
        pp = getappdata(fig,'pair');
        if isempty(St.R) || pp.k < 1 || pp.k > height(St.R.pairs) || ~all(isgraphics(pp.ax)), return; end
        row = St.R.pairs(pp.k,:);
        idA = row.trackA; idB = row.trackB;

        st = St.R.steps(St.R.steps.trackA == idA & St.R.steps.trackB == idB, :);
        if isempty(st), return; end
        st = sortrows(st,'tp');
        tps = st.tp;
        tpNow = min(max(pp.tp, min(tps)), max(tps));
        w = pp.win.Value;
        sw = st(tps >= tpNow - w/2 & tps <= tpNow + w/2, :);

        pgA = pageFor(1, tpNow); pgB = pageFor(2, tpNow);
        if isempty(pgA) || isempty(pgB), return; end
        A = readPg(1, pgA);  B = readPg(2, pgB);
        px = St.pxUm; rPx = (prm.diamUm/px)/2;

        % the box the three panels share, so they are comparable at a glance
        xl = [0.5 size(A,2)+0.5]; yl = [0.5 size(A,1)+0.5];
        if pp.zoom.Value && ~isempty(sw)
            pad = max(0.8/px, 3*rPx);
            xs = [sw.xa; sw.xb]/px; ys = [sw.ya; sw.yb]/px;
            xl = [min(xs)-pad, max(xs)+pad]; yl = [min(ys)-pad, max(ys)+pad];
        end

        names = {sprintf('%s only', chKey(1)), sprintf('%s only', chKey(2)), 'merged'};
        cols  = {[1 0.45 1], [0.45 1 0.5], []};
        for i = 1:3
            ax = pp.ax(i); cla(ax);
            switch i
                case 1, im = gray3f(A, 0.8, dispLim('A'));    % each colour in GREY on its own:
                case 2, im = gray3f(B, 0.8, dispLim('B'));    % judging a spot should not fight a hue
                case 3, im = dc_composite(A, B, dispOpts(0.8));
            end
            image(ax, im);
            set(ax,'DataAspectRatio',[1 1 1],'YDir','reverse','XLim',xl,'YLim',yl);
            hold(ax,'on');
            if i == 1 || i == 3
                dc_draw(ax,'tracks',drawSrc(),tpNow,struct('ch',chKey(1),'colour',cols{1}, ...
                    'pxUm',px,'rPx',rPx,'tail',w,'tracks',idA));
            end
            if i == 2 || i == 3
                dc_draw(ax,'tracks',drawSrc(),tpNow,struct('ch',chKey(2),'colour',cols{2}, ...
                    'pxUm',px,'rPx',rPx,'tail',w,'tracks',idB));
            end
            if i == 3
                hp = dc_draw(ax,'pair',drawSrc(),idA,idB,tpNow, ...
                    struct('pxUm',px,'rPx',rPx,'tail',w));
                if pp.quiv.Value && ~isempty(sw)
                    % Scale 0, autoscale off: the arrows are the real displacements in the image's
                    % own units. MATLAB would otherwise resize them to look tidy, which makes two
                    % arrows of very different length look comparable.
                    quiver(ax, sw.xa/px, sw.ya/px, sw.uax/px, sw.uay/px, 0, ...
                        'Color',cols{1}, 'LineWidth',0.9, 'MaxHeadSize',0.5);
                    quiver(ax, sw.xb/px, sw.yb/px, sw.ubx/px, sw.uby/px, 0, ...
                        'Color',cols{2}, 'LineWidth',0.9, 'MaxHeadSize',0.5);
                end
            end
            hold(ax,'off');
            title(ax, names{i}, 'FontSize', 9);
        end

        % the fourth panel: the same pair with the offset between them removed
        m4 = pp.view4.Value;
        vi = dc_comotion_view(pp.ax(4), sw, strrep(m4,'unitrose','rose'), ...
            struct('unit', strcmp(m4,'unitrose')));
        St.lastView = vi;

        if isgraphics(pp.sld), pp.sld.Limits = [min(tps) max(tps)]; pp.sld.Value = tpNow; end
        pp.tp = tpNow; setappdata(fig,'pair',pp);
        rNow = NaN;
        j = find(st.tp == tpNow, 1); if ~isempty(j), rNow = 1000*st.r(j); end
        pp.tpLbl.Text = sprintf('tp %d · %.0f nm', tpNow, rNow);
    end

    function o = dispOpts(gam)
        o = struct('gamma', gam);
        if ~isempty(St.disp)
            o.loA = St.disp.loA; o.hiA = St.disp.hiA;
            o.loB = St.disp.loB; o.hiB = St.disp.hiB;
        end
    end
    function lh = dispLim(which)
        lh = [];
        if isempty(St.disp), return; end
        if which == 'A', lh = [St.disp.loA St.disp.hiA]; else, lh = [St.disp.loB St.disp.hiB]; end
    end
    function g3 = gray3f(X, gam, lh)
        % Fixed limits when the cell has them, so brightness on screen tracks brightness in the data.
        if isempty(lh), lh = [prctile(X(:),50) prctile(X(:),99.9)]; end
        lo = lh(1); hi = lh(2); if ~(hi>lo), hi = lo+1; end
        y = min(max((X-lo)/(hi-lo),0),1) .^ gam;
        g3 = cat(3,y,y,y);
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

        end
        St.C = dc_channels('manual', specs);
        if any(isfinite(dtAll)), St.dtS = median(dtAll(isfinite(dtAll))); St.dtSrc = 'metadata'; end

        openReaders();
        A = dc_align(St.C(1), St.C(2));
        St.align = A;
        if iCell < 1
            addCell(St.folder); iCell = numel(cells);
        end
        stash(); refreshCells(); showChannels();
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
        % nearUm is the measurement's cutoff: the pair list is the pairs that were close enough to
        % be worth looking at. Steps at every separation are still collected, because the far ones
        % are the null and the thing that says how many steps a correlation needs.
        R = dc_comotion(S, struct('rMaxUm',co.rMaxUm, 'nMin',co.nMin, ...
                                  'nearUm',co.rNearUm, 'classes',"cross"));
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

        % --- the curve, with the null as a BAND rather than a line -------------------------------
        cla(c.axC); hold(c.axC,'on');
        if ~isempty(N)
            % A band, because "is this point above the null" is the only question being asked of
            % this plot, and a dashed line makes the reader do the comparison by eye against a
            % number that has its own uncertainty.
            sdFar = N.far.sdCos / sqrt(max(N.far.nStepPairs,1));
            xb = [0 max(b.rHi)*1000];
            fill(c.axC, [xb fliplr(xb)], ...
                 [N.far.meanCos-2*sdFar, N.far.meanCos-2*sdFar, ...
                  N.far.meanCos+2*sdFar, N.far.meanCos+2*sdFar], ...
                 [0.85 0.88 0.92], 'EdgeColor','none', 'DisplayName','far-field null');
        end
        errorbar(c.axC, b.rMid*1000, b.meanCos, 2*b.seCos, 'o-', 'Color',[0.80 0.20 0.20], ...
            'MarkerFaceColor',[0.80 0.20 0.20], 'MarkerSize',4, 'LineWidth',1.3, ...
            'DisplayName','cross-colour');
        xline(c.axC, co.rNearUm*1000, ':', 'Color',[0.3 0.5 0.3], 'HandleVisibility','off');
        yline(c.axC, 0, '-', 'Color',[0.8 0.8 0.8], 'HandleVisibility','off');
        hold(c.axC,'off'); grid(c.axC,'on');
        % Show the near field and the start of the null, not the whole outer cap: rMaxUm exists to
        % collect distant pairs for the null, and plotting out to it squeezes everything that
        % matters into the first few pixels.
        xlim(c.axC, [0 min(co.rMaxUm, 2*co.rFarUm)*1000]);
        ylim(c.axC, 'auto');
        xlabel(c.axC,'separation (nm)'); ylabel(c.axC,'mean cos\theta');
        title(c.axC,'co-motion vs separation','FontSize',10);
        legend(c.axC,'Location','northeast','Box','off','FontSize',8);

        % --- near vs far, as distributions of the SAME quantity -----------------------------------
        cla(c.axH); hold(c.axH,'on');
        near = R.steps.cos(R.steps.r <= co.rNearUm);
        far  = R.steps.cos(R.steps.r >  co.rFarUm);
        ed = -1:0.1:1;
        if ~isempty(far)
            histogram(c.axH, far, ed, 'Normalization','probability', ...
                'FaceColor',[0.62 0.66 0.70], 'EdgeColor','none', 'DisplayName','far (null)');
        end
        if ~isempty(near)
            histogram(c.axH, near, ed, 'Normalization','probability', ...
                'FaceColor',[0.80 0.25 0.25], 'FaceAlpha',0.65, 'EdgeColor','none', ...
                'DisplayName',sprintf('near (\\leq%.2f µm)', co.rNearUm));
        end
        if ~isempty(far),  xline(c.axH, mean(far,'omitnan'),  '-', 'Color',[0.35 0.4 0.45], 'HandleVisibility','off'); end
        if ~isempty(near), xline(c.axH, mean(near,'omitnan'), '-', 'Color',[0.75 0.15 0.15], 'LineWidth',1.4, 'HandleVisibility','off'); end
        hold(c.axH,'off'); grid(c.axH,'on');
        xlabel(c.axH,'cos\theta per step'); ylabel(c.axH,'fraction');
        title(c.axH,'near vs far — shifted right means co-motion','FontSize',10);
        legend(c.axH,'Location','northwest','Box','off','FontSize',8);

        P = R.pairs;
        % Formatted as text: a uitable renders a double as "8.0000" and then truncates it, so a
        % column of track ids reads as measurements with four decimal places.
        c.tbl.Data = [cellstr(string(P.trackA)), cellstr(string(P.trackB)), ...
                      cellstr(string(P.n)), cellstr(compose('%.0f', P.rMedian*1000)), ...
                      cellstr(compose('%+.3f', P.meanCos)), cellstr(compose('%.1f', P.z))];

        % --- the verdict, in words -----------------------------------------------------------------
        if isempty(N)
            c.verdict.Text = sprintf('%d pairs within %.2f µm — no null could be built', ...
                height(P), co.rNearUm);
            c.verdict.FontColor = [0.45 0.45 0.5];
        else
            nSE = N.far.sdCos / sqrt(max(nnz(R.steps.r <= co.rNearUm),1));
            zz = N.excess / max(nSE, eps);
            if zz >= 3
                c.verdict.Text = sprintf(['Nearby tracks DO move together: cos %+.3f within %.2f µm ' ...
                    'against %+.3f far away (%.0f sigma).'], N.near.meanCos, co.rNearUm, ...
                    N.far.meanCos, zz);
                c.verdict.FontColor = [0.65 0.13 0.13];
            else
                c.verdict.Text = sprintf(['No co-motion above the null: cos %+.3f within %.2f µm ' ...
                    'against %+.3f far away (%.1f sigma).'], N.near.meanCos, co.rNearUm, ...
                    N.far.meanCos, zz);
                c.verdict.FontColor = [0.3 0.35 0.4];
            end
            nNeed = N.nNeeded.steps(end);
            c.lbl.Text = sprintf(['%d pairs listed (within %.2f µm, %d+ shared steps). A pair needs ' ...
                'about %d steps to resolve a correlation of %.2f, so %d of them are long enough.'], ...
                height(P), co.rNearUm, co.nMin, nNeed, N.nNeeded.correlation(end), nnz(P.n >= nNeed));
        end
        if ~isempty(N), sayCo('%s', N.text); end
        if ~isempty(P), showPair(1); end
    end

    function showPair(k)
        if isempty(St.R) || isempty(St.R.pairs) || k < 1 || k > height(St.R.pairs), return; end
        pp = getappdata(fig,'pair'); pp.k = k;
        st = St.R.steps(St.R.steps.trackA == St.R.pairs.trackA(k) & ...
                        St.R.steps.trackB == St.R.pairs.trackB(k), :);
        if ~isempty(st), pp.tp = round(median(st.tp)); end
        setappdata(fig,'pair',pp);
        dtv = []; if isfinite(St.dtS), dtv = St.dtS; end
        Hp = dc_pair_panel(St.D, St.R, k, struct('parent',pp.host,'null',St.N,'dtS',dtv));
        drawPairImage();
        pp = getappdata(fig,'pair');
        sh = '';
        if isfield(St,'lastView') && ~isempty(St.lastView) && isfinite(St.lastView.shared)
            % `shared` is the correlation itself; the cosine estimates (pi/4) times it. Both are
            % shown because the cosine is what the pair list is ranked by.
            sh = sprintf(', shared %+.2f', St.lastView.shared);
        end
        pp.lbl.Text = sprintf(['pair %d of %d — tracks %g and %g, %d steps, cos %+.3f%s, ' ...
            '%.0f nm apart'], k, height(St.R.pairs), Hp.trackA, Hp.trackB, Hp.n, Hp.meanCos, ...
            sh, 1000*Hp.rMedian);
        tg.SelectedTab = t5;      % the Pair tab — t4 is Co-motion since the Track tab arrived
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
