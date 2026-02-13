function plot_freevolume_EOL_AOS_2x4_unique()
% Spatial free-volume maps for two conditions (rows) × four alloys per row (columns)
% Row 1 = EOL.xlsx (4 alloys), Row 2 = AOS.xlsx (4 different alloys)
% Modes: 'raw' (default), 'zscore' (baseline), 'delta' (kept for parity)

%% ---------------- USER SETTINGS ----------------
% Workbooks by row
preFiles   = {'EOL.xlsx','AOS.xlsx'};         % row 1, row 2
rowLabels  = {'EOL','AOS'};

% Sheets per row (set the exact sheet names in each workbook)
rowSheetNames = { ...
    {'Sheet1','Sheet2','Sheet3','Sheet4'}, ... % EOL's 4 alloys (columns)
    {'Sheet1','Sheet2','Sheet3','Sheet4'}  ... % AOS's 4 (different) alloys
};

% Alloy labels per row (LaTeX-friendly; must match the sheet order above)
alloyLblRows = { ...
    { 'Cu_{50}Zr_{50}', ...
      'Cu_{45}Zr_{48}Al_{7}', ...
      'Cu_{39}Zr_{46}Al_{7}Ag_{8}', ...
      'Cu_{37}Zr_{42}Al_{7}Ag_{8}Ti_{6}' }, ...   % EOL row
    { 'Cu_{45}Zr_{47}Ag_{8}', ...
      'Cu_{48}Zr_{38}Ag_{8}Ti_{6}', ...
      'Cu_{50}Zr_{44}Ti_{6}', ...
      'Cu_{40}Zr_{47}Ti_{6}Ag_{7}' }              % AOS row
};

% Free-volume column auto-detection keys
valueKeys  = {'free_volume','free','vfree','voronoi','volume','v_free','voronoi_volume'};

% Plot mode (as-cast): use 'raw' by default
plotMode   = 'raw';  % 'raw' | 'zscore' | 'delta' | 'auto'

% Orientation & slice
plane      = 'YZ';     % 'XZ','XY','YZ'
sliceFrac  = 0.40;     % thickness fraction along dropped axis
gridN      = 400;      % per axis
smooth_px  = 0.3;      % Gaussian blur (pixels), 0 = off

% (Only for zscore far-field baseline)
R_tip  = 30;   % Å
h_max  = 20;   % Å
far_k_list = [3.0 2.5 2.0 1.8 1.5 1.3 1.2 1.0];

% Color scaling
symPercentile = 95;     % for diverging (delta/zscore)
seqPercentile = [5 95]; % for sequential (raw)

% Palettes
divergingPalette  = 'coolwarm';
sequentialPalette = 'parula';
%% ------------------------------------------------

% Basic checks
assert(numel(preFiles)==2 && numel(rowLabels)==2, 'Expect exactly 2 rows (EOL, AOS).');
assert(numel(rowSheetNames)==2 && numel(alloyLblRows)==2, 'Provide sheets and labels for both rows.');
assert(numel(rowSheetNames{1})==4 && numel(rowSheetNames{2})==4, 'Each row must have 4 sheets.');
assert(numel(alloyLblRows{1})==4 && numel(alloyLblRows{2})==4, 'Each row must have 4 alloy labels.');

if strcmpi(plotMode,'auto')
    % No post files here → use zscore for comparability
    plotMode = 'zscore';
end

nRows = 2; nCols = 4;

% ---- load all datasets; track global bounds ----
D = cell(nCols, nRows);
xlimG=[inf,-inf]; ylimG=[inf,-inf]; zlimG=[inf,-inf];

for r = 1:nRows          % rows: EOL then AOS
  for c = 1:nCols        % cols: 4 alloys per row
    Tpre = readtable(preFiles{r}, "Sheet", rowSheetNames{r}{c});
    x = Tpre.(pickCol(Tpre,'x')); 
    y = Tpre.(pickCol(Tpre,'y')); 
    z = Tpre.(pickCol(Tpre,'z'));
    v = Tpre.(pickValCol(Tpre,valueKeys));
    m = isfinite(x)&isfinite(y)&isfinite(z)&isfinite(v);
    x=x(m); y=y(m); z=z(m); v=v(m);
    idpre = []; if hasCol(Tpre,'id'), idpre = Tpre.(pickCol(Tpre,'id')); idpre = idpre(m); end

    S = struct('x',x,'y',y,'z',z,'vpre',v,'vpost',[], ...
               'hasPost',false,'idpre',idpre,'idpost',[]);

    modeNow = plotMode;
    if strcmpi(plotMode,'delta')
        warning('No POST files configured; switching %s/%s to zscore.', rowLabels{r}, rowSheetNames{r}{c});
        modeNow = 'zscore';
    end

    D{c,r} = struct('S',S,'modeNow',modeNow);
    xlimG = [min(xlimG(1),min(S.x)) max(xlimG(2),max(S.x))];
    ylimG = [min(ylimG(1),min(S.y)) max(ylimG(2),max(S.y))];
    zlimG = [min(zlimG(1),min(S.z)) max(zlimG(2),max(S.z))];
  end
end

cx = mean(xlimG); cy = mean(ylimG); cz = mean(zlimG);
Lx = diff(xlimG); Ly = diff(ylimG); Lz = diff(zlimG);
a_contact = sqrt(max(0, 2*R_tip*h_max - h_max^2));

% build common grid
switch upper(plane)
  case 'XZ'
    xi = linspace(xlimG(1),xlimG(2),gridN);
    zi = linspace(zlimG(1),zlimG(2),gridN);
    [XI,ZI] = meshgrid(xi,zi);
  case 'XY'
    xi = linspace(xlimG(1),xlimG(2),gridN);
    yi = linspace(ylimG(1),ylimG(2),gridN);
    [XI,ZI] = meshgrid(xi,yi); % ZI used as Y
  case 'YZ'
    yi = linspace(ylimG(1),ylimG(2),gridN);
    zi = linspace(zlimG(1),zlimG(2),gridN);
    [XI,ZI] = meshgrid(yi,zi); % XI=Y, ZI=Z
  otherwise
    error('plane must be XZ, XY, or YZ');
end

% compute panel fields
valsAll = [];
maps = cell(nCols, nRows);
labelForCB = '';

for r = 1:nRows
  for c = 1:nCols
    S = D{c,r}.S;
    modeNow = D{c,r}.modeNow;

    % slice mask
    switch upper(plane)
      case 'XZ', half = 0.5*sliceFrac*Ly; in = S.y>=cy-half & S.y<=cy+half;
      case 'XY', half = 0.5*sliceFrac*Lz; in = S.z>=cz-half & S.z<=cz+half;
      case 'YZ', half = 0.5*sliceFrac*Lx; in = S.x>=cx-half & S.x<=cx+half;
    end

    switch lower(modeNow)
      case 'delta'
        val = S.vpost - S.vpre;
        labelForCB = '\DeltaV_{free} (Å^3)';
      case 'zscore'
        [mu, sg, rule] = chooseBaseline(S, cx, cy, a_contact, far_k_list);
        val = (S.vpre - mu) ./ max(sg,eps);
        labelForCB = sprintf('z(V_{free}) — %s', rule);
      case 'raw'
        val = S.vpre;
        labelForCB = 'V_{free} (Å^3)';
      otherwise
        error('Unknown plotMode: %s', modeNow);
    end

    % interpolate onto grid
    switch upper(plane)
      case 'XZ'
        F = scatteredInterpolant(S.x(in),S.z(in),val(in),'natural','nearest');
        V = F(XI,ZI);
      case 'XY'
        F = scatteredInterpolant(S.x(in),S.y(in),val(in),'natural','nearest');
        V = F(XI,ZI);
      case 'YZ'
        F = scatteredInterpolant(S.y(in),S.z(in),val(in),'natural','nearest');
        V = F(XI,ZI);
    end
    if smooth_px>0, V = gaussBlur2D(V, smooth_px); end
    maps{c,r} = V;
    valsAll = [valsAll; V(:)]; %#ok<AGROW>
  end
end

% color limits
if any(strcmpi(plotMode,{'delta','zscore'}))
  absv = abs(valsAll(isfinite(valsAll)));
  v = 1; if ~isempty(absv), v = prctile(absv, symPercentile); end
  clims = [-v, +v];
  palette = divergingPalette;
else
  v1 = prctile(valsAll, seqPercentile(1));
  v2 = prctile(valsAll, seqPercentile(2));
  if ~isfinite(v1) || ~isfinite(v2) || v1==v2
      v1 = min(valsAll); v2 = max(valsAll);
  end
  clims = [v1, v2];
  palette = sequentialPalette;
end

% plot 2×4
f = figure('Color','w','Position',[100 100 1600 800]);
tiledlayout(nRows,nCols,'Padding','compact','TileSpacing','compact');

for r = 1:nRows
  for c = 1:nCols
    nexttile
    imagesc(maps{c,r}); set(gca,'YDir','normal'); axis image
    switch lower(palette)
      case 'coolwarm', colormap(coolwarm(256));
      case 'parula',   colormap(parula(256));
      case 'turbo',    colormap(turbo(256));
      otherwise,       colormap(parula(256));
    end
    caxis(clims);
    switch upper(plane)
      case 'XZ', xlabel('x (Å)'); ylabel('z (Å)');
      case 'XY', xlabel('x (Å)'); ylabel('y (Å)');
      case 'YZ', xlabel('y (Å)'); ylabel('z (Å)');
    end
    title(sprintf('%s — %s', alloyLblRows{r}{c}, rowLabels{r}), 'FontWeight','normal')
    set(gca,'LineWidth',1,'FontSize',10)

    if r==nRows && c==nCols
      cb = colorbar; cb.Label.String = labelForCB; cb.Box='off';
    end
  end
end

outfile = sprintf('FreeVolume_%s_2x4_%s.png', upper(plane), lower(palette));
print(gcf, outfile, '-dpng','-r300');
fprintf('Saved: %s\n', outfile);

end % main

%% ----------------- helpers (same as before) -----------------
function name = pickCol(T, prefer)
names = T.Properties.VariableNames; low = lower(names);
p = lower(prefer);
idx = find(strcmp(low,p),1);
if ~isempty(idx), name = names{idx}; return; end
idx = find(contains(low,p),1);
if ~isempty(idx), name = names{idx}; return; end
error('Column "%s" not found. Available: %s', prefer, strjoin(names, ', '));
end

function name = pickValCol(T, keys)
names = T.Properties.VariableNames; low = lower(names);
for k = 1:numel(keys)
    hit = find(contains(low, lower(keys{k})),1);
    if ~isempty(hit), name = names{hit}; return; end
end
error('Could not find a free-volume column by keys: %s', strjoin(keys,', '));
end

function tf = hasCol(T, prefer)
tf = any(strcmpi(T.Properties.VariableNames, prefer)) || ...
     any(contains(lower(T.Properties.VariableNames), lower(prefer)));
end

function M = gaussBlur2D(M, sigma)
if sigma<=0, return; end
sz = max(3, round(6*sigma)); if mod(sz,2)==0, sz=sz+1; end
x = linspace(-((sz-1)/2), ((sz-1)/2), sz);
g = exp(-0.5*(x/sigma).^2); g = g/sum(g);
M = conv2(conv2(M,g,'same'), g','same');
end

function [mu, sg, rule] = chooseBaseline(S, cx, cy, a_contact, klist)
rr = hypot(S.x - cx, S.y - cy);
mu = NaN; sg = NaN; rule = '';
for k = klist
    far = rr > k * a_contact;
    ref = S.vpre(far);
    if nnz(isfinite(ref)) >= 0.10 * numel(S.vpre)
        mu = mean(ref,'omitnan');
        sg = std(ref,'omitnan');
        rule = sprintf('far-field (k=%.1f)', k);
        break
    end
end
if ~isfinite(mu) || ~isfinite(sg) || sg == 0
    med = median(S.vpre,'omitnan');
    mad = 1.4826 * median(abs(S.vpre - med),'omitnan');
    mu = med;
    sg = max(mad, eps);
    rule = 'global robust (median/MAD)';
end
end

function cmap = coolwarm(n)
if nargin<1, n=256; end
anchors = [...
    0.230, 0.299, 0.754;
    0.472, 0.621, 0.871;
    0.780, 0.914, 0.973;
    0.968, 0.968, 0.968;
    0.992, 0.839, 0.733;
    0.875, 0.512, 0.318;
    0.706, 0.016, 0.150];
t0 = linspace(0,1,size(anchors,1))';
t  = linspace(0,1,n)';
cmap = [interp1(t0,anchors(:,1),t,'pchip'), ...
        interp1(t0,anchors(:,2),t,'pchip'), ...
        interp1(t0,anchors(:,3),t,'pchip')];
cmap = max(0,min(1,cmap));
end