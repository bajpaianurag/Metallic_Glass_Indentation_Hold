function plot_D2min_four_alloys_three_states()
% Plot 4×3 D^2_min maps for planes XZ, XY, YZ
% Alloys (rows): CuZr, CuZrAl, CuZrAlAg, CuZrAlAgTi
% States (cols): as-cast, 200C, 400C
% Also computes robust STZ metrics via 3D voxel clustering (once per panel).

%% ---------------- USER SETTINGS ----------------
% Excel files (one per alloy; each has 3 sheets = thermal states)
alloyFiles = { ...
    'CuZr.xlsx', ...         % Cu50Zr50
    'CuZrAl.xlsx', ...       % Cu47.5Zr47.5Al5 (example)
    'CuZrAlAg.xlsx', ...     % add your file
    'CuZrAlAgTi.xlsx' ...    % add your file
};

% Sheet names inside each workbook (order must match states below)
stateSheets = {'Sheet1','Sheet2','Sheet3'};   % exact sheet names recommended

% Alloy display names (row labels)
alloyNames = {'Cu_{50}Zr_{50}', 'Cu_{45}Zr_{48}Al_{7}', ...
              'Cu_{39}Zr_{46}Al_{7}Al_{8}', 'Cu_{37}Zr_{42}Al_{7}Ag_{8}Ti_{6}'};

% Indenter / contact geometry (for far-field baseline only)
R_tip  = 30;    % Å (3 nm)
h_max  = 20;    % Å (2 nm)
far_k  = 3.0;   % far-field radius r > k * a_contact, with a = sqrt(2Rh - h^2)

% Slice & grid for plane maps
sliceFrac = 0.30;   % slab thickness as fraction of box length
gridN     = 400;    % grid resolution (per axis)
sigma_px  = 0.3;    % Gaussian blur (pixels) to de-speckle

usePlasma = true;

% Robust color scaling (per figure: all 12 panels share limits)
prcLo = 1; prcHi = 99;

% Thresholding options for STZ (priority order)
useFixedThreshold = false;   fixedThr = 0.02;   % Å^2, if useFixedThreshold = true
use_P95_threshold = false;   % true => P95 of far-field; false => μ+3σ

% 3D clustering (voxelization)
voxel_h         = 2.5;     % Å
min_vox_cluster = 10;      % drop clusters < 10 voxels (noise)

%% ---------------- LOAD DATA ----------------
nAlloys = numel(alloyFiles);
nStates = numel(stateSheets);

D = cell(nAlloys, nStates);
xlims = [inf,-inf]; ylims = [inf,-inf]; zlims = [inf,-inf];

for iA = 1:nAlloys
    for jS = 1:nStates
        T = readtable(alloyFiles{iA}, "Sheet", stateSheets{jS});
        x  = T.(pickCol(T,'x'));
        y  = T.(pickCol(T,'y'));
        z  = T.(pickCol(T,'z'));
        d2 = T.(pickCol(T,'D2min'));
        m = isfinite(x) & isfinite(y) & isfinite(z) & isfinite(d2);
        x=x(m); y=y(m); z=z(m); d2=d2(m);
        D{iA,jS} = struct('x',x,'y',y,'z',z,'d2',d2);
        xlims = [min(xlims(1),min(x)) max(xlims(2),max(x))];
        ylims = [min(ylims(1),min(y)) max(ylims(2),max(y))];
        zlims = [min(zlims(1),min(z)) max(zlims(2),max(z))];
    end
end

% box & contact geometry
cx = mean(xlims); cy = mean(ylims); cz = mean(zlims);
Lx = diff(xlims); Ly = diff(ylims); Lz = diff(zlims);
a_contact = sqrt(max(0, 2*R_tip*h_max - h_max^2));
if isnan(a_contact) || a_contact==0, a_contact = min([Lx,Ly,Lz])/6; end % safe fallback

%% ---------------- STZ METRICS (3D clustering) ----------------
fprintf('\n STZ metrics (3D voxel clustering, 26-connectivity) \n');
for iA = 1:nAlloys
    for jS = 1:nStates
        S = D{iA,jS};
        rr = hypot(S.x - cx, S.y - cy);
        far = rr > far_k * a_contact;
        ref = S.d2(far & isfinite(S.d2));
        if isempty(ref)
            warning('Far-field empty at k=%.1f; relaxing to k=2.0 (Alloy %d, State %s).', ...
                far_k, iA, stateSheets{jS});
            far = rr > 2.0 * a_contact; ref = S.d2(far & isfinite(S.d2));
        end

        if useFixedThreshold
            thr = fixedThr; rule = sprintf('fixed=%.3g',fixedThr);
        elseif use_P95_threshold
            thr = prctile(ref,95); rule = 'P95(far-field)';
        else
            thr = mean(ref) + 3*std(ref); rule = 'μ+3σ (far-field)';
        end

        active = S.d2 >= thr;

        [labels, nComp, sizesVox, volVox, idxLargest] = ...
            labelClusters3D(S.x, S.y, S.z, active, voxel_h, min_vox_cluster);

        stzCount = nComp;
        totVol   = sum(volVox);
        maxVol   = (idxLargest>0) * volVox(max(1,idxLargest));
        eqR      = ((3*totVol)/(4*pi))^(1/3);

        ani = NaN;
        if idxLargest>0
            [cxL,cyL,czL] = voxelCenters(labels==idxLargest, voxel_h, S.x, S.y, S.z);
            if ~isempty(cxL)
                C = cov([cxL,cyL,czL]); ev = sort(eig(C),'descend');
                if numel(ev)==3 && ev(3)>0, ani = ev(1)/ev(3); end
            end
        end

        fprintf('%s | %s: thr(%s)=%.3g  STZ=%d  V_tot=%.2f Å^3  V_max=%.2f Å^3  R_eq=%.2f Å  anis.=%.2f\n', ...
            alloyNames{iA}, stateSheets{jS}, rule, thr, stzCount, totVol, maxVol, eqR, ani);
    end
end

%% ---------------- FIGURES: XZ, XY, YZ ----------------
make_plane_figure('XZ', D, xlims, ylims, zlims, cx, cy, cz, Ly, sliceFrac, gridN, sigma_px, ...
                  usePlasma, prcLo, prcHi, alloyNames, stateSheets);
make_plane_figure('XY', D, xlims, ylims, zlims, cx, cy, cz, Lz, sliceFrac, gridN, sigma_px, ...
                  usePlasma, prcLo, prcHi, alloyNames, stateSheets);
make_plane_figure('YZ', D, xlims, ylims, zlims, cx, cy, cz, Lx, sliceFrac, gridN, sigma_px, ...
                  usePlasma, prcLo, prcHi, alloyNames, stateSheets);
end

%% =============== plane figure helper =================
function make_plane_figure(plane, D, xlims, ylims, zlims, cx, cy, cz, Ldrop, sliceFrac, gridN, sigma_px, ...
                           usePlasma, prcLo, prcHi, alloyNames, stateSheets)

nAlloys = size(D,1);
nStates = size(D,2);

% common grid for this plane
switch upper(plane)
    case 'XZ'
        xi = linspace(xlims(1),xlims(2),gridN);
        zi = linspace(zlims(1),zlims(2),gridN);
        [XI,ZI] = meshgrid(xi, zi);
    case 'XY'
        xi = linspace(xlims(1),xlims(2),gridN);
        zi = linspace(ylims(1),ylims(2),gridN);
        [XI,ZI] = meshgrid(xi, zi);
    case 'YZ'
        xi = linspace(ylims(1),ylims(2),gridN);
        zi = linspace(zlims(1),zlims(2),gridN);
        [XI,ZI] = meshgrid(xi, zi);
    otherwise
        error('plane must be XZ, XY, or YZ');
end

maps = cell(nAlloys,nStates);  allVals = [];

for iA = 1:nAlloys
    for jS = 1:nStates
        S = D{iA,jS};
        switch upper(plane)
            case 'XZ'
                half = 0.5*sliceFrac*Ldrop;
                in = (S.y >= cy - half) & (S.y <= cy + half);
                F = scatteredInterpolant(S.x(in), S.z(in), S.d2(in), 'natural', 'nearest');
                V = F(XI, ZI);
                xlab = 'x (Å)'; ylab = 'z (Å)';
                y0_line = (0 >= zlims(1) && 0 <= zlims(2));
            case 'XY'
                half = 0.5*sliceFrac*Ldrop;
                in = (S.z >= cz - half) & (S.z <= cz + half);
                F = scatteredInterpolant(S.x(in), S.y(in), S.d2(in), 'natural', 'nearest');
                V = F(XI, ZI);
                xlab = 'x (Å)'; ylab = 'y (Å)';
                y0_line = false;
            case 'YZ'
                half = 0.5*sliceFrac*Ldrop;
                in = (S.x >= cx - half) & (S.x <= cx + half);
                F = scatteredInterpolant(S.y(in), S.z(in), S.d2(in), 'natural', 'nearest');
                V = F(XI, ZI);
                xlab = 'y (Å)'; ylab = 'z (Å)';
                y0_line = (0 >= zlims(1) && 0 <= zlims(2));
        end

        if sigma_px > 0, V = gaussBlur2D(V, sigma_px); end
        maps{iA,jS} = V;
        allVals = [allVals; V(:)];
    end
end

% shared color limits for this figure
vmin = prctile(allVals, prcLo);
vmax = prctile(allVals, prcHi);
if ~(isfinite(vmin) && isfinite(vmax)) || vmin==vmax
    vmin = min(allVals); vmax = max(allVals);
end

% colormap
if usePlasma, cmap = sky(256); else, cmap = turbo(256); end

% draw 4×3 panel figure
f = figure('Color','w','Position',[120 80 1320 1180]);
tiledlayout(nAlloys,nStates,'Padding','compact','TileSpacing','compact');

for iA = 1:nAlloys
    for jS = 1:nStates
        nexttile
        imagesc(maps{iA,jS});
        set(gca,'YDir','normal'); axis image
        colormap(cmap); caxis([vmin vmax]); hold on
        if y0_line && strcmpi(plane,'XZ')
            xlim([1 size(maps{iA,jS},2)]);
        end
        title(sprintf('%s  |  %s', alloyNames{iA}, stateSheets{jS}), ...
              'FontWeight','normal')
        xlabel(xlab); ylabel(ylab);
        set(gca,'LineWidth',1,'FontSize',10)

        if iA==1 && jS==nStates
            cb = colorbar; cb.Label.String = 'D^2_{min}'; cb.Box = 'off';
        end
    end
end

outfile = sprintf('D2min_%s_4x3_plasma.png', upper(plane));
exportgraphics(f, outfile, 'Resolution', 300);
fprintf('Saved: %s\n', outfile);
end

%% -------- util & colormaps (unchanged from your version) --------
function name = pickCol(T, prefer)
names = T.Properties.VariableNames; low = lower(names);
p = lower(prefer);
idx = find(strcmp(low,p),1);
if ~isempty(idx), name = names{idx}; return; end
idx = find(contains(low,p),1);
if ~isempty(idx), name = names{idx}; return; end
error('Column "%s" not found. Available: %s', prefer, strjoin(names, ', '));
end

function M = gaussBlur2D(M, sigma)
if sigma <= 0, return; end
sz = max(3, round(6*sigma)); if mod(sz,2)==0, sz=sz+1; end
x = linspace(-((sz-1)/2), ((sz-1)/2), sz);
g = exp(-0.5*(x/sigma).^2); g = g/sum(g);
M = conv2(conv2(M, g, 'same'), g', 'same');
end

function cmap = plasma(n)
if nargin < 1, n = 256; end
anchors = [...
    0.050383, 0.029803, 0.527975;
    0.305868, 0.102646, 0.681348;
    0.542052, 0.220124, 0.750531;
    0.760956, 0.332256, 0.750383;
    0.916242, 0.511643, 0.631478;
    0.977856, 0.744140, 0.444018;
    0.940015, 0.975158, 0.131326];
t0 = linspace(0,1,size(anchors,1));
t  = linspace(0,1,n)';
cmap = [interp1(t0,anchors(:,1),t,'linear'), ...
        interp1(t0,anchors(:,2),t,'linear'), ...
        interp1(t0,anchors(:,3),t,'linear')];
cmap = max(0,min(1,cmap));
end

function [labels, nComp, sizesVox, volVox, idxLargest] = labelClusters3D(x,y,z,active,vox, minSize)
if ~any(active)
    labels = uint32(0); nComp=0; sizesVox=[]; volVox=[]; idxLargest=0; return;
end
XR = [min(x) max(x)]; YR = [min(y) max(y)]; ZR = [min(z) max(z)];
Nx = max(2,ceil((XR(2)-XR(1))/vox));
Ny = max(2,ceil((YR(2)-YR(1))/vox));
Nz = max(2,ceil((ZR(2)-ZR(1))/vox));
ix = min(Nx, max(1, 1+floor((x-XR(1))/vox)));
iy = min(Ny, max(1, 1+floor((y-YR(1))/vox)));
iz = min(Nz, max(1, 1+floor((z-ZR(1))/vox)));
BW = false(Nx,Ny,Nz);
lin = sub2ind([Nx,Ny,Nz], ix(active),iy(active),iz(active));
BW(lin) = true;

labels = uint32(zeros(Nx,Ny,Nz,'uint32'));
lab = uint32(0);
sizesList = [];

neighbors = int32([ ...
 -1 -1 -1; -1 -1 0; -1 -1 1; -1 0 -1; -1 0 0; -1 0 1; -1 1 -1; -1 1 0; -1 1 1; ...
  0 -1 -1;  0 -1 0;  0 -1 1;  0 0 -1;           0 0 1;  0 1 -1;  0 1 0;  0 1 1; ...
  1 -1 -1;  1 -1 0;  1 -1 1;  1 0 -1;  1 0 0;  1 0 1;  1 1 -1;  1 1 0;  1 1 1]);

totalVox = numel(BW);
Q = zeros(totalVox,1,'uint32');

for idx = 1:totalVox
    if ~BW(idx) || labels(idx)~=0, continue; end
    lab = lab + 1;
    qh=1; qt=1; Q(qt) = uint32(idx);
    labels(idx) = lab;
    sz = 1;
    while qh <= qt
        idc = Q(qh); qh = qh+1;
        [ix0,iy0,iz0] = ind2sub([Nx,Ny,Nz], idc);
        for k=1:26
            nx = ix0 + neighbors(k,1);
            ny = iy0 + neighbors(k,2);
            nz = iz0 + neighbors(k,3);
            if nx<1||ny<1||nz<1||nx>Nx||ny>Ny||nz>Nz, continue; end
            if BW(nx,ny,nz) && labels(nx,ny,nz)==0
                qt = qt+1; nid = sub2ind([Nx,Ny,Nz], nx,ny,nz);
                Q(qt) = uint32(nid);
                labels(nx,ny,nz) = lab;
                sz = sz + 1;
            end
        end
    end
    sizesList(double(lab)) = sz;
end

if isempty(sizesList)
    labels(:)=0; nComp=0; sizesVox=[]; volVox=[]; idxLargest=0; return;
end
keep = find(sizesList >= minSize);
if isempty(keep)
    labels(:)=0; nComp=0; sizesVox=[]; volVox=[]; idxLargest=0; return;
end
nComp    = numel(keep);
sizesVox = double(sizesList(keep));
volVox   = sizesVox * (vox^3);

map = zeros(double(lab),1,'uint32'); map(keep) = 1:nComp;
labels = mapLabels(labels, map);
[~,idxLargest] = max(sizesVox);
end

function L2 = mapLabels(L, map)
L2 = uint32(zeros(size(L),'uint32'));
nz = L ~= 0;  L2(nz) = map(double(L(nz)));
end

function [cx,cy,cz] = voxelCenters(mask, vox, x, y, z)
XR = [min(x) max(x)]; YR = [min(y) max(y)]; ZR = [min(z) max(z)];
[Nx,Ny,Nz] = size(mask);
[ix,iy,iz] = ind2sub([Nx,Ny,Nz], find(mask));
cx = XR(1) + (double(ix)-0.5)*vox;
cy = YR(1) + (double(iy)-0.5)*vox;
cz = ZR(1) + (double(iz)-0.5)*vox;
end