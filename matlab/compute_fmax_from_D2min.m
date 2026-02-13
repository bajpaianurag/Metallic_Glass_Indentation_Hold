function metrics = compute_fmax_from_D2min_alloys_t50ps()
% Compute NROI, Npl, Nmax, fmax from D2min snapshots at t = 50 ps.
% One Excel file per alloy; each file has 3 sheets for thermal states.
%
% OUTPUT:
%   - Writes CSV: STZ_connectivity_metrics_t50ps.csv
%   - Returns table 'metrics'
%
% REQUIREMENTS:
%   - MATLAB R2018b+ (for graph/conncomp) and Statistics Toolbox (createns/rangesearch)
%   - Each sheet contains columns for x, y, z, D2min (header names can vary; auto-matched)
%
% PHYSICAL DEFINITIONS (locked):
%   ROI (under-indenter): cylinder centered at (cx,cy)=box XY center
%       r <= rROI, and z in [z_surf - zROIdepth, z_surf]
%   D* threshold: computed from far-field annulus r > far_k * a_contact
%       D* = mu + 3*sigma  (default; can switch to P95 or fixed)
%   Clustering: distance-based on ACTIVE atoms within ROI using cutoff rc (Å)
%       Nmax = size of largest connected component
%       fmax = Nmax / Npl

%% ---------------- USER INPUTS ----------------
% Files (edit paths if needed)
alloyFiles = { ...
    'CuZr.xlsx', ...
    'CuZrAl.xlsx', ...
    'CuZrAlAg.xlsx', ...
    'CuZrAlAgTi.xlsx' ...
};

% Names for reporting
alloyNames = { ...
    'Cu50Zr50', ...
    'Cu45Zr48Al7', ...
    'Cu39Zr46Al7Ag8', ...
    'Cu37Zr42Al7Ag8Ti6' ...
};

% Sheet names (confirmed for CuZr.xlsx: 'as-cast','200','400')
% If any workbook differs, this code will auto-detect alternatives via find_state_sheet().
stateSheetsWanted = {'Sheet1','Sheet2','Sheet3'};
stateNames        = {'As-cast','Annealed 200C','Annealed 400C'};

% Indenter geometry (given)
R_tip = 40;   % Å
h_max = 20;   % Å

% ROI under indenter (recommended)
rROI      = 1.5 * R_tip;  % Å
zROIdepth = 2.0 * R_tip;  % Å

% Threshold D* mode (pick one)
thrMode  = 'mu+3sigma';   % 'mu+3sigma' | 'p95' | 'fixed'
fixedThr = 0.02;          % Å^2 (only used if thrMode='fixed')

% Far-field annulus for D* calibration
a_contact = sqrt(max(0, 2*R_tip*h_max - h_max^2)); % Å
far_k = 3.0;

% Clustering distance cutoff (Å)
% Ideally set to the first minimum of RDF. If unknown, start ~3.6–4.0 Å and refine.
rc = 3.8;

% Output
outCSV = 'STZ_connectivity_metrics_t50ps.csv';

%% ---------------- COMPUTE ----------------
varNames = { ...
    'Alloy','State','File','Sheet', ...
    'R_tip_A','h_max_A','rROI_A','zROIdepth_A','a_contact_A','far_k', ...
    'thrMode','thrRule','Dstar','rc_A', ...
    'NROI','Npl','Nmax','fmax','nClusters' ...
};
nCols = numel(varNames);

nAlloys = numel(alloyFiles);
nStates = numel(stateSheetsWanted);

rows = cell(nAlloys*nStates, nCols);
ridx = 0;

for iA = 1:nAlloys
    file = alloyFiles{iA};

    % list sheets in this workbook
    sh = sheetnames(file);

    % map desired state sheets to actual sheets in this workbook
    sheetMap = cell(1, nStates);
    for jS = 1:nStates
        sheetMap{jS} = find_state_sheet(sh, stateSheetsWanted{jS});
    end

    for jS = 1:nStates
        ridx = ridx + 1;
        sheetName = sheetMap{jS};

        T = readtable(file, 'Sheet', sheetName);

        x  = T.(pickCol(T,'x'));
        y  = T.(pickCol(T,'y'));
        z  = T.(pickCol(T,'z'));
        d2 = T.(pickCol(T,'d2min'));

        m = isfinite(x) & isfinite(y) & isfinite(z) & isfinite(d2);
        x=x(m); y=y(m); z=z(m); d2=d2(m);

        % Indenter axis (assumed at box center in XY)
        cx = 0.5*(min(x)+max(x));
        cy = 0.5*(min(y)+max(y));

        % Free surface (assumed at max(z))
        z_surf = max(z);

        % ROI under indenter
        r = hypot(x-cx, y-cy);
        inROI = (r <= rROI) & (z >= (z_surf - zROIdepth)) & (z <= z_surf);
        NROI = nnz(inROI);

        % Far-field for threshold calibration
        far = (r > far_k * a_contact);
        ref = d2(far); ref = ref(isfinite(ref));
        if isempty(ref)
            far2 = (r > 2.0 * a_contact);
            ref = d2(far2); ref = ref(isfinite(ref));
        end
        if isempty(ref)
            error('Far-field empty for %s | %s. Reduce far_k or check axis.', alloyNames{iA}, stateNames{jS});
        end

        % Threshold D*
        switch lower(thrMode)
            case 'fixed'
                thr = fixedThr;
                thrRule = sprintf('fixed(%.3g)', fixedThr);
            case 'p95'
                thr = prctile(ref,95);
                thrRule = 'P95(far-field)';
            case 'mu+3sigma'
                thr = mean(ref) + 3*std(ref);
                thrRule = 'mu+3sigma(far-field)';
            otherwise
                error('thrMode must be: mu+3sigma | p95 | fixed');
        end

        % Active atoms in ROI
        activeROI = inROI & (d2 >= thr);
        Npl = nnz(activeROI);

        % Cluster active atoms (distance-based)
        if Npl < 2
            Nmax = Npl;
            fmax = double(Npl > 0);
            nClusters = double(Npl > 0);
        else
            Xa = [x(activeROI), y(activeROI), z(activeROI)];
            [Nmax, nClusters] = largest_cluster_size(Xa, rc);
            fmax = Nmax / Npl;
        end

        rows(ridx,:) = { ...
            alloyNames{iA}, stateNames{jS}, file, sheetName, ...
            R_tip, h_max, rROI, zROIdepth, a_contact, far_k, ...
            thrMode, thrRule, thr, rc, ...
            NROI, Npl, Nmax, fmax, nClusters ...
        };
    end
end

rows = rows(1:ridx,:);
metrics = cell2table(rows, 'VariableNames', varNames);
writetable(metrics, outCSV);
fprintf('Saved: %s\n', outCSV);

end

%% ---------------- Helper: robust column picking ----------------
function name = pickCol(T, prefer)
names = T.Properties.VariableNames;
low   = lower(names);
p     = lower(prefer);

% aliases
if strcmp(p,'d2min')
    aliases = {'d2min','d_2min','d2_min','d2','dmin2','d2min_','d2minvalue','d2minval'};
else
    aliases = {p};
end

% exact match first
for k = 1:numel(aliases)
    idx = find(strcmp(low, aliases{k}), 1);
    if ~isempty(idx)
        name = names{idx}; return;
    end
end

% contains match
idx = find(contains(low, p), 1);
if ~isempty(idx)
    name = names{idx}; return;
end

error('Column "%s" not found. Available: %s', prefer, strjoin(names, ', '));
end

%% ---------------- Helper: map desired sheet to actual sheet ----------------
function sheetOut = find_state_sheet(sheetList, desired)
% Finds a sheet name in sheetList matching desired, robust to case and minor formatting.
% Examples: desired='as-cast' may match 'as-cast' or 'AsCast' or 'as cast'.
desired0 = lower(strrep(strrep(desired,'_',''),' ',''));
S = string(sheetList);
S0 = lower(strrep(strrep(S,'_',''),' ',''));

% direct match
idx = find(S0 == desired0, 1);
if ~isempty(idx)
    sheetOut = char(S(idx)); return;
end

% if desired is numeric like '200' or '400', allow match within sheet name
if all(isstrprop(desired,'digit'))
    idx = find(contains(S0, desired0), 1);
    if ~isempty(idx)
        sheetOut = char(S(idx)); return;
    end
end

% fallback: contains match
idx = find(contains(S0, desired0), 1);
if ~isempty(idx)
    sheetOut = char(S(idx)); return;
end

error('Could not find a sheet matching "%s". Available sheets: %s', desired, strjoin(cellstr(S), ', '));
end

%% ---------------- Helper: largest connected cluster size ----------------
function [Nmax, nClusters] = largest_cluster_size(X, rc)
% X: N×3 coordinates of active atoms
% rc: distance cutoff (Å)
N = size(X,1);

Mdl = createns(X, 'NSMethod','kdtree');
idx = rangesearch(Mdl, X, rc);

ii = []; jj = [];
for i = 1:N
    nbr = idx{i};
    nbr(nbr==i) = [];
    if ~isempty(nbr)
        ii = [ii; repmat(i, numel(nbr), 1)];
        jj = [jj; nbr(:)];
    end
end

A = sparse([ii; jj], [jj; ii], 1, N, N);  % symmetrized adjacency
G = graph(A);

bins = conncomp(G);
nClusters = max(bins);

counts = accumarray(bins(:), 1);
Nmax = max(counts);
end