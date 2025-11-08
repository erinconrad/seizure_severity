%% ============================================================
%  Spearman correlations + matching figure style (2x2 panels)
%  Panels:
%    A) All epilepsy (overall) — r, p (unadjusted)
%    B) Frontal — r, p_bonf    C) Temporal — r, p_bonf    D) General — r, p_bonf
%  Style: log–log scatter, dotted zero guides, zero-origin stars (same color)
%  ------------------------------------------------------------
%clear; clc;

only_amb = 0; % 1 = ONLY ALLOW AMBULATORIES; 2 = No ambulatories; 0 = everything

%% ===== Paths (EDIT) =====
outCsv     = '../data/SN_counts/spike_counts_summary.csv';
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-10-20_1418.csv';
outFigure  = '../figures/spearman_spikerate_szfreq_panels.png';
resultsCsv = '../output/spearman_spikerate_szfreq_results.csv';   % [] to skip saving

%% ===== Load spike counts; per-patient mean spike rate =====
S = readtable(outCsv,'TextType','string');
S.SpikeRate_Hz      = S.Total_Spikes ./ S.Duration_sec;
S.SpikeRate_perHour = S.SpikeRate_Hz * 3600;

if only_amb == 1
    S(S.Duration_sec < 3600*12,:) = [];
elseif only_amb == 2
    S(S.Duration_sec > 3600*12,:) = [];
end

Sg = groupsummary(S,'Patient','mean','SpikeRate_perHour');
Sg.Properties.VariableNames = {'Patient','GroupCount','MeanSpikeRate_perHour'};
Sg.GroupCount = [];

%% ===== Load report; keep epilepsy only =====
R = readtable(reportFile,'TextType','string');
if ~isnumeric(R.patient_id), R.patient_id = double(str2double(R.patient_id)); end
if ~isstring(R.epilepsy_type), R.epilepsy_type = string(R.epilepsy_type); end
assert(ismember('sz_freqs',R.Properties.VariableNames),'Column "sz_freqs" not found.');

badTypes = lower(["Non-Epileptic Seizure Disorder","Uncertain if Epilepsy","Unknown or MRN not found",""]);
isEmpty  = ismissing(R.epilepsy_type) | strlength(strtrim(R.epilepsy_type))==0;
Rt = sortrows(R(~isEmpty, {'patient_id','epilepsy_type'}), 'patient_id');
[uniq_pid, ia] = unique(Rt.patient_id, 'stable');
PerPatType = table(uniq_pid, Rt.epilepsy_type(ia), 'VariableNames', {'Patient','EpilepsyType'});

epi_norm = lower(strtrim(PerPatType.EpilepsyType));
isBad    = ismember(epi_norm, badTypes);
validPatients = PerPatType.Patient(~isBad);
R = R(ismember(R.patient_id, validPatients), :);

fprintf('Total patients in report: %d\n', numel(unique(R.patient_id(~isnan(R.patient_id)))));
fprintf('Patients with ANY epilepsy type kept: %d\n', numel(unique(R.patient_id)));
fprintf('Excluded by epilepsy_type: %d\n\n', sum(isBad));

%% ===== Per-patient seizure frequency (robust parsing + rescue) =====
pids = unique(R.patient_id(~isnan(R.patient_id)));
MeanSzFreq = nan(numel(pids),1);

for k = 1:numel(pids)
    pid = pids(k);
    rr  = R(R.patient_id==pid,:);
    vals = [];

    % Parse sz_freqs arrays (tolerate "null")
    for j = 1:height(rr)
        raw = strtrim(rr.sz_freqs(j));
        if strlength(raw)==0 || raw=="[]" || raw==""; continue; end
        s = regexprep(raw, 'null', 'NaN', 'ignorecase');
        try
            v = jsondecode(char(s)); v = v(:);
        catch
            nums = regexp(s,'[-+]?\d+(\.\d+)?([eE][-+]?\d+)?','match');
            v = []; if ~isempty(nums), v = str2double(string(nums(:))); end
        end
        if ~isempty(v)
            v(~isfinite(v)) = NaN;
            v(v < 0)        = NaN;   % negatives treated as missing placeholders
            vals = [vals; v]; %#ok<AGROW>
        end
    end

    % Rescue: if all missing & visit_hasSz exists & all zeros -> 0
    if (isempty(vals) || all(isnan(vals))) && ismember('visit_hasSz', rr.Properties.VariableNames)
        rawHas = string(rr.visit_hasSz);
        allHas = [];
        for jj = 1:numel(rawHas)
            sHas = strtrim(rawHas(jj));
            if strlength(sHas)==0 || sHas=="[]" || sHas==""; continue; end
            try
                vHas = jsondecode(char(sHas)); vHas = vHas(:);
            catch
                nums = regexp(sHas,'\d+','match'); vHas = [];
                if ~isempty(nums), vHas = str2double(string(nums(:))); end
            end
            allHas = [allHas; vHas]; %#ok<AGROW>
        end
        nonnanHas = allHas(~isnan(allHas));
        if ~isempty(nonnanHas) && all(nonnanHas==0)
            vals = 0;
        end
    end

    MeanSzFreq(k) = mean(vals,'omitnan');
end

Rg = table(pids, MeanSzFreq, 'VariableNames', {'Patient','MeanSzFreq'});

%% ===== Attach EpilepsySpecific and backfill "General" =====
if ismember('epilepsy_specific', R.Properties.VariableNames)
    if ~isstring(R.epilepsy_specific), R.epilepsy_specific = string(R.epilepsy_specific); end
    spec_ok = ~ismissing(R.epilepsy_specific) & strlength(strtrim(R.epilepsy_specific))>0;
    RtSpec = sortrows(R(spec_ok, {'patient_id','epilepsy_specific'}), 'patient_id');
    [uniq_pid_s, ia_s] = unique(RtSpec.patient_id, 'stable');
    PerPatSpecific = table(uniq_pid_s, RtSpec.epilepsy_specific(ia_s), ...
        'VariableNames', {'Patient','EpilepsySpecific'});

    Rg = outerjoin(Rg, PerPatSpecific, 'Keys','Patient', 'MergeKeys', true);

    % Coalesce potential duplicates from outerjoin
    specMask = startsWith(Rg.Properties.VariableNames, 'EpilepsySpecific','IgnoreCase',true);
    specCols = Rg.Properties.VariableNames(specMask);
    if ~isempty(specCols)
        baseSpec = 'EpilepsySpecific';
        if ~any(strcmp(specCols, baseSpec)), Rg = renamevars(Rg, specCols{1}, baseSpec); end
        others = setdiff(specCols, baseSpec, 'stable');
        if ~isstring(Rg.(baseSpec)), Rg.(baseSpec) = string(Rg.(baseSpec)); end
        Rg.(baseSpec) = strtrim(Rg.(baseSpec));
        for ii = 1:numel(others)
            v = Rg.(others{ii}); if ~isstring(v), v = string(v); end; v = strtrim(v);
            fillIdx = (ismissing(Rg.(baseSpec)) | strlength(Rg.(baseSpec))==0) & ...
                      (~ismissing(v) & strlength(v)>0);
            Rg.(baseSpec)(fillIdx) = v(fillIdx);
        end
        if ~isempty(others), Rg = removevars(Rg, others); end
    else
        Rg.EpilepsySpecific = strings(height(Rg),1);
    end
else
    Rg.EpilepsySpecific = strings(height(Rg),1);
end

% Backfill "General" from high-level generalized type
Rg = innerjoin(Rg, PerPatType(~isBad,:), 'Keys','Patient');
spec_norm = lower(strtrim(Rg.EpilepsySpecific));
type_norm = lower(strtrim(Rg.EpilepsyType));
isTemporal    = contains(spec_norm,"temporal");
isFrontal     = contains(spec_norm,"frontal");
isGeneralType = contains(type_norm,"general");

SpecCanon = strings(size(spec_norm));
SpecCanon(isTemporal) = "Temporal Lobe";
SpecCanon(isFrontal)  = "Frontal Lobe";
SpecCanon((strlength(SpecCanon)==0) & isGeneralType) = "General";
Rg.EpilepsySpecific = SpecCanon;

%% ===== Merge with per-patient spike rate; keep complete cases =====
P_all = innerjoin(Sg, Rg, 'Keys','Patient');  % all epilepsy (valid types)
P_all = P_all(isfinite(P_all.MeanSzFreq) & isfinite(P_all.MeanSpikeRate_perHour), :);

keep3 = ismember(P_all.EpilepsySpecific, ["General","Temporal Lobe","Frontal Lobe"]);
P     = P_all(keep3, :);

fprintf('Patients with both measures (all epilepsy): %d\n', height(P_all));
fprintf('Patients in modeled 3 groups: %d\n\n', height(P));

%% ===== Spearman correlations (raw scale) =====
% Overall
x_all = P_all.MeanSpikeRate_perHour;
y_all = P_all.MeanSzFreq;
mask_all = isfinite(x_all) & isfinite(y_all);
[rs_all, p_all] = corr(x_all(mask_all), y_all(mask_all), 'Type','Spearman','Rows','complete');
n_all = sum(mask_all);

% Subgroups
groupsWanted = ["Frontal Lobe","Temporal Lobe","General"];
rows = {};
for g = groupsWanted
    m = (P.EpilepsySpecific == g);
    x = P.MeanSpikeRate_perHour(m);
    y = P.MeanSzFreq(m);
    mask = isfinite(x) & isfinite(y);
    n = sum(mask);
    if n >= 3
        [rs, p] = corr(x(mask), y(mask), 'Type','Spearman','Rows','complete');
    else
        rs = NaN; p = NaN;
    end
    rows(end+1,:) = {char(g), n, rs, p}; %#ok<SAGROW>
end

Results = cell2table(rows, 'VariableNames', {'Group','N','Spearman_r','p_raw'});
k = height(Results);               % should be 3
Results.p_bonf = min(Results.p_raw * k, 1);

% Print results
disp('=== Spearman correlations: SpikeRate_perHour vs MeanSzFreq ===');
disp([table("Overall (all epilepsy)", n_all, rs_all, p_all, NaN, ...
      'VariableNames', {'Group','N','Spearman_r','p_raw','p_bonf'}); Results])

%% ===== Plotting prep (log display; Spearman stats unchanged by log) =====
% Epsilons for zero-friendly logs
minpos_rate = min(P_all.MeanSpikeRate_perHour(P_all.MeanSpikeRate_perHour>0)); if isempty(minpos_rate), minpos_rate=1e-6; end
minpos_sz   = min(P_all.MeanSzFreq(P_all.MeanSzFreq>0));                       if isempty(minpos_sz),   minpos_sz=1e-6; end
eps_rate = 0.5*minpos_rate;
eps_sz   = 0.5*minpos_sz;

% Log-projected coordinates for display
Tall = table;
Tall.logSpikeRate = log(P_all.MeanSpikeRate_perHour + (P_all.MeanSpikeRate_perHour<=0).*eps_rate);
Tall.logSzFreq    = log(P_all.MeanSzFreq           + (P_all.MeanSzFreq<=0).*eps_sz);

T = table;
T.EpilepsySpecific = categorical(P.EpilepsySpecific);
T.EpilepsySpecific = reordercats(removecats(T.EpilepsySpecific), {'Frontal Lobe','Temporal Lobe','General'});
T.logSpikeRate = log(P.MeanSpikeRate_perHour + (P.MeanSpikeRate_perHour<=0).*eps_rate);
T.logSzFreq    = log(P.MeanSzFreq           + (P.MeanSzFreq<=0).*eps_sz);
ok = isfinite(T.logSpikeRate) & isfinite(T.logSzFreq) & ~ismissing(T.EpilepsySpecific);
T = T(ok,:);
present = categories(T.EpilepsySpecific);

% For zero-origin detection on modeled subset
isZeroSz_T   = (P.MeanSzFreq==0);                isZeroSz_T   = isZeroSz_T(ok);
isZeroRate_T = (P.MeanSpikeRate_perHour==0);     isZeroRate_T = isZeroRate_T(ok);

% Zero masks for overall
isZeroSz_all   = (P_all.MeanSzFreq==0);
isZeroRate_all = (P_all.MeanSpikeRate_perHour==0);
onlySz_all     =  isZeroSz_all & ~isZeroRate_all;
onlyRate_all   = ~isZeroSz_all &  isZeroRate_all;
bothZero_all   =  isZeroSz_all &  isZeroRate_all;
nonZero_all    = ~(isZeroSz_all | isZeroRate_all);

% Guide locations in log-space
xZero = log(eps_sz);  % seizure=0 placed here
yZero = log(eps_rate);% spike=0 placed here

%% ===== Figure (2x2) in prior style, but Spearman stats in annotations =====
f = figure('Color','w','Position',[60 60 1200 900]);
tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
pal   = lines(3);
fontL = 20;

% Tiny nudges so stars don’t sit under the guide lines
xrange = range(Tall.logSzFreq);    if xrange==0, xrange=1; end
yrange = range(Tall.logSpikeRate); if yrange==0, yrange=1; end
dx = 0.008*xrange; 
dy = 0.008*yrange;

% -------- Panel A: Overall (Spearman r, p) --------
axA = nexttile(1); hold(axA,'on'); grid(axA,'on'); box(axA,'off');
xline(axA, xZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
yline(axA, yZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

baseColor = [0.4 0.4 0.4];
scatter(axA, Tall.logSzFreq(nonZero_all), Tall.logSpikeRate(nonZero_all), ...
        14, baseColor, 'filled', 'MarkerFaceAlpha', 0.25);
plot(axA, Tall.logSzFreq(onlySz_all)+dx,   Tall.logSpikeRate(onlySz_all),   '*','Color',baseColor,'MarkerSize',7,'LineWidth',1);
plot(axA, Tall.logSzFreq(onlyRate_all),    Tall.logSpikeRate(onlyRate_all)+dy,'*','Color',baseColor,'MarkerSize',7,'LineWidth',1);
plot(axA, Tall.logSzFreq(bothZero_all)+dx, Tall.logSpikeRate(bothZero_all)+dy,'*','Color',baseColor,'MarkerSize',8,'LineWidth',1.2);

% Visual aid: OLS line (stats are Spearman-only)
X = [ones(sum(nonZero_all),1), Tall.logSzFreq(nonZero_all)];
b = X \ Tall.logSpikeRate(nonZero_all);
xgrid = linspace(min(Tall.logSzFreq), max(Tall.logSzFreq), 300)';
plot(axA, xgrid, b(1) + b(2)*xgrid, 'k-', 'LineWidth', 2);

% Place zero-guide labels centered after limits finalize
xl = xlim(axA); yl = ylim(axA);
text(axA, xZero + 0.03*(xl(2)-xl(1)), mean(yl), 'Seizure frequency = 0', ...
     'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
     'FontSize',fontL-2, 'Color',[0.25 0.25 0.25]);
text(axA, mean(xl), yZero + 0.06*(yl(2)-yl(1)), 'Spike rate = 0', ...
     'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
     'FontSize',fontL-2, 'Color',[0.25 0.25 0.25]);

xlabel(axA,'log Seizure Frequency','FontSize',fontL);
ylabel(axA,'log Spike Rate','FontSize',fontL);
title(axA, sprintf('A. All epilepsy  (N=%d)', height(Tall)), 'FontSize',fontL, 'FontWeight','bold');

txtA = sprintf('Spearman r=%.3f, p=%.3g', rs_all, p_all);
text(axA, 0.98, 0.95, txtA, 'Units','normalized', ...
     'HorizontalAlignment','right', 'VerticalAlignment','top', ...
     'FontSize',fontL-2, 'FontWeight','bold');
set(axA,'FontSize',fontL);

% -------- Panels B/C/D: Frontal / Temporal / General --------
order      = {'Frontal Lobe','Temporal Lobe','General'};
panelTitle = {'B. Frontal','C. Temporal','D. General'};

for p = 1:3
    ax = nexttile(p+1); hold(ax,'on'); grid(ax,'on'); box(ax,'off');

    if ~ismember(order{p}, present)
        axis(ax,'off'); continue;
    end

    g   = order{p};
    gi  = find(strcmp(present, g), 1);
    col = pal(gi,:);
    m   = (T.EpilepsySpecific == g);

    % zero masks (modeled subset)
    onlySz_g   =  m &  isZeroSz_T & ~isZeroRate_T;
    onlyRate_g =  m & ~isZeroSz_T &  isZeroRate_T;
    bothZero_g =  m &  isZeroSz_T &  isZeroRate_T;
    nonZero_g  =  m & ~(isZeroSz_T | isZeroRate_T);

    xline(ax, xZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
    yline(ax, yZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

    scatter(ax, T.logSzFreq(nonZero_g), T.logSpikeRate(nonZero_g), 18, col, 'filled', 'MarkerFaceAlpha', 0.35);
    if any(onlySz_g),   plot(ax, T.logSzFreq(onlySz_g)+dx,   T.logSpikeRate(onlySz_g),   '*','Color',col,'MarkerSize',8,'LineWidth',1.1); end
    if any(onlyRate_g), plot(ax, T.logSzFreq(onlyRate_g),    T.logSpikeRate(onlyRate_g)+dy,'*','Color',col,'MarkerSize',8,'LineWidth',1.1); end
    if any(bothZero_g), plot(ax, T.logSzFreq(bothZero_g)+dx, T.logSpikeRate(bothZero_g)+dy,'*','Color',col,'MarkerSize',9,'LineWidth',1.2); end

    % Visual OLS trend for this group
    if nnz(nonZero_g) >= 3
        Xg = [ones(nnz(nonZero_g),1), T.logSzFreq(nonZero_g)];
        bg = Xg \ T.logSpikeRate(nonZero_g);
        xg = linspace(min(T.logSzFreq(nonZero_g)), max(T.logSzFreq(nonZero_g)), 250)';
        plot(ax, xg, bg(1)+bg(2)*xg, '-', 'Color', col, 'LineWidth', 2);
    end

    % Centered zero-guide labels after limits finalize
    xl = xlim(ax); yl = ylim(ax);
    text(ax, xZero + 0.03*(xl(2)-xl(1)), mean(yl), 'Seizure frequency = 0', ...
         'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
         'FontSize',fontL-2, 'Color',[0.25 0.25 0.25]);
    text(ax, mean(xl), yZero + 0.06*(yl(2)-yl(1)), 'Spike rate = 0', ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
         'FontSize',fontL-2, 'Color',[0.25 0.25 0.25]);

    xlabel(ax,'log Seizure Frequency','FontSize',fontL);
    ylabel(ax,'log Spike Rate','FontSize',fontL);

    % Pull subgroup stats from Results (Spearman; Bonferroni-corrected p)
    row = Results(strcmp(Results.Group, string(g)), :);
    if ~isempty(row) && row.N >= 3 && isfinite(row.Spearman_r)
        txt = sprintf('Spearman r=%.3f, p_{bonf}=%.3g', row.Spearman_r, row.p_bonf);
    else
        txt = 'Insufficient data';
    end

    nNow = sum(T.EpilepsySpecific==g);
    title(ax, sprintf('%s (N=%d)', panelTitle{p}, nNow), 'FontSize', fontL, 'FontWeight','bold');

    % Right-aligned annotation
    text(ax, 0.98, 0.95, txt, 'Units','normalized', ...
         'HorizontalAlignment','right', 'VerticalAlignment','top', ...
         'FontSize', fontL-3, 'FontWeight','bold');

    set(ax,'FontSize',fontL);
end

%% ===== Save figure =====
if ~isempty(outFigure)
    if ~exist(fileparts(outFigure),'dir'), mkdir(fileparts(outFigure)); end
    exportgraphics(f, outFigure, 'Resolution', 300);
    fprintf('Saved figure: %s\n', outFigure);
end

%% ===== Save results (optional) =====
if ~isempty(resultsCsv)
    if ~exist(fileparts(resultsCsv),'dir'), mkdir(fileparts(resultsCsv)); end
    % Add overall row to table for saving
    ResultsOut = [table("Overall (all epilepsy)", n_all, rs_all, p_all, NaN, ...
        'VariableNames', {'Group','N','Spearman_r','p_raw','p_bonf'}); Results];
    writetable(ResultsOut, resultsCsv);
    fprintf('Saved results to: %s\n', resultsCsv);
end
